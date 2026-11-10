"""Equation-based conditional reconstruction; see METHODS.md for approximations."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['VECLIB_MAXIMUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,itertools,math,time
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.integrate import cumulative_trapezoid
from scipy.linalg import solve_triangular,qr
from validation import checked_cholesky, diagnostics, profile_diagnostics
from scipy.optimize import least_squares
from scipy.stats import multivariate_t
from scipy.special import logsumexp
R=Path(__file__).resolve().parents[1];O=Path(os.environ['SD_OUTPUT']);D=Path(os.environ['SD_DATA']);C=299792.458
knots=np.array([0,.35,.512,.782,1.236,1.626,2.33]);zg=np.unique(np.r_[np.linspace(0,2.33,467),.47]);zi=np.linspace(0,2.4,48001)
oi=.3153; integ=cumulative_trapezoid(1/np.sqrt(oi*(1+zi)**3+1-oi),zi,initial=0)
def baseline(z):
 e=np.sqrt(oi*(1+z)**3+1-oi);return np.interp(z,zi,integ),1/e,-1.5*oi*(1+z)**2/e**3
rows=np.genfromtxt(D/'desi-dr2/desi_gaussian_bao_ALL_GCcomb_mean.txt',dtype=None,encoding='utf8');zb=np.array([x[0] for x in rows]);yb=np.array([x[1] for x in rows]);typ=np.array([x[2] for x in rows]);cb=np.loadtxt(D/'desi-dr2/desi_gaussian_bao_ALL_GCcomb_cov.txt');wb=solve_triangular(checked_cholesky(cb),np.eye(len(zb)),lower=True)
sn=np.array([[float(v) for v in l.split()[3:7]] for l in (D/'des-dovekie/DES-Dovekie_HD.csv').read_text().splitlines() if l.startswith('SN:')]);zs,zh,mu,mue=sn.T
class Model:
 def __init__(self,closure='I4',stat=False,fixed_mag=False,zhel=True,fixed_h=None,pressure=True):
  self.pressure=pressure
  bc={'I4':([(1,np.zeros(6)),(3,np.zeros(6))],[(2,np.zeros(6))]),'I2':([(1,np.zeros(6)),(2,np.zeros(6))],[(2,np.zeros(6))])}
  kw=dict(k=4)
  if closure=='minimal':kw.update(t=np.r_[np.repeat(0.,5),knots[2:-2],np.repeat(2.33,5)],bc_type=([(1,np.zeros(6))],None))
  else:kw['bc_type']=bc[closure]
  self.spl=make_interp_spline(knots,np.vstack([np.zeros(6),np.eye(6)]),**kw);self.fixed_h=fixed_h;self.fixed_mag=fixed_mag
  self.bs=self.spl(zs)*.01;self.bb=self.spl(zb)*.01;self.bd=self.spl(zb,1)*.01;self.gb=self.spl(zg,1)*.01;self.g2=self.spl(zg,2)*.01;self.bg=baseline(zg)
  self.snbase=5*np.log10(C/70*baseline(zs)[0]*(1+(zh if zhel else zs)))+25
  f=np.load(D/'des-dovekie'/('STATONLY.npz' if stat else 'STAT+SYS.npz'));n=int(f['nsn'][0]);P=np.zeros((n,n));P[np.triu_indices(n)]=f['cov'];P=P+P.T-np.diag(P.diagonal());self.L=checked_cholesky(P)
  self.one=self.L.T@np.ones(n);self.unit=self.one/np.linalg.norm(self.one)
  def white(x):
   w=self.L.T@x
   return w if fixed_mag else w-self.unit[:,None]*(self.unit@w)[None,:]
  self.white=white;self.obs=white((mu-self.snbase)[:,None])[:,0]
  self.u=self.bs/baseline(zs)[0][:,None];self.powers=[];cols=[]
  for degree in (1,2,3):
   for inds in itertools.combinations_with_replacement(range(6),degree):
    coeff=math.factorial(degree)
    for v in np.bincount(inds,minlength=6):coeff/=math.factorial(int(v))
    cols.append(5/np.log(10)*(-1)**(degree+1)/degree*coeff*np.prod(self.u[:,inds],axis=1));self.powers.append(inds)
  F=white(np.array(cols).T)
  if fixed_mag:F=np.column_stack([F,self.one])
  Q,self.J=qr(F,mode='economic');self.target=Q.T@self.obs;self.constant=self.obs@self.obs-self.target@self.target
  self.free=8 if fixed_mag else 7
 def features(self,x):
  x=np.atleast_2d(x);p=np.array([np.prod(x[:,inds],axis=1) for inds in self.powers]).T
  if self.fixed_mag:p=np.column_stack([p,-5*np.log10((100+x[:,6])/x[:,7]/.7)])
  return p
 def residuals(self,x):
  x=np.atleast_2d(x);dm=baseline(zb)[0]+x[:,:6]@self.bb.T;dh=baseline(zb)[1]+x[:,:6]@self.bd.T
  pred=dm.copy();pred[:,typ=='DH_over_rs']=dh[:,typ=='DH_over_rs'];iv=typ=='DV_over_rs';pred[:,iv]=np.cbrt(zb[iv]*dm[:,iv]**2*dh[:,iv]);pred*=C/(100*(100+x[:,6,None]));b=(pred-yb)@wb.T
  snres=self.features(x)@self.J.T-self.target;out=np.column_stack([b,snres])
  if self.fixed_mag:out=np.column_stack([out,(x[:,7]-147.09)/.27])
  if self.fixed_h is not None:out=np.column_stack([out,((100+x[:,6])/self.fixed_h-147.09)/.27])
  return out
 def physical(self,x):
  d=self.bg[1]+x[:,:6]@self.gb.T;dd=self.bg[2]+x[:,:6]@self.g2.T
  return (abs(x[:,:6])<50).all(1)&(d>0).all(1)&((d+2/3*(1+zg)*dd>0).all(1) if self.pressure else True)&((100+x[:,6])>73.5)&((100+x[:,6])<176.5)
 def integrate(self,name,n=50000,seed=731009):
  rng=np.random.default_rng(seed);start=np.r_[np.zeros(6),-.3,147.09] if self.fixed_mag else np.r_[np.zeros(6),-.3]
  fit=least_squares(lambda p:self.residuals(p)[0],start,jac='3-point',xtol=1e-11,ftol=1e-11,gtol=1e-9);cv=np.linalg.inv(fit.jac.T@fit.jac);rv=multivariate_t(fit.x,cv,df=8,seed=rng);samples=rv.rvs(n);lw=[]
  for b in np.array_split(samples,max(1,n//500)):
   r=self.residuals(b);v=-.5*np.sum(r*r,axis=1)-rv.logpdf(b);v[~self.physical(b)]=-np.inf
   if self.fixed_mag:v-=np.log(b[:,7]) # flat h and Gaussian rd in q,rd coordinates
   lw.extend(v)
  lw=np.array(lw);w=np.exp(lw-logsumexp(lw));ess=1/np.sum(w*w);maxw=w.max()
  if self.fixed_h is not None:h=np.full(n,self.fixed_h)
  else:
   rd=samples[:,7] if self.fixed_mag else rng.normal(147.09-.27**2/147.09,.27,n);h=(100+samples[:,6])/rd
  q=100+samples[:,6];d1=baseline(np.array([.47]))[1][0]+samples[:,:6]@(self.spl(.47,1)*.01);H=100*h/d1/1.47
  def moments(a):
   ave=w@a;return [float(ave),float(np.sqrt(w@((a-ave)**2)))]
  def quant(a):
   ix=np.argsort(a);return np.interp([.16,.84],np.cumsum(w[ix]),a[ix])
  mean,sd=moments(H);qq=quant(H);ix=np.argsort(H);cw=np.cumsum(w[ix]);left=np.linspace(0,.32,1000);width=np.interp(left+.68,cw,H[ix])-np.interp(left,cw,H[ix]);hp=[float(width.min()/2)]
  idx=rng.choice(n,40,p=w);diff=[]
  for p in samples[idx]:
   theory=5*np.log10(1+self.u@p[:6]);
   if self.fixed_mag:theory-=5*np.log10((100+p[6])/p[7]/.7)
   wr=self.obs-self.white(theory[:,None])[:,0];ap=self.features(p)@self.J.T-self.target;diff.append(float(wr@wr-self.constant-np.sum(ap*ap)))
  out=dict(name=name,seed=seed,n=n,ESS=float(ess),max_weight=float(maxw),H047_mean=mean,H047_sd=sd,half68=float((qq[1]-qq[0])/2),shortest68_half=hp[0],H0=moments(100*h),rd_fixed=moments(100*q/147.09/d1/1.47),h_plugin=moments(100*(w@h)/d1/1.47),corr_h_H=float(w@((h-w@h)*(H-mean))/np.sqrt(w@(h-w@h)**2)/sd),exact_delta_chi2_max=float(np.max(np.abs(diff))),map=fit.x.tolist(),map_cov=cv.tolist(),chi2=float(np.sum(fit.fun**2)+self.constant))
  # Direct weighted pointwise percentiles; no redshift averaging or resampling.
  upper=[]; median=[]; lower=[]
  for z in zg:
   d=baseline(np.array([z]))[1][0]+samples[:,:6]@(self.spl(z,1)*.01)
   dd=baseline(np.array([z]))[2][0]+samples[:,:6]@(self.spl(z,2)*.01)
   nec=-2/3*(1+z)*dd/d**3
   ix=np.argsort(nec); qn=np.interp([.025,.5,.975],np.cumsum(w[ix]),nec[ix])
   lower.append(float(qn[0]));median.append(float(qn[1]));upper.append(float(qn[2]))
  out.update(closure=name.split('_')[0],pressure_prior=self.pressure,NEC_min_upper975=min(upper),NEC_no_negative_upper975=bool(min(upper)>=0),NEC_grid_points=len(zg))
  profile=dict(z=zg.tolist(),lower025=lower,median=median,upper975=upper)
  if not diagnostics(out) or not profile_diagnostics(out,profile):raise ValueError('Invalid numerical result/profile: '+name)
  (O/(name+'_nec.json')).write_text(json.dumps(profile,indent=2,allow_nan=False)+'\n')
  (O/(name+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');return out
