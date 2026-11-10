"""Exact dimensionless pressure/NEC identities and falsification controls."""
import sympy as s

def run():
 z=s.symbols('z',nonnegative=True); h0=s.symbols('H0',positive=True); om=s.symbols('Omega_m'); H=s.Function('H')(z); D=s.Function('D')(z)
 p=(s.Rational(2,3)*(1+z)*H*s.diff(H,z)-H**2)/h0**2
 pd=s.simplify(p.subs(H,h0/s.diff(D,z)).doit())
 expected=-s.Rational(2,3)*(1+z)*s.diff(D,z,2)/s.diff(D,z)**3-1/s.diff(D,z)**2
 nec=s.simplify(pd+1/s.diff(D,z)**2)
 lcd=s.simplify(p.subs(H,h0*s.sqrt(om*(1+z)**3+1-om)).doit())
 phantom=s.simplify(nec.subs(D,s.log(1+z)).doit()) # positive NEC control
 negative=s.simplify(nec.subs(D,(1+z)**2-1).doit())
 stiff=s.simplify(p.subs(H,h0*(1+z)**3).doit())
 checks={'distance_pressure_identity':s.simplify(pd-expected)==0,'LCDM_pressure':lcd==om-1,'matter_pressure':lcd.subs(om,1)==0,'deSitter_pressure':lcd.subs(om,0)==-1,'NEC_identity':s.simplify(nec+s.Rational(2,3)*(1+z)*s.diff(D,z,2)/s.diff(D,z)**3)==0,'positive_NEC_control':phantom.subs(z,1)>0,'negative_NEC_control':negative.subs(z,1)<0,'positive_pressure_control':stiff.subs(z,1)>0}
 return {'claim':'A1','checks':{k:bool(v) for k,v in checks.items()},'passed':all(checks.values()),'pressure':str(pd),'NEC':str(nec),'LCDM_pressure':str(lcd)}
