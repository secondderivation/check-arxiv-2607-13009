import hashlib, math
import numpy as np
from scipy.linalg import cholesky
from pathlib import Path

RUNS=[f'{c}_{p}_{s}' for c in ['I4','I2'] for p in ['RA','REIRA'] for s in [731009,20261009]]

def digest(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def verify_manifest(root):
 root=Path(root); entries={}
 for line in (root/'MANIFEST.txt').read_text().splitlines():
  h,name=line.split('  ',1)
  p=Path(name)
  if p.is_absolute() or '..' in p.parts or name in entries: raise ValueError('Unsafe or duplicate manifest entry')
  if (root/name).is_symlink() or digest(root/name)!=h: raise ValueError('Hash mismatch: '+name)
  entries[name]=h
 actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and not any(x in ('.git','.venv','_runs','_verify','__pycache__') for x in p.relative_to(root).parts)}
 if actual != set(entries)|{'MANIFEST.txt'}: raise ValueError('Unexpected or missing package files: '+str(actual.symmetric_difference(set(entries)|{'MANIFEST.txt'})))
 return entries

def diagnostics(row):
 keys=['ESS','max_weight','exact_delta_chi2_max','H047_mean','H047_sd','NEC_min_upper975']
 return finite_tree(row) and all(k in row for k in keys) and all(math.isfinite(row[k]) for k in keys) and row['ESS']>=1000 and 0<row['max_weight']<=.01 and 0<=row['exact_delta_chi2_max']<=.01 and row['NEC_min_upper975']>=0 and row['NEC_no_negative_upper975'] and row['H047_sd']>0

def compare(got,ref):
 if not diagnostics(got):return False
 limits={'H047_mean':.002,'H047_sd':.002,'NEC_min_upper975':.02,'ESS':.05*ref['ESS'],'max_weight':.10*ref['max_weight']}
 return all(abs(got[k]-ref[k])<=limit for k,limit in limits.items())

def seed_agreement(rows):
 return all(abs(rows[f'{c}_{p}_731009'][k]-rows[f'{c}_{p}_20261009'][k])<=.03 for c in ['I4','I2'] for p in ['RA','REIRA'] for k in ['H047_mean','H047_sd'])


def finite_tree(value):
 if isinstance(value,dict):return all(finite_tree(v) for v in value.values())
 if isinstance(value,(list,tuple)):return all(finite_tree(v) for v in value)
 if isinstance(value,(int,float)):return math.isfinite(value)
 return isinstance(value,(str,bool))


def checked_cholesky(matrix):
 """Preserve the intended quadratic form on the actual input matrix.

 Some library builds return unused-triangle contents. Never use those entries
 in a dense multiplication. These are factorization accuracy checks, not new
 scientific acceptance tolerances.
 """
 matrix=np.asarray(matrix,dtype=float)
 if matrix.ndim!=2 or matrix.shape[0]!=matrix.shape[1] or not np.isfinite(matrix).all():
  raise ValueError('Invalid precision/covariance matrix')
 if not np.allclose(matrix,matrix.T,rtol=1e-12,atol=1e-12):
  raise ValueError('Nonsymmetric precision/covariance matrix')
 factor=np.tril(cholesky(matrix,lower=True))
 if not np.isfinite(factor).all() or not np.allclose(factor@factor.T,matrix,rtol=1e-10,atol=1e-10):
  raise ValueError('Cholesky reconstruction failed')
 return factor


def profile_diagnostics(row,profile):
 """Validate full curves before accepting a summary of those curves."""
 try:
  if not finite_tree(row) or not finite_tree(profile):return False
  keys={'z','lower025','median','upper975'}
  if set(profile)!=keys:return False
  arrays={key:np.asarray(profile[key],dtype=float) for key in keys}
  count=row['NEC_grid_points']
  if isinstance(count,bool) or not isinstance(count,int) or count<2:return False
  if any(a.shape!=(count,) or not np.isfinite(a).all() for a in arrays.values()):return False
  if not np.all(np.diff(arrays['z'])>0):return False
  if not np.all(arrays['lower025']<=arrays['median']) or not np.all(arrays['median']<=arrays['upper975']):return False
  minimum=float(arrays['upper975'].min())
  return minimum==row['NEC_min_upper975'] and isinstance(row['NEC_no_negative_upper975'],bool) and row['NEC_no_negative_upper975']==(minimum>=0)
 except (KeyError,TypeError,ValueError,OverflowError):return False
