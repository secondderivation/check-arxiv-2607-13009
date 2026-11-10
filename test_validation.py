import unittest,tempfile,json
from pathlib import Path
import numpy as np
from scipy.linalg import cholesky,solve_triangular
from validation import digest,verify_manifest,compare,diagnostics
from fetch_data import fetch
class Controls(unittest.TestCase):
 def row(self):return dict(ESS=5000.,max_weight=.001,exact_delta_chi2_max=.002,H047_mean=60.,H047_sd=.5,NEC_min_upper975=.55,NEC_no_negative_upper975=True)
 def test_reject_scientific_and_replay_failures(self):
  r=self.row();self.assertTrue(compare(r,r))
  for key,value in dict(ESS=500.,max_weight=.1,exact_delta_chi2_max=.02,H047_mean=61.,H047_sd=.6,NEC_min_upper975=-.1,NEC_no_negative_upper975=False).items():
   self.assertFalse(compare({**r,key:value},r),key)
  self.assertFalse(diagnostics({**r,'ESS':float('nan')}))
 def test_covariance_precision_and_marginalization(self):
  C=np.array([[2.,.4],[.4,1.]]);P=np.linalg.inv(C);v=np.array([1.,-2.]);one=np.ones(2)
  a=solve_triangular(cholesky(C,lower=True),v,lower=True);L=cholesky(P,lower=True);b=L.T@v
  self.assertAlmostEqual(a@a,b@b,places=12)
  w=L.T@one;projected=b-w*(w@b)/(w@w);closed=v@P@v-(one@P@v)**2/(one@P@one)
  self.assertAlmostEqual(projected@projected,closed,places=12)
  # Wrong precision orientation must be detectable with this correlated input.
  self.assertGreater(abs((L@v)@(L@v)-v@P@v),.01)
 def test_damaged_and_extra_payload(self):
  with tempfile.TemporaryDirectory() as t:
   root=Path(t);p=root/'reference.json';p.write_text('{}');(root/'MANIFEST.txt').write_text(digest(p)+'  reference.json\n');verify_manifest(root)
   p.write_text('{"changed":true}')
   with self.assertRaises(ValueError):verify_manifest(root)
   p.write_text('{}');(root/'private-report.md').write_text('private')
   with self.assertRaises(ValueError):verify_manifest(root)
 def test_damaged_data(self):
  with tempfile.TemporaryDirectory() as t:
   files=json.loads((Path(__file__).parent/'datasets.json').read_text())['files'];p=Path(t)/files[0]['path'];p.parent.mkdir(parents=True);p.write_text('corrupt')
   with self.assertRaises(ValueError):fetch(t)

class RepairControls(unittest.TestCase):
 def test_factor_triangle_and_fail_closed(self):
  from unittest.mock import patch
  from validation import checked_cholesky
  P=np.array([[2.,.4],[.4,1.]])
  true=np.linalg.cholesky(P);dirty=true.copy();dirty[0,1]=123.
  with patch('validation.cholesky',return_value=dirty):
   L=checked_cholesky(P);np.testing.assert_allclose(L@L.T,P,rtol=1e-12,atol=1e-12)
   self.assertEqual(L[0,1],0.)
  with patch('validation.cholesky',return_value=np.eye(2)):
   with self.assertRaises(ValueError):checked_cholesky(P)
  for bad in [np.array([[2.,float('nan')],[float('nan'),1.]]),np.array([[2.,.4],[.5,1.]])]:
   with self.assertRaises(ValueError):checked_cholesky(bad)
 def test_all_numeric_fields_finite(self):
  import copy
  row=json.loads((Path(__file__).parent/'reference/I4_RA_731009.json').read_text())
  for key in ['H0','rd_fixed','h_plugin','map','map_cov','half68','shortest68_half','corr_h_H','chi2','NEC_grid_points']:
   for bad in [float('nan'),float('inf'),-float('inf')]:
    changed=copy.deepcopy(row);changed[key]=[[bad]] if key=='map_cov' else [bad,1.] if key in ['H0','rd_fixed','h_plugin','map'] else bad
    self.assertFalse(diagnostics(changed),(key,bad))
  self.assertFalse(diagnostics({**row,'exact_delta_chi2_max':-1.}))
 def test_full_profile_rejection(self):
  import copy
  from validation import profile_diagnostics
  root=Path(__file__).parent/'reference';row=json.loads((root/'I4_RA_731009.json').read_text());profile=json.loads((root/'I4_RA_731009_nec.json').read_text())
  self.assertTrue(profile_diagnostics(row,profile))
  for field in profile:
   for bad in [float('nan'),float('inf'),-float('inf')]:
    changed=copy.deepcopy(profile);changed[field][1]=bad
    self.assertFalse(profile_diagnostics(row,changed),(field,bad))
  changed=copy.deepcopy(profile);changed['z'][1]=changed['z'][0];self.assertFalse(profile_diagnostics(row,changed))
  changed=copy.deepcopy(profile);changed['upper975'].pop();self.assertFalse(profile_diagnostics(row,changed))
  changed=copy.deepcopy(profile);changed['lower025'][2]=changed['upper975'][2]+1.;self.assertFalse(profile_diagnostics(row,changed))
  changed=copy.deepcopy(profile);del changed['median'];self.assertFalse(profile_diagnostics(row,changed))
  self.assertFalse(profile_diagnostics({**row,'NEC_min_upper975':row['NEC_min_upper975']+1.},profile))

if __name__=='__main__':unittest.main()
