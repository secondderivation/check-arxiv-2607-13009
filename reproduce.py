#!/usr/bin/env python3
"""Full immutable replay; exits nonzero for any failed acceptance check."""
import os,sys
os.environ['PYTHONDONTWRITEBYTECODE']='1';sys.dont_write_bytecode=True
for key in ['OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','OMP_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import json, argparse, datetime, platform, subprocess
from validation import verify_manifest,RUNS,diagnostics,compare,seed_agreement,profile_diagnostics
from fetch_data import fetch
ROOT=Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser();p.add_argument('--all',action='store_true',required=True);a=p.parse_args()
 verify_manifest(ROOT)
 stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ');out=ROOT/'_runs'/stamp;out.mkdir(parents=True)
 os.environ['SD_OUTPUT']=str(out);os.environ['SD_DATA']=str(ROOT/'_runs/data')
 count=fetch(os.environ['SD_DATA'])
 subprocess.run([sys.executable,'-B','-m','unittest','test_validation.py'],cwd=ROOT,check=True)
 from checks.algebra import run
 from checks.reconstruction import Model
 algebra=run();(out/'algebra.json').write_text(json.dumps(algebra,indent=2)+'\n')
 rows={}
 for name in RUNS:
  closure,prior,seed=name.split('_')
  rows[name]=Model(closure=closure,pressure=prior=='REIRA').integrate(name,seed=int(seed))
  print(name, 'ESS',round(rows[name]['ESS']), 'NEC minimum upper bound',round(rows[name]['NEC_min_upper975'],6),flush=True)
 references={name:json.loads((ROOT/'reference'/f'{name}.json').read_text()) for name in RUNS}
 profiles={name:json.loads((out/f'{name}_nec.json').read_text()) for name in RUNS}
 reference_profiles={name:json.loads((ROOT/'reference'/f'{name}_nec.json').read_text()) for name in RUNS}
 checks={'full_profiles':all(profile_diagnostics(rows[n],profiles[n]) and profile_diagnostics(references[n],reference_profiles[n]) for n in RUNS),'algebra':algebra['passed'],'numerical_diagnostics':all(diagnostics(r) for r in rows.values()),'reference_replay':all(compare(rows[n],references[n]) for n in RUNS),'independent_seed_agreement':seed_agreement(rows)}
 verify_manifest(ROOT)
 import numpy,scipy,sympy
 result={'passed':all(checks.values()),'checks':checks,'data_files_verified':count,'runs':len(rows),'python':platform.python_version(),'platform':platform.platform(),'numpy':numpy.__version__,'scipy':scipy.__version__,'sympy':sympy.__version__,'approved':False,'published':False}
 (out/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));print('Output:',out)
 if not result['passed']:raise SystemExit(1)
if __name__=='__main__':main()
