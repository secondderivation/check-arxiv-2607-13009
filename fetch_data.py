"""Fetch public upstream inputs; fail closed on a changed file."""
from pathlib import Path
import json, urllib.request, argparse
from validation import digest
ROOT=Path(__file__).resolve().parent

def fetch(destination):
 destination=Path(destination)
 for f in json.loads((ROOT/'datasets.json').read_text())['files']:
  p=destination/f['path'];p.parent.mkdir(parents=True,exist_ok=True)
  if not p.exists():
   temporary=p.with_name(p.name+'.download')
   with urllib.request.urlopen(f['url'],timeout=120) as response: temporary.write_bytes(response.read())
   if digest(temporary)!=f['sha256']:
    temporary.unlink();raise ValueError('Downloaded data hash mismatch: '+f['path'])
   temporary.replace(p)
  if digest(p)!=f['sha256']:raise ValueError('Data hash mismatch: '+f['path'])
 return len(json.loads((ROOT/'datasets.json').read_text())['files'])
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--destination',default='_runs/data');a=p.parse_args();print(fetch(a.destination),'verified files')
