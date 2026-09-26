#!/usr/bin/env python3
"""Reject malformed certificate headers and tree roots before accepting a cover."""
from pathlib import Path
import gzip,json,tempfile,shutil
from replay_certificate import replay
HERE=Path(__file__).resolve().parent
if not __debug__:raise RuntimeError('Do not use Python -O')
with gzip.open(HERE/'radius_000995/COVER_TREE.jsonl.gz','rt') as f:
 header=json.loads(next(f));root=json.loads(next(f))
cases=[]
for name in ['wrong_challenger','radius_below_declared','precision_mismatch','invalid_root_id']:
 with tempfile.TemporaryDirectory(prefix='qrof_negative_') as tmp:
  d=Path(tmp);shutil.copy2(HERE/'radius_000995/RADIUS_SUMMARY.json',d/'RADIUS_SUMMARY.json')
  h=header.copy();r=root.copy()
  if name=='wrong_challenger':h['challenger']=1
  if name=='radius_below_declared':h['radius_dyadic_hex']=float(0.0009).hex()
  if name=='precision_mismatch':h['bits']=64
  if name=='invalid_root_id':r['id']=1
  with gzip.open(d/'COVER_TREE.jsonl.gz','wt') as f:
   f.write(json.dumps(h)+'\n'+json.dumps(r)+'\n')
  rejected=False
  try:replay(d)
  except (AssertionError,ValueError,KeyError):rejected=True
  cases.append({'name':name,'rejected':rejected})
assert all(c['rejected'] for c in cases)
o={'status':'PASS','rejected':len(cases),'tested':len(cases),'tests':cases,'scope':'Malformed header/root rejection; not exhaustive adversarial software verification'}
(HERE/'CHECKER_NEGATIVE_VERIFICATION.json').write_text(json.dumps(o,indent=2)+'\n');print(json.dumps(o,indent=2))
