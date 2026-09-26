"""Finite implementation diagnostics, not a proof of transcendental correctness."""
from pathlib import Path
from fractions import Fraction as F
import copy,json,random,sys,time
HERE=Path(__file__).resolve().parent
try:import mpmath as mp
except ImportError as exc:
 raise RuntimeError('Preinstalled mpmath is required for diagnostics; add its directory to PYTHONPATH or use reproduce_v17.py --python-deps PATH') from exc
import verify as v
mp.mp.dps=110
def m(x):
 x=F(x);return mp.mpf(x.numerator)/x.denominator
def contained(value,bounds):return m(bounds[0])<=value<=m(bounds[1])
def source_scores(raw,q):
 # Original source-form transport and qROFAWA; deliberately not transformed
 # weighted logarithms used by verify.py.
 q=m(q);x=[m(z)**(4/q) for z in raw]
 entropy=[1+sum(x[(i*5+j)*2+c]*mp.log(x[(i*5+j)*2+c]) for i in range(4) for c in (0,1))/4 for j in range(5)]
 w=[a/sum(entropy) for a in entropy];scores=[]
 for i in range(4):
  membership=(1-mp.fprod((1-x[(i*5+j)*2]**q)**w[j] for j in range(5)))**(1/q)
  nonmembership=mp.fprod(x[(i*5+j)*2+1]**w[j] for j in range(5))
  scores.append(membership**q-nonmembership**q)
 return scores
def main():
 start=time.perf_counter();checks=[]
 # Both monotone branches, crossing the turning point, critical singleton,
 # zero limiting endpoint, and a tiny interval around the critical point.
 cases=[(F(-4),F(-2)),(F(-2),F(-1,2)),(F(-9,10),F(0)),(F(-1),F(-1)),(F(0),F(0)),(F(-1000001,1000000),F(-999999,1000000))]
 for a,b in cases:
  bounds=v.yexp_range(a,b)
  for k in range(65):
   y=m(a+(b-a)*F(k,64));v.require(contained(y*mp.exp(y),bounds),'yexp range miss')
  checks.append({'check':'yexp_branch','input':[str(a),str(b)],'samples':65,'enclosure':[str(x) for x in bounds]})
 for a,b in [(F(1,10),F(1,5)),(F(0),F(-1)),(F(-1),F(1))]:
  try:v.yexp_range(a,b)
  except (ValueError,OverflowError):pass
  else:raise ValueError('invalid yexp domain accepted')
 raw=v.baseline();rng=random.Random(20260926)
 for q in [F(4),F(8),F(16),F(7,1),F(1000),F(100000000)]:
  for attempt in range(3):
   sample=[x+F(rng.randint(-100,100),10000) for x in raw] if attempt else list(raw)
   intervals,_=v.box_scores(1/q,1/q,[(x,x) for x in sample]);direct=source_scores(sample,q)
   v.require(all(contained(x,b) for x,b in zip(direct,intervals)),'source-form point mismatch')
  checks.append({'check':'source_form_sample','q':str(q),'samples':3,'precision_decimal_digits':mp.mp.dps})
 for q in (F(4),F(8),F(16)):
  intervals,_=v.scores(1/q,1/q,F(1,100))
  for attempt in range(8):
   sample=[x+F(rng.randint(-100,100),10000) for x in raw]
   direct=source_scores(sample,q)
   v.require(all(contained(x,b) for x,b in zip(direct,intervals)),'source-box sample mismatch')
  checks.append({'check':'source_box_sample','q':str(q),'halfwidth':'1/100','samples':8})
 # Make a small independently computed valid record, then corrupt each of its
 # structural/arithmetic fields in isolation. Final integer endpoint replay is
 # a separate responsibility, not asserted by these diagnostic checks.
 sc,_=v.scores(F(1,4),F(1,4),F(0));margins=[sc[1][0]-sc[j][1] for j in (0,2,3)]
 n={'id':'r','t':['1/4','1/4'],'scores':[[str(x),str(y)] for x,y in sc],
    'margins':list(map(str,margins)),'positive':True}
 record={'status':'PASS_COMPUTATIONAL','version':17,'enclosure':'scalar_yexp_endpoint_critical+weight_simplex_LP',
         'baseline_rung':4,'winner':'Y2','root_t':['1/4','1/4'],'q_exact':'4',
         'normalized_baseline_raw_halfwidth':'0','input_model_sha256':v.INPUT_SHA,
         'unresolved_count':0,'unresolved':[],'node_count':1,'leaf_count':1,
         'nodes':[n],'leaves':[copy.deepcopy(n)],'margin_floor':str(min(margins))}
 v.check_record(record)
 mutations={
 'version':lambda d:d.update(version=16),
 'enclosure':lambda d:d.update(enclosure='unknown'),
 'q_exact':lambda d:d.update(q_exact='8'),
 'uniform_domain':lambda d:d.update(q_exact=None),
 'status':lambda d:d.update(status='UNRESOLVED'),
 'baseline':lambda d:d.update(baseline_rung=3),
 'winner':lambda d:d.update(winner='Y1'),
 'model':lambda d:d.update(input_model_sha256='0'*64),
 'halfwidth':lambda d:d.update(normalized_baseline_raw_halfwidth='-1/100'),
 'root':lambda d:d.update(root_t=['0','1/2']),
 'unresolved':lambda d:d.update(unresolved_count=1),
 'node_count':lambda d:d.update(node_count=2),
 'missing_leaf':lambda d:d.update(leaves=[],leaf_count=0),
 'node_t':lambda d:d['nodes'][0].update(t=['0','1/4']),
 'positivity_type':lambda d:d['nodes'][0].update(positive=1),
 'positivity':lambda d:d['nodes'][0].update(positive=False),
 'margin':lambda d:d['nodes'][0].update(margins=['1','1','1']),
 'leaf_copy':lambda d:d['leaves'][0].update(positive=False),
 'floor':lambda d:d.update(margin_floor='1'),
 'weights':lambda d:d['nodes'][0].update(weights=[['1','1']]*5)}
 for name,change in mutations.items():
  d=copy.deepcopy(record);change(d)
  try:v.check_record(d)
  except (ValueError,KeyError):pass
  else:raise ValueError('corrupted record accepted: '+name)
 checks.append({'check':'negative_record_tests','rejected':list(mutations),'count':len(mutations)})
 result={'status':'PASS','scope':'Finite branch/sample/negative-record diagnostics; not an exhaustive proof',
         'mpfr_version':v.library().kernel_version().decode(),'mpfr_bits':256,
         'mpmath_version':mp.__version__,'elapsed_seconds':time.perf_counter()-start,'checks':checks}
 (HERE/'TESTS.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({k:r for k,r in result.items() if k!='checks'},indent=2))
if __name__=='__main__':main()
