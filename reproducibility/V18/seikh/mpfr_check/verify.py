"""Second arithmetic for V17, independent of Python-integer transcendental code.

MPFR 256-bit primitives; binary64 endpoints always directed outwards; exact
Fraction LP vertex enumeration. Shared model definitions and source inputs mean
this is a second implementation, not independent peer review or formal proof.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,ctypes as C,hashlib,itertools,json,math,os
HERE=Path(__file__).resolve().parent
INPUT_SHA='71351416bb73b22f3781bad962a023b6ea7362e8020429aab97b47d3f36d1797'
P=C.POINTER(C.c_double)
_LIB=None;_DLL_DIRECTORY=None
def require(value,message):
 if not value:raise ValueError(message)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def baseline():
 p=HERE/'model_inputs.json';require(sha(p)==INPUT_SHA,'input model pin mismatch')
 d=json.loads(p.read_text())['seikh'];require(d['r']==4 and d['shape']==[4,5,2] and d['denominator']==100,'wrong model geometry')
 return tuple(F(n,d['denominator']) for n in d['numerators'])
def library():
 global _LIB,_DLL_DIRECTORY
 if _LIB is None:
  if os.name=='nt':
   configured=os.environ.get('QROF_MPFR_TOOLCHAIN')
   require(configured,'Windows requires QROF_MPFR_TOOLCHAIN pointing to the preinstalled MPFR/GMP prefix')
   toolchain=Path(configured).resolve()
   _DLL_DIRECTORY=os.add_dll_directory(str(toolchain/'bin'))
  _LIB=C.CDLL(str(HERE/('seikh_v17_components.dll' if os.name=='nt' else 'seikh_v17_components.so')))
  _LIB.kernel_precision.argtypes=[C.c_uint];_LIB.kernel_precision(256)
  _LIB.kernel_error.restype=C.c_char_p;_LIB.kernel_version.restype=C.c_char_p
  _LIB.v17_components.argtypes=[P,P,C.c_double,C.c_double,P,P];_LIB.v17_components.restype=C.c_int
  for name in ('v17_yexp','v17_exp'):
   getattr(_LIB,name).argtypes=[C.c_double,C.c_double,P];getattr(_LIB,name).restype=C.c_int
 return _LIB
def outward(x,upper):
 x=F(x);f=float(x);require(math.isfinite(f),'nonfinite conversion')
 if (upper and F(f)<x) or (not upper and F(f)>x):f=math.nextafter(f,math.inf if upper else -math.inf)
 return f
def unary_range(name,lo,hi):
 lib=library();out=(C.c_double*2)();rc=getattr(lib,name)(outward(lo,False),outward(hi,True),out)
 require(rc==0,name+' failed: '+lib.kernel_error().decode());return F(out[0]),F(out[1])
def yexp_range(lo,hi):return unary_range('v17_yexp',lo,hi)
def exp_range(lo,hi):return unary_range('v17_exp',lo,hi)
def vertices(bounds):
 result=set();n=len(bounds)
 require(all(0<l<=u<1 for l,u in bounds),'positive weight bounds required')
 for free in range(n):
  rest=[i for i in range(n) if i!=free]
  for mask in itertools.product((0,1),repeat=n-1):
   w=[F(0)]*n
   for i,k in zip(rest,mask):w[i]=bounds[i][k]
   w[free]=1-sum(w)
   if bounds[free][0]<=w[free]<=bounds[free][1]:result.add(tuple(w))
 require(result,'empty boxed simplex');return sorted(result)
def box_scores(a,b,raw_bounds):
 a,b=F(a),F(b);require(0<=a<=b<=F(1,4),'invalid t interval')
 require(len(raw_bounds)==40,'expected 40 normalized raw input components')
 raw_bounds=[(F(x),F(y)) for x,y in raw_bounds]
 require(all(0<x<=y<1 for x,y in raw_bounds),'strictly interior source bounds required')
 z=[(x**4,y**4) for x,y in raw_bounds]
 require(all(z[i][1]+z[i+1][1]<=1 for i in range(0,40,2)),'source box not wholly admissible')
 lo=(C.c_double*40)(*[outward(x,False) for x,y in z]);hi=(C.c_double*40)(*[outward(y,True) for x,y in z])
 eo=(C.c_double*10)();logs=(C.c_double*80)();lib=library()
 rc=lib.v17_components(lo,hi,outward(a,False),outward(b,True),eo,logs)
 require(rc==0,'components failed: '+lib.kernel_error().decode())
 e=[(F(eo[2*j]),F(eo[2*j+1])) for j in range(5)]
 w=[(l/(l+sum(e[k][1] for k in range(5) if k!=j)),u/(u+sum(e[k][0] for k in range(5) if k!=j))) for j,(l,u) in enumerate(e)]
 vs=vertices(w);out=[]
 for i in range(4):
  terms=[]
  for c in (0,1):
   coeff=[(F(logs[2*((i*5+j)*2+c)]),F(logs[2*((i*5+j)*2+c)+1])) for j in range(5)]
   lower=min(sum(x*coeff[j][0] for j,x in enumerate(v)) for v in vs)
   upper=max(sum(x*coeff[j][1] for j,x in enumerate(v)) for v in vs)
   terms.append(exp_range(lower,upper))
  out.append((1-terms[0][1]-terms[1][1],1-terms[0][0]-terms[1][0]))
 return out,{'weight_vertices':len(vs),'entropy':e,'weights':w}
def scores(a,b,h):
 h=F(h);require(h>=0,'negative halfwidth')
 result,info=box_scores(a,b,[(x-h,x+h) for x in baseline()]);return result,info['weight_vertices']
def witness(raw,q):
 raw=tuple(F(x) for x in raw);q=F(q);require(q>=4,'witness q<4')
 sc,info=box_scores(1/q,1/q,[(x,x) for x in raw])
 winners=[i for i in range(4) if sc[i][0]>max(sc[j][1] for j in range(4) if j!=i)]
 require(len(winners)==1 and winners[0]!=1,'failure witness has no certified different unique winner')
 dist=max(abs(x-c) for x,c in zip(raw,baseline()))
 return {'status':'PASS_MPFR','q':str(q),'raw_inputs':[str(x) for x in raw],
  'actual_linf_distance':str(dist),'scores':[[str(x),str(y)] for x,y in sc],
  'margin_Y2_minus_Y1':[str(sc[1][0]-sc[0][1]),str(sc[1][1]-sc[0][0])],
  'new_winner':'Y'+str(winners[0]+1),'all_input_pairs_admissible':True,
  'scope':'one admissible source-raw Linf witness at the stated finite q; not an optimum'}
def rational_pair(value,label):
 require(isinstance(value,list) and len(value)==2,label+' must be an endpoint pair')
 a,b=map(F,value);require(a<=b,label+' reversed');return a,b
def check_record(d):
 require(d['status']=='PASS_COMPUTATIONAL','incomplete certificate')
 require(d['version']==17 and d['enclosure']=='scalar_yexp_endpoint_critical+weight_simplex_LP','wrong version/enclosure')
 require(d['baseline_rung']==4 and d['winner']=='Y2','wrong target')
 require(d['input_model_sha256']==INPUT_SHA,'wrong pinned input model')
 h=F(d['normalized_baseline_raw_halfwidth']);require(h>=0,'negative halfwidth')
 root=rational_pair(d['root_t'],'root_t');require(0<=root[0]<=root[1]<=F(1,4),'root out of domain')
 if d['q_exact'] is None:require(root==(F(0),F(1,4)),'uniform certificate must cover all [0,1/4]')
 else:
  q=F(d['q_exact']);require(q>=4 and root==(1/q,1/q),'fixed-q domain mismatch')
 require(d['unresolved_count']==0 and d['unresolved']==[],'unresolved leaves')
 nodes={n['id']:n for n in d['nodes']};leaves={n['id']:n for n in d['leaves']}
 require(len(nodes)==len(d['nodes'])==d['node_count'] and len(leaves)==len(d['leaves'])==d['leaf_count'],'counts/duplicates')
 require(leaves and set(leaves)<=set(nodes),'missing leaf nodes')
 pending=[('r',*root)];seen=set();stored_floor=None
 while pending:
  ident,a,b=pending.pop();require(ident in nodes and ident not in seen,'missing/duplicate node');n=nodes[ident];seen.add(ident)
  require(rational_pair(n['t'],'node t')==(a,b),'wrong cover interval')
  require(type(n['positive']) is bool,'positive must be Boolean')
  sc=[rational_pair(s,'stored score') for s in n['scores']];require(len(sc)==4,'score count')
  margins=[sc[1][0]-sc[j][1] for j in (0,2,3)]
  require(list(map(F,n['margins']))==margins,'stored score/margin mismatch')
  require(n['positive']==(min(margins)>0),'stored positivity mismatch')
  if 'weights' in n:
   w=[rational_pair(x,'stored weight') for x in n['weights']]
   require(len(w)==5 and all(0<l<=u<=1 for l,u in w) and sum(l for l,u in w)<=1<=sum(u for l,u in w),'invalid stored weight enclosure')
  if ident in leaves:
   require(n==leaves[ident] and n['positive'],'leaf identity or positivity')
   stored_floor=min(margins) if stored_floor is None else min(stored_floor,*margins)
  else:
   require(not n['positive'] and a<b,'invalid split')
   m=(a+b)/2;pending.extend([(ident+'0',a,m),(ident+'1',m,b)])
 require(seen==set(nodes) and len(nodes)==2*len(leaves)-1,'invalid complete tree')
 require(F(d['margin_floor'])==stored_floor,'stored floor mismatch')
 rows=[];floor=None
 for leaf in sorted(leaves.values(),key=lambda n:F(n['t'][0])):
  a,b=map(F,leaf['t']);sc,count=scores(a,b,h);margins=[sc[1][0]-sc[j][1] for j in (0,2,3)]
  require(min(margins)>0,'MPFR failed to re-prove positive margin at '+leaf['id'])
  floor=min(margins) if floor is None else min(floor,*margins)
  rows.append({'id':leaf['id'],'t':leaf['t'],'scores':[[str(x),str(y)] for x,y in sc],
               'margins':[str(x) for x in margins],'weight_vertices':count})
 return {'status':'PASS','arithmetic':'MPFR 256-bit directed endpoints + exact Fraction vertex enumeration',
  'mpfr_version':library().kernel_version().decode(),'root_t':[str(x) for x in root],
  'normalized_baseline_raw_halfwidth':str(h),'input_model_sha256':INPUT_SHA,
  'leaf_count':len(rows),'floor_exact':str(floor),'floor_display':float(floor),'leaf_checks':rows,
  'validation_scope':'Strict tree and internal record consistency; independent MPFR positivity on every leaf. Stored integer endpoint enclosures require their own deterministic replay.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('record',type=Path);p.add_argument('--output',type=Path);a=p.parse_args()
 data=json.loads(a.record.read_text());result=check_record(data)
 if a.output:a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({k:v for k,v in result.items() if k!='leaf_checks'},indent=2))
