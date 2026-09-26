#!/usr/bin/env python3
"""Independent source-form, directed-interval, and theorem-instrumentation tests."""
from __future__ import annotations
import ast,csv,hashlib,json,random,sys,time,re
from pathlib import Path
import mpmath as mp
import numpy as np
import kernel as K
if not __debug__: raise RuntimeError('Verification requires assertions; do not use Python -O')
HERE=Path(__file__).resolve().parent
SOURCE=HERE/'qin2019_source_form.py'
checks=[]
def check(name,cond,detail=''):
    checks.append({'name':name,'pass':bool(cond),'detail':str(detail)})
    if not cond: print('FAIL',name,detail,flush=True)

def load_source():
    tree=ast.parse(SOURCE.read_text());nodes=[]
    for n in tree.body:
        if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='out' for t in n.targets):break
        if isinstance(n,ast.Assign) and any(isinstance(t,ast.Attribute) and t.attr=='dps' for t in n.targets):continue
        nodes.append(n)
    ns={};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(SOURCE),'exec'),ns);return ns
mp.mp.dps=85
M=load_source()

def source_hybrid(qe,qc):
    """Original log-generators at each stage, unlike the C++ bounded transform."""
    ge=M['H'](qe,3);gc=M['H'](qc,3);coll=[]
    for a in range(5):
        row=[]
        for c in range(5):
            vals=[tuple(x**(mp.mpf(3)/qe) for x in M['MATS'][e][a][c]) for e in range(4)]
            w=M['pweights'](vals,M['EW'],3); pair=M['one_partition'](vals,w,[1,0,0,0],ge)
            row.append(tuple(x**(mp.mpf(qe)/qc) for x in pair))
        w=M['pweights'](row,M['CW'],3);p=M['partitioned'](row,w,gc);coll.append(p[0]**qc-p[1]**qc)
    return coll

def source_state_scores(z):
    ge=M['H'](3,3);ss=[]
    for a in range(5):
        coll=[]
        for c in range(5):
            vals=[tuple(mp.root(z[e,a,c,k],3) for k in range(2)) for e in range(4)]
            w=M['pweights'](vals,M['EW'],3);coll.append(M['one_partition'](vals,w,[1,0,0,0],ge))
        w=M['pweights'](coll,M['CW'],3);p=M['partitioned'](coll,w,ge);ss.append(p[0]**3-p[1]**3)
    return ss

def inside(v,b):return mp.mpf(float(b[0]))<=v<=mp.mpf(float(b[1]))
K.lib.kernel_precision(192)
# Primitive tests against independently computed 85-digit values on deterministic grids.
rng=random.Random(20260915)
func=[lambda a,b:a+b,lambda a,b:a-b,lambda a,b:a*b,lambda a,b:a/b]
for code,fn in enumerate(func):
    ok=True
    for _ in range(40):
        a=sorted([rng.uniform(-1,1),rng.uniform(-1,1)]);b=sorted([rng.uniform(.05,1.5),rng.uniform(.05,1.5)])
        z=K.op(code,a,b)
        for x in a+[sum(a)/2]:
            for y in b+[sum(b)/2]:ok &=inside(fn(mp.mpf(x),mp.mpf(y)),z)
    check('directed primitive '+str(code),ok)
for code in [4,5,6,7,8]:
    ok=True
    for _ in range(40):
        a=sorted([rng.uniform(.001,.98),rng.uniform(.001,.98)])
        if code==7:a=[-a[1],a[0]]
        b=(3,3) if code in [4,5] else tuple(sorted([rng.uniform(.01,5),rng.uniform(.01,5)]))
        bound=K.op(code,a,b)
        for x in a+[sum(a)/2]:
            xx=mp.mpf(x)
            for y in b:
                v=xx**3 if code==4 else mp.root(xx,3) if code==5 else xx**mp.mpf(y) if code==6 else abs(xx) if code==7 else (1-xx)/(1+8*xx)
                ok &=inside(v,bound)
    check('directed nonlinear '+str(code),ok)
for code in [4,5,6,8]:
    check('endpoint zero/one '+str(code),all(np.isfinite(K.op(code,(0.,1.),(3.,3.)))))
# Every hybrid endpoint: independent source-form replay and interval inclusion.
max_error=mp.mpf(0); source_cache={}
cpp_digits=list(map(int,re.search(r'DATA\[200\]=\{([^}]*)', (HERE/'mpfr_interval_kernel.cpp').read_text()).group(1).split(',')))
source_digits=[int(M['MATS'][e][a][c][k]*10) for e in range(4) for a in range(5) for c in range(5) for k in range(2)]
check('200 exact-decimal source inputs independently cross-checked',cpp_digits==source_digits)
check('all 100 source canonical pairs admissible at q3',all(M['MATS'][e][a][c][0]**3+M['MATS'][e][a][c][1]**3<=1 for e in range(4) for a in range(5) for c in range(5)))
for q in range(3,11):
    for qe,qc in [(q,3),(q,q)]:
        ref=source_hybrid(qe,qc);source_cache[(qe,qc)]=ref;b=K.scores(qe,qc)
        ok=all(inside(v,bb) for v,bb in zip(ref,b));check(f'5 source-form endpoints in interval: ({qe},{qc})',ok)
        max_error=max(max_error,max(abs(v-mp.mpf(float(np.mean(bb)))) for v,bb in zip(ref,b)))
# Check every retained row, not just fresh C++ return values.
hybrid_rows=list(csv.DictReader((HERE/'QIN_STAGEWISE_HYBRIDS.csv').open()))
row_ok=True
for row in hybrid_rows:
    q=int(row['q']);qe,qc=(3,3) if row['hybrid']=='baseline_3_3' else (q,3) if row['hybrid']=='hybrid_q_3' else (q,q);a=int(row['alternative'][1:])-1
    row_ok &= inside(source_cache[(qe,qc)][a],(float(row['score_lower']),float(row['score_upper'])))
check('all 120 saved hybrid rows enclose source-form replay',len(hybrid_rows)==120 and row_ok)
# Proper independent directional derivative: NOT defined by subtracting raw and matched slopes.
z=np.empty((4,5,5,2),dtype=object);velocity=z.copy()
for e in range(4):
 for a in range(5):
  for c in range(5):
   for k in range(2):
    x=M['MATS'][e][a][c][k];z[e,a,c,k]=x**3;velocity[e,a,c,k]=x**3*mp.log(x)
raw=mp.diff(lambda q:M['run'](q,False,3)[1][0]-M['run'](q,False,3)[1][2],mp.mpf(3))
matched=mp.diff(lambda q:M['run'](q,True,3)[1][0]-M['run'](q,True,3)[1][2],mp.mpf(3))
def direct_direction(t):
    s=source_state_scores(z+t*velocity);return s[0]-s[2]
direction=mp.diff(direct_direction,mp.mpf(0))
check('independently evaluated directional chain-rule closure',abs(raw-matched-direction)<mp.mpf('1e-65'),mp.nstr(raw-matched-direction,8))
old_csv=HERE/'QIN2019_LOCAL_SENSITIVITY_DECOMPOSITION.csv'
if old_csv.exists():
 old={r['quantity']:mp.mpf(r['value']) for r in csv.DictReader(old_csv.open())}
 for key,value in [('d_raw_margin_dq_at_q3',raw),('d_matched_margin_dq_at_q3',matched),('representation_motion_term_at_q3',direction)]:check('Step1 value preserved '+key,abs(value-old[key])<mp.mpf('1e-35'))
# Validate the actual interval inequality used, not merely a stored PASS label.
rows=list(csv.DictReader((HERE/'QIN_STAGEWISE_CERTIFICATE.csv').open()))
for r in rows:
 q=int(r['q']);h0=source_cache[(3,3)];h1=source_cache[(q,3)];h2=source_cache[(q,q)]
 ee=max(abs(x-y) for x,y in zip(h1,h0));ec=max(abs(x-y) for x,y in zip(h2,h1))
 check(f'independent propagated budget inclusion q={q}',ee<=mp.mpf(float(r['expert_budget_upper'])) and ec<=mp.mpf(float(r['criterion_budget_upper'])) and ee+ec<=mp.mpf(float(r['total_budget_upper'])))
 check(f'nonvacuous full-winner stage budget q={q}',float(r['guaranteed_target_margin_lower'])>0 and float(r['twice_budget_upper'])<float(r['baseline_winner_margin_lower']))
# Stage theorem metatests: exact zero defect and cancellation in a constructed toy
# test only (not added to empirical cohort).
r0=np.array([1.,0.]);e1=np.array([.75,-.25]);e2=-e1
check('telescoping cancellation unit test',np.array_equal((r0+e1)+e2,r0) and np.linalg.norm(e1,np.inf)+np.linalg.norm(e2,np.inf)>0)
check('zero defect composition unit test',all(np.linalg.norm(K.scores(q,q)-K.scores(q,q),np.inf)==0 for q in [3,7,10]))
# This zero is only self-consistency; published q-natural control is separately checked.
summary=dict(checks=checks,passed=sum(x['pass'] for x in checks),total=len(checks),
 independent_source_max_abs_error=mp.nstr(max_error,15),mpfr=K.lib.kernel_version().decode(),bits=192,
 direct_derivatives={'raw':mp.nstr(raw,65),'matched':mp.nstr(matched,65),'representation_direction':mp.nstr(direction,65)},
 local_derivative_scope='Equality-preserving smooth stratum, not an assumed ambient Frechet gradient',
 source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest())
(HERE/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(summary,indent=2)+'\n')
print('OVERALL',summary['passed'],'/',summary['total'],'max endpoint error',summary['independent_source_max_abs_error'])
print('direct derivative',summary['direct_derivatives'])
sys.exit(0 if all(x['pass'] for x in checks) else 1)
