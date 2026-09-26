#!/usr/bin/env python3
"""Existing He control, chart relabeling, and Qin differentiability-domain audit."""
from __future__ import annotations
import ast,csv,json,hashlib,sys
from pathlib import Path
import numpy as np
import mpmath as mp
import he2019_archived_evaluator as H
if not __debug__: raise RuntimeError('Verification requires assertions; do not use Python -O')
HERE=Path(__file__).resolve().parent
checks=[]
def ck(name,ok,detail=''):
 checks.append({'name':name,'pass':bool(ok),'detail':str(detail)})
 if not ok:print('FAIL',name,detail,flush=True)
base=H.full_run(3,True);rows=[]
for q in range(1,11):
 coll,fin,s,rank=H.full_run(q,True)
 at_r=coll**(q/3)
 hybrid=np.array([H.rpfdwhm(at_r[a],H.CW,3) for a in range(5)])
 hs=hybrid[:,0]**3+1-hybrid[:,2]**3
 dE=hs-base[2];dC=s-hs
 rows.append({'q':q,'matched_rank':rank,'expert_canonical_drift':float(np.max(abs(coll**q-base[0]**3))),
  'propagated_expert_score_drift':float(np.max(abs(dE))),'criterion_score_drift':float(np.max(abs(dC))),
  'max_matched_score_drift':float(np.max(abs(s-base[2])))})
 ck(f'He existing q-natural control q={q}',rank==H.BASE_RANK and max(rows[-1][x] for x in ['expert_canonical_drift','propagated_expert_score_drift','criterion_score_drift','max_matched_score_drift'])<5e-14)
with (HERE/'HE_STAGEWISE_CONTROL.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
# Chart relabeling preserves transport, but not an arbitrary newly imposed metric.
x=np.array([.7,.2]);r,q=3,7;phi=lambda a,p:a**p;psi=lambda a:a**2
original=phi(x,r)**(1/q);relabelled=np.sqrt(psi(phi(x,r)))**(1/q)
ck('fixed common relabeling transport unchanged',np.max(abs(original-relabelled))<1e-15)
z=np.array([.3,.1]);zz=np.array([.4,.2])
ck('transport does not imply arbitrary chart Euclidean metric invariance',abs(np.linalg.norm(z-zz)-np.linalg.norm(psi(z)-psi(zz)))>.01)
# Source-form functions loaded without running archived file's top-level report code.
mp.mp.dps=70;tree=ast.parse((HERE/'qin2019_source_form.py').read_text());nodes=[]
for n in tree.body:
 if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='out' for t in n.targets):break
 if isinstance(n,ast.Assign) and any(isinstance(t,ast.Attribute) and t.attr=='dps' for t in n.targets):continue
 nodes.append(n)
M={};exec(compile(ast.Module(body=nodes,type_ignores=[]),'qin2019_source_form.py','exec'),M)
X=np.array(M['MATS'],dtype=object);Z=X**3;V=np.empty_like(Z)
for index in np.ndindex(X.shape):V[index]=Z[index]*mp.log(X[index])
coincident=[]
for a in range(5):
 for c in range(5):
  for e in range(4):
   for f in range(e+1,4):
    if all(X[e,a,c,k]==X[f,a,c,k] for k in range(2)):
     coincident.append({'alternative':a+1,'criterion':c+1,'experts':[e+1,f+1]})
     ck(f'Qin canonical velocity tangent to equality A{a+1}C{c+1}E{e+1}E{f+1}',all(V[e,a,c,k]==V[f,a,c,k] for k in range(2)))
ck('six source input coincidences explicitly retained',len(coincident)==6)
gen=M['H'](3,3)
def score(z,a,return_coll=False):
 coll=[]
 for c in range(5):
  vals=[tuple(mp.root(z[e,a,c,k],3) for k in range(2)) for e in range(4)]
  w=M['pweights'](vals,M['EW'],3);coll.append(M['one_partition'](vals,w,[1,0,0,0],gen))
 if return_coll:return coll
 w=M['pweights'](coll,M['CW'],3);p=M['partitioned'](coll,w,gen);return p[0]**3-p[1]**3
for a in [0,2]:
 coll=score(Z,a,True)
 separation=min(sum(abs(coll[c][k]-coll[d][k])**3 for k in range(2)) for c in range(5) for d in range(c+1,5))
 ck(f'Qin criterion-support distinctness near baseline A{a+1}',separation>0,mp.nstr(separation,20))
# Diagnostic of a normal direction breaking one equality; not a proof by numerics.
g0=score(Z,0)-score(Z,2);one_sided=[]
for h in [mp.mpf('1e-5'),mp.mpf('1e-7'),mp.mpf('1e-9')]:
 zp=Z.copy();zm=Z.copy();zp[2,0,1,0]+=h;zm[2,0,1,0]-=h
 right=(score(zp,0)-score(Z,2)-g0)/h
 left=(g0-(score(zm,0)-score(Z,2)))/h
 one_sided.append({'h':str(h),'right_slope':mp.nstr(right,25),'left_slope':mp.nstr(left,25),'difference':mp.nstr(right-left,25)})
summary={'passed':sum(x['pass'] for x in checks),'total':len(checks),'checks':checks,
 'He_scope':'The archived table-consistent implementation, including its documented source Table-4 discrepancy; no new global-optimality claim.',
 'Qin_equality_pairs':coincident,'normal_direction_one_sided_diagnostic':one_sided,
 'local_derivative_scope':'Raw/matched q paths and the canonical representation-motion direction preserve all source equalities. Differentiate on this equality-preserving stratum. An ambient Frechet gradient is not asserted.',
 'input_radius_scope':'The MPFR radius calculation is derivative-free and does not restrict perturbations to this stratum.'}
(HERE/'STRUCTURAL_CONTROL_VERIFICATION.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2));sys.exit(not all(x['pass'] for x in checks))
