#!/usr/bin/env python3
"""Validated finite-rung, fixed-input stage-replacement audit (not an input radius)."""
import csv,json
from pathlib import Path
import numpy as np
import kernel as K
if not __debug__: raise RuntimeError('Verification requires assertions; do not use Python -O')
HERE=Path(__file__).resolve().parent
K.lib.kernel_precision(192)
base=K.scores(3,3)

def delta(a,b):
    return np.array([K.op(1,aa,bb) for aa,bb in zip(a,b)])
def norm_hi(x): return max(K.op(7,y)[1] for y in x)
def difflo(x,y): return K.op(1,(x,x),(y,y))[0]
margin0=min(K.op(1,base[0],base[j])[0] for j in range(1,5))
rows=[];hybrids=[]
for q in range(3,11):
    h1=K.scores(q,3); h2=K.scores(q,q)
    dE=delta(h1,base);dC=delta(h2,h1)
    if q==3: # structurally identical function calls, rather than dependency overestimate
        dE[:]=0;dC[:]=0
    eE=norm_hi(dE);eC=norm_hi(dC);budget=K.op(0,(eE,eE),(eC,eC))[1]
    floor=K.op(1,(margin0,margin0),K.op(2,(2,2),(budget,budget)))[0]
    dqE=K.op(1,dE[0],dE[2]);dqC=K.op(1,dC[0],dC[2])
    rows.append(dict(q=q,expert_budget_upper=eE,criterion_budget_upper=eC,total_budget_upper=budget,
                     twice_budget_upper=K.op(2,(2,2),(budget,budget))[1],baseline_winner_margin_lower=margin0,
                     guaranteed_target_margin_lower=floor,expert_A1A3_delta_lower=dqE[0],expert_A1A3_delta_upper=dqE[1],
                     criterion_A1A3_delta_lower=dqC[0],criterion_A1A3_delta_upper=dqC[1],
                     actual_A1A3_margin_lower=K.op(1,h2[0],h2[2])[0],winner_certified=floor>0))
    for name,s in [('baseline_3_3',base),('hybrid_q_3',h1),('matched_q_q',h2)]:
        for a,(l,h) in enumerate(s):hybrids.append(dict(q=q,hybrid=name,alternative=f'A{a+1}',score_lower=l,score_upper=h))
for name,data in [('QIN_STAGEWISE_CERTIFICATE.csv',rows),('QIN_STAGEWISE_HYBRIDS.csv',hybrids)]:
    with (HERE/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=data[0].keys());w.writeheader();w.writerows(data)
print('MPFR',K.lib.kernel_version().decode(),'precision=192 bits')
for r in rows: print(r['q'], 'budget',r['total_budget_upper'],'guaranteed margin',r['guaranteed_target_margin_lower'],'E margin delta',r['expert_A1A3_delta_lower'],'C',r['criterion_A1A3_delta_lower'])
assert all(r['winner_certified'] for r in rows)
summary=dict(status='PASS',scope='q=3,...,10 only; fixed exact printed inputs; source-order expert-then-criterion hybrid path',rungs=8,
             max_budget_upper=max(r['total_budget_upper'] for r in rows),min_margin_lower=min(r['guaranteed_target_margin_lower'] for r in rows),
             baseline_margin_lower=margin0,library='MPFR '+K.lib.kernel_version().decode(),bits=192,
             exact_global_input_radius='NOT_CLAIMED',general_propagation_novelty='NOT_CLAIMED')
(HERE/'STAGEWISE_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
