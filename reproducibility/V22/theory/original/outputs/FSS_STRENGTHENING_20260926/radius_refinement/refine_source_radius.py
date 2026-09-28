"""Experimental monotonicity-assisted certificate, unchanged Seikh model.

All proof decisions use the V18 hash-pinned integer interval kernel. A successful
local optimizer is never used as a lower-bound proof. Source coordinates and
L-infinity geometry are identical to V17/V18. Standard interval monotonicity and
centered forms are numerical techniques, not a claim of theoretical novelty.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1] / "FSS_V18_RELEASE_20260926/public_artifact/seikh"
sys.path.insert(0, str(MODEL))
from seikh_certificates import IV, SCALE, PREC, entropy_term, scores, INPUT_SHA
from exact_models import INPUTS, tighter_weights, simplex_linear_bounds

DATA = INPUTS['seikh']
CENTER = [F(n, DATA['denominator']) for n in DATA['numerators']]


def evaluate(box, q, rival, *, t=None):
    """Containing scores and all forty analytical partial derivatives of gap."""
    x = [IV.frac(a,b) for a,b in box]
    p = IV.frac(F(4)/q) if t is None else 4*IV.frac(*t)
    logs = [a.log() for a in x]
    z = [a.pow_int(4) for a in x]
    factors = []
    for j in range(5):
        e = 1 + sum((entropy_term(p * logs[(i*5+j)*2+c])
                     for i in range(4) for c in range(2)), IV.integer(0)).div_int(4)
        factors.append(IV(max(e.lo,SCALE//4),min(e.hi,SCALE)))
    weights = tighter_weights(factors)
    total = sum(factors)
    lp = [[(1-z[(i*5+j)*2]).log() for j in range(5)] for i in range(4)]
    ln = [[z[(i*5+j)*2+1].log() for j in range(5)] for i in range(4)]
    pp = [simplex_linear_bounds(weights,row).exp() for row in lp]
    pn = [simplex_linear_bounds(weights,row).exp() for row in ln]
    sc = [1-a-b for a,b in zip(pp,pn)]
    gap = sc[1]-sc[rival]
    a = [-pp[1]*lp[1][j]-pn[1]*ln[1][j]
         +pp[rival]*lp[rival][j]+pn[rival]*ln[rival][j] for j in range(5)]
    aw = simplex_linear_bounds(weights,a)
    grad = []
    for idx,xi in enumerate(x):
        i,j,c = idx//10,(idx%10)//2,idx%2
        de = (p*((p-1)*logs[idx]).exp()*(p*logs[idx]+1)).div_int(4)
        g = de*(a[j]-aw)/total
        sign = 1 if i==1 else -1 if i==rival else 0
        if sign and c==0:
            g = g + sign*4*xi.pow_int(3)*weights[j]*pp[i]/(1-z[idx])
        elif sign:
            g = g - sign*4*weights[j]*pn[i]/xi
        grad.append(g)
    return gap,grad


def point_gap(point,q,rival,*,t=None):
    ti=IV.frac(F(1)/F(q)) if t is None else IV.frac(*t)
    sc,_=scores(ti,[IV.frac(x**4) for x in point])
    return sc[1]-sc[rival]


def reduced_node(box,q,rival,*,t=None):
    """Replace globally monotone coordinates by minimizing endpoints.

    Each replacement is justified on the pre-replacement box. Sequential or
    simultaneous replacement leaves the box minimum unchanged. At termination,
    the centered enclosure follows the integral mean-value identity on a box.
    """
    box=list(box)
    reductions=[]
    while True:
        gap,gradient=evaluate(box,q,rival,t=t)
        moves=[]
        for k,((lo,hi),g) in enumerate(zip(box,gradient)):
            if lo==hi: continue
            if g.lo>=0: moves.append((k,lo,g))
            elif g.hi<=0: moves.append((k,hi,g))
        if not moves: break
        for k,value,g in moves:
            reductions.append({'coordinate':k,'endpoint':str(value),'gradient':g.endpoints()})
            box[k]=(value,value)
    midpoint=[(a+b)/2 for a,b in box]
    middle=point_gap(midpoint,q,rival,t=t)
    mv=middle
    for (a,b),m,g in zip(box,midpoint,gradient):
        mv=mv+g*IV.frac(a-m,b-m)
    lower=max(gap.lo,mv.lo)
    return box,gradient,lower,{'reductions':reductions,
                             'remaining_coordinates':sum(a!=b for a,b in box),
                             'gap_natural':gap.endpoints(),'gap_centered':mv.endpoints(),
                             'lower':str(F(lower,SCALE))}


def certify(radius,q,max_nodes=2047):
    radius,q=F(radius),F(q)
    if radius<=0 or q<4: raise ValueError('positive radius and q>=4 required')
    root=[(x-radius,x+radius) for x in CENTER]
    if any(a<=0 or b>=1 for a,b in root): raise ValueError('noninterior root')
    if any(root[i][1]**4+root[i+1][1]**4>1 for i in range(0,40,2)):
        raise ValueError('inadmissible source root')
    cases=[]
    for rival in (0,2,3):
        stack=[('r',root)]
        nodes=[]
        unresolved=[]
        while stack:
            ident,box=stack.pop()
            reduced,gradient,lower,rec=reduced_node(box,q,rival)
            rec.update(id=ident)
            if lower>0:
                rec['result']='POSITIVE'
            elif len(nodes)+len(stack)+2>max_nodes:
                rec['result']='UNRESOLVED_BUDGET'
                unresolved.append(ident)
            else:
                candidates=[((b-a)*(g.hi-g.lo),k) for k,((a,b),g) in enumerate(zip(reduced,gradient)) if a<b]
                if not candidates:
                    rec['result']='NONPOSITIVE_POINT';unresolved.append(ident)
                else:
                    _,k=max(candidates)
                    a,b=reduced[k];m=(a+b)/2
                    left=list(reduced);right=list(reduced)
                    left[k]=(a,m);right[k]=(m,b)
                    stack.extend([(ident+'1',right),(ident+'0',left)])
                    rec.update(result='SPLIT',coordinate=k,midpoint=str(m))
            nodes.append(rec)
        cases.append({'rival':rival,'node_count':len(nodes),'unresolved':unresolved,'nodes':nodes})
    return {'status':'PASS_COMPUTATIONAL' if all(not c['unresolved'] for c in cases) else 'UNRESOLVED',
            'model':'Seikh-Mandal entropy + product qROFAWA, unchanged V18',
            'model_input_sha256':INPUT_SHA,'bits':PREC,'q':str(q),'radius':str(radius),
            'metric':'40 normalized baseline-raw components, L-infinity, baseline rung 4',
            'method':'integer interval analytical derivatives, monotone-face reduction, centered form, bisection',
            'claim_limit':'one fixed rung; no all-q monotonicity, no exact optimum, no novelty claim for interval methods',
            'cases':cases}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--radius',required=True)
    p.add_argument('--q',required=True)
    p.add_argument('--max-nodes',type=int,default=255)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    start=time.perf_counter();result=certify(a.radius,a.q,a.max_nodes)
    result['elapsed_seconds']=time.perf_counter()-start
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='cases'}))
    print(json.dumps([{'rival':c['rival'],'nodes':c['node_count'],'unresolved':len(c['unresolved']),
                      'first_remaining':c['nodes'][0]['remaining_coordinates']} for c in result['cases']]))


if __name__=='__main__':main()
