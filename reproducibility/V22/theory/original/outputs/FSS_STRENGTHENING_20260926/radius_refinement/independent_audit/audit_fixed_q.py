"""Read-only independent fixed-q analytical/numerical audit and strict replay.

The mpmath oracle uses direct sums/products, not interval LP score bounds.
mp.diff differentiates that oracle independently of the solver's gradient.
Sampling is a diagnostic, not an interval proof. Proof validity is audited
separately in PROOF_AUDIT.md, and retained decisions are replayed exactly.
"""
from fractions import Fraction as F
from pathlib import Path
import copy
import argparse
import hashlib
import importlib.util
import json
import random
import sys
import time
sys.path.insert(0,str(Path(__file__).resolve().parent/'private_dependencies'))
import mpmath as mp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
source=ROOT/'refine_source_radius.py'
source_hash=hashlib.sha256(source.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('radius_solver_under_audit',source)
solver=importlib.util.module_from_spec(spec);spec.loader.exec_module(solver)
mp.mp.dps=100

def m(x):
    x=F(x)
    return mp.mpf(x.numerator)/x.denominator

def oracle(x,q,rival):
    p=4/m(q)
    ef=[]
    for j in range(5):
        ef.append(1+sum(x[(i*5+j)*2+c]**p*mp.log(x[(i*5+j)*2+c]**p)
                        for i in range(4) for c in range(2))/4)
    total=sum(ef);w=[e/total for e in ef]
    out=[]
    for i in range(4):
        pp=mp.exp(sum(w[j]*mp.log(1-x[(i*5+j)*2]**4) for j in range(5)))
        pn=mp.exp(sum(w[j]*mp.log(x[(i*5+j)*2+1]**4) for j in range(5)))
        out.append(1-pp-pn)
    return out[1]-out[rival]

def derivatives(x,q,rival):
    answer=[]
    for k in range(40):
        def f(v):
            z=list(x);z[k]=v
            return oracle(z,q,rival)
        answer.append(mp.diff(f,x[k]))
    return answer

def contains(iv,value):
    return m(F(iv.lo,solver.SCALE)) <= value <= m(F(iv.hi,solver.SCALE))

def replay(record):
    """Freshly recompute every proof node, require exact tree and all margins.

    Unlike a trusting status-reader this rejects altered gradients/endpoints,
    deleted rival/tree branches, forged lower bounds, and false PASS flags.
    This replay shares the integer evaluator; the oracle audit does not.
    """
    if record.get('status')!='PASS_COMPUTATIONAL':raise ValueError('not PASS')
    if record.get('model_input_sha256')!=solver.INPUT_SHA:raise ValueError('input hash')
    if record.get('bits')!=solver.PREC:raise ValueError('precision')
    q,radius=F(record['q']),F(record['radius'])
    if q<4 or radius<=0:raise ValueError('domain')
    root=[(x-radius,x+radius) for x in solver.CENTER]
    if any(a<=0 or b>=1 for a,b in root):raise ValueError('root boundary')
    if any(root[i][1]**4+root[i+1][1]**4>1 for i in range(0,40,2)):raise ValueError('root admissibility')
    if [c['rival'] for c in record['cases']]!=[0,2,3]:raise ValueError('rival coverage')
    nn=0
    for case in record['cases']:
        if case['unresolved']:raise ValueError('unresolved')
        if case['node_count']!=len(case['nodes']):raise ValueError('node_count')
        nodes={n['id']:n for n in case['nodes']}
        if len(nodes)!=len(case['nodes']):raise ValueError('duplicate node')
        pending=[('r',root)];seen=set()
        while pending:
            ident,box=pending.pop()
            if ident not in nodes or ident in seen:raise ValueError('tree gap/duplicate')
            seen.add(ident);node=nodes[ident]
            reduced,grad,lower,expected=solver.reduced_node(box,q,case['rival'])
            for key in ('reductions','remaining_coordinates','gap_natural','gap_centered','lower'):
                if node.get(key)!=expected[key]:raise ValueError('recomputed '+key)
            if node['result']=='POSITIVE':
                if lower<=0:raise ValueError('false positive')
            elif node['result']=='SPLIT':
                if lower>0:raise ValueError('unexpected split')
                k=node['coordinate'];a,b=reduced[k];mid=F(node['midpoint'])
                if a==b or mid!=(a+b)/2:raise ValueError('split midpoint')
                left=list(reduced);right=list(reduced);left[k]=(a,mid);right[k]=(mid,b)
                pending.extend([(ident+'1',right),(ident+'0',left)])
            else:raise ValueError('unresolved or unknown node')
        if seen!=set(nodes):raise ValueError('unused nodes')
        nn+=len(seen)
    return nn

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',default='FIXED_Q_AUDIT.json')
    args=parser.parse_args()
    target=HERE/args.output
    if target.exists():raise FileExistsError(target)
    start=time.perf_counter();rng=random.Random(20260926)
    tests=[];grad_count=0;point_grad_count=0;gap_count=0
    cases=[(4,F('0.0232195')),(8,F('0.0250515')),(16,F('0.0258973'))]
    # Root-box extrema, centre and deterministic interior test points.
    for q,radius in cases:
        root=[(x-radius,x+radius) for x in solver.CENTER]
        for rival in (0,2,3):
            interval,gradient=solver.evaluate(root,F(q),rival)
            points=[list(solver.CENTER),[a for a,b in root],[b for a,b in root]]
            points.append([a+(b-a)*F(rng.randrange(10001),10000) for a,b in root])
            for pi,point in enumerate(points):
                xx=[m(x) for x in point]
                value=oracle(xx,q,rival)
                if not contains(interval,value):raise AssertionError(('root gap',q,rival,pi))
                gap_count+=1
                dd=derivatives(xx,q,rival)
                for k,d in enumerate(dd):
                    if not contains(gradient[k],d):raise AssertionError(('root gradient',q,rival,pi,k))
                    grad_count+=1
                # Point intervals are extremely narrow; numerical differences
                # with 100 decimal precision therefore test formula equality.
                if pi==0:
                    vg,gg=solver.evaluate([(v,v) for v in point],F(q),rival)
                    if not contains(vg,value):raise AssertionError(('point gap',q,rival))
                    for k,d in enumerate(dd):
                        if not contains(gg[k],d):raise AssertionError(('point derivative',q,rival,k,str(d),gg[k].endpoints()))
                        point_grad_count+=1
            tests.append({'q':q,'radius':str(radius),'rival':rival,'points':len(points),'status':'PASS'})
    certificates=[]
    for name in ('Q4_H_0_0232195.json','Q8_H_0_0250515.json','Q16_H_0_0258973.json'):
        obj=json.loads((ROOT/name).read_text());n=replay(obj)
        certificates.append({'name':name,'sha256':hashlib.sha256((ROOT/name).read_bytes()).hexdigest(),'nodes':n,'status':'PASS'})
    seed=json.loads((ROOT/'Q4_H_0_0232195.json').read_text())
    mutations=[]
    def bad(label,fn):
        x=copy.deepcopy(seed);fn(x)
        try:replay(x)
        except (ValueError,KeyError,IndexError,TypeError) as e:mutations.append({'mutation':label,'rejected':True,'reason':str(e)});return
        raise AssertionError('accepted mutation: '+label)
    bad('delete rival',lambda x:x['cases'].pop())
    bad('duplicate rival',lambda x:x['cases'].__setitem__(1,copy.deepcopy(x['cases'][0])))
    bad('delete node',lambda x:x['cases'][0]['nodes'].pop())
    bad('forge lower',lambda x:x['cases'][0]['nodes'][0].__setitem__('lower','1'))
    bad('alter reduction endpoint',lambda x:x['cases'][0]['nodes'][0]['reductions'][0].__setitem__('endpoint','1/2'))
    bad('forge gradient',lambda x:x['cases'][0]['nodes'][0]['reductions'][0].__setitem__('gradient',['1','2']))
    bad('alter radius',lambda x:x.__setitem__('radius','1/50'))
    bad('alter q',lambda x:x.__setitem__('q','8'))
    bad('alter model hash',lambda x:x.__setitem__('model_input_sha256','0'*64))
    bad('alter precision',lambda x:x.__setitem__('bits',53))
    bad('unresolved flag',lambda x:x['cases'][0].__setitem__('unresolved',['r']))
    bad('false status',lambda x:x.__setitem__('status','UNRESOLVED'))
    result={'status':'PASS_FIXED_Q_AUDIT','solver_sha256_at_load':source_hash,
      'solver_sha256_at_finish':hashlib.sha256(source.read_bytes()).hexdigest(),
      'oracle_precision_decimal':mp.mp.dps,'gradient_root_containment_checks':grad_count,
      'gradient_point_containment_checks':point_grad_count,'gap_root_containment_checks':gap_count,
      'cases':tests,'certificate_replay':certificates,'negative_mutations':mutations,
      'scope':'Only fixed q=4,8,16. Numerical sampling supports implementation checking but is not the containing proof.',
      'elapsed_seconds':time.perf_counter()-start}
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('cases','negative_mutations')},indent=2))

if __name__=='__main__':main()
