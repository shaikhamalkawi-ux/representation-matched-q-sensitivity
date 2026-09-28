"""Integer replay checker for uniform joint t/source-face certificates.

Reconstructs every source box and t interval from the root. All reductions,
derivative bounds, natural/centered enclosures and positive margins are freshly
computed. T_SPLIT and X_SPLIT children cover the reconstructed parent face.
No proof decisions are accepted on the basis of stored positive/status flags.
"""
import argparse
import copy
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import importlib.util
import json
from pathlib import Path
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'refine_source_radius.py'
SOURCE_HASH=hashlib.sha256(SOURCE.read_bytes()).hexdigest()
spec=importlib.util.spec_from_file_location('joint_kernel_under_audit',SOURCE)
solver=importlib.util.module_from_spec(spec);spec.loader.exec_module(solver)

@lru_cache(maxsize=None)
def recompute(box,t,rival):
    # Cache only computations, never certificate values or acceptance verdicts.
    return solver.reduced_node(list(box),None,rival,t=t)

def check(record,progress=False):
    if record.get('status')!='PASS_COMPUTATIONAL':raise ValueError('status')
    if record.get('model_input_sha256')!=solver.INPUT_SHA:raise ValueError('model input hash')
    if record.get('bits')!=solver.PREC:raise ValueError('precision')
    if record.get('root_t')!=['0','1/4']:raise ValueError('root t coverage')
    if record.get('metric')!='40 normalized baseline-raw components, L-infinity, baseline rung 4':raise ValueError('metric')
    h=F(record['radius'])
    if h<=0:raise ValueError('radius')
    root=tuple((x-h,x+h) for x in solver.CENTER)
    if any(a<=0 or b>=1 for a,b in root):raise ValueError('source interior')
    if any(root[k][1]**4+root[k+1][1]**4>1 for k in range(0,40,2)):raise ValueError('source admissibility')
    if [c['rival'] for c in record['cases']]!=[0,2,3]:raise ValueError('rival coverage')
    total=0;stats=[]
    for case in record['cases']:
        if case['unresolved']:raise ValueError('unresolved rival')
        if case['node_count']!=len(case['nodes']):raise ValueError('node count')
        nodes={node['id']:node for node in case['nodes']}
        if len(nodes)!=len(case['nodes']):raise ValueError('duplicate identifier')
        pending=[('r',root,(F(0),F(1,4)))];seen=set()
        counts={'POSITIVE':0,'T_SPLIT':0,'X_SPLIT':0};leaf_lowers=[];reductions=0
        while pending:
            ident,box,t=pending.pop()
            if ident in seen or ident not in nodes:raise ValueError('missing/duplicate branch '+ident)
            seen.add(ident);node=nodes[ident]
            if node.get('t')!=[str(a) for a in t]:raise ValueError('node t interval '+ident)
            face,gradient,lower,expected=recompute(box,t,case['rival'])
            for key in ('reductions','remaining_coordinates','gap_natural','gap_centered','lower'):
                if node.get(key)!=expected[key]:raise ValueError('recomputed '+key+' '+ident)
            reductions+=len(expected['reductions'])
            result=node.get('result')
            if result not in counts:raise ValueError('unknown/unresolved result '+ident)
            counts[result]+=1
            if result=='POSITIVE':
                if lower<=0:raise ValueError('forged positive leaf '+ident)
                leaf_lowers.append(F(lower,solver.SCALE))
                if 'coordinate' in node or 'midpoint' in node:raise ValueError('leaf split metadata')
            elif result=='T_SPLIT':
                a,b=t;mid=F(node['midpoint'])
                if a>=b or mid!=(a+b)/2:raise ValueError('T split midpoint')
                if 'coordinate' in node:raise ValueError('T split source coordinate')
                pending.extend([(ident+'1',tuple(face),(mid,b)),(ident+'0',tuple(face),(a,mid))])
            else:
                k=node['coordinate']
                if type(k) is not int or not 0<=k<40:raise ValueError('X split coordinate')
                a,b=face[k];mid=F(node['midpoint'])
                if a>=b or mid!=(a+b)/2:raise ValueError('X split midpoint')
                left=list(face);right=list(face);left[k]=(a,mid);right[k]=(mid,b)
                pending.extend([(ident+'1',tuple(right),t),(ident+'0',tuple(left),t)])
            if progress and len(seen)%100==0:
                print(json.dumps({'rival':case['rival'],'checked':len(seen),'total_case_nodes':len(nodes)}),flush=True)
        if set(nodes)!=seen:raise ValueError('unreachable/extra nodes')
        if not leaf_lowers:raise ValueError('no positive leaves')
        total+=len(seen)
        stats.append({'rival':case['rival'],'nodes':len(seen),'decisions':counts,'endpoint_reductions_recomputed':reductions,
                      'minimum_positive_lower':str(min(leaf_lowers))})
    return {'nodes':total,'cases':stats}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input',default='UNIFORM_JOINT_H_0_0232195.json')
    parser.add_argument('--output',default='UNIFORM_JOINT_INTEGER_AUDIT.json');args=parser.parse_args()
    path=ROOT/args.input;out=HERE/args.output
    if out.exists():raise FileExistsError(out)
    start=time.perf_counter();record=json.loads(path.read_text());stats=check(record,progress=True)
    mutations=[]
    def first(kind,obj):
        return next(n for c in obj['cases'] for n in c['nodes'] if n['result']==kind)
    def reject(label,fn):
        x=copy.deepcopy(record);fn(x)
        try:check(x)
        except (ValueError,KeyError,IndexError,TypeError) as e:
            mutations.append({'mutation':label,'rejected':True,'reason':str(e)});return
        raise AssertionError('accepted malformed certificate: '+label)
    reject('T to X split axis',lambda x:first('T_SPLIT',x).__setitem__('result','X_SPLIT'))
    reject('X to T split axis',lambda x:first('X_SPLIT',x).__setitem__('result','T_SPLIT'))
    reject('T split midpoint',lambda x:first('T_SPLIT',x).__setitem__('midpoint','1/7'))
    reject('X split midpoint',lambda x:first('X_SPLIT',x).__setitem__('midpoint','1/2'))
    reject('X split coordinate',lambda x:first('X_SPLIT',x).__setitem__('coordinate',40))
    reject('source minimizing endpoint',lambda x:x['cases'][0]['nodes'][0]['reductions'][0].__setitem__('endpoint','1/2'))
    reject('source derivative bound',lambda x:x['cases'][0]['nodes'][0]['reductions'][0].__setitem__('gradient',['0','0']))
    reject('root t interval',lambda x:x.__setitem__('root_t',['0','1/8']))
    reject('child t interval',lambda x:x['cases'][0]['nodes'][1].__setitem__('t',['0','1/16']))
    def omit(x):
        x['cases'][0]['nodes'].pop(1);x['cases'][0]['node_count']-=1
    reject('omit branch and adjust count',omit)
    reject('forge positive at root',lambda x:x['cases'][0]['nodes'][0].__setitem__('result','POSITIVE'))
    reject('forge positive lower',lambda x:first('POSITIVE',x).__setitem__('lower','1'))
    reject('omit rival',lambda x:x['cases'].pop())
    reject('duplicate rival',lambda x:x['cases'].__setitem__(1,copy.deepcopy(x['cases'][0])))
    reject('wrong precision',lambda x:x.__setitem__('bits',53))
    reject('wrong model hash',lambda x:x.__setitem__('model_input_sha256','0'*64))
    result={'status':'PASS_UNIFORM_JOINT_INTEGER_REPLAY','input':path.name,'input_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
      'solver_sha256_at_load':SOURCE_HASH,'solver_sha256_at_finish':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
      'radius':record['radius'],'root_t':record['root_t'],**stats,'negative_mutations':mutations,
      'fresh_computations':recompute.cache_info().misses,'elapsed_seconds':time.perf_counter()-start,
      'scope':'Every finite real q>=4 and complete normalized baseline source cube. Uses original integer enclosure kernel; analytic derivatives audited separately.'}
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('negative_mutations','cases')},indent=2))

if __name__=='__main__':main()
