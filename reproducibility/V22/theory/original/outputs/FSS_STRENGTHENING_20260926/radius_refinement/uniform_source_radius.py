"""Uniform t=1/q cover using the fixed source-box reduction.

Positive leaves prove the complete source box for the stated t interval.
Source boxes reduced for a parent are valid minimizing faces for every child.
The t=0 endpoint is the continuous positive-input limit, not a finite rung.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import time
from refine_source_radius import CENTER,INPUT_SHA,PREC,reduced_node


def certify_uniform(radius,max_nodes=4095,max_depth=28,joint=False):
    h=F(radius)
    root=[(x-h,x+h) for x in CENTER]
    if h<=0 or any(a<=0 or b>=1 for a,b in root):raise ValueError('noninterior root')
    if any(root[k][1]**4+root[k+1][1]**4>1 for k in range(0,40,2)):
        raise ValueError('inadmissible source root')
    cases=[]
    for rival in (0,2,3):
        stack=[('r',F(0),F(1,4),root,0)]
        nodes=[];unresolved=[]
        while stack:
            ident,a,b,box,depth=stack.pop()
            face,grad,lower,rec=reduced_node(box,None,rival,t=(a,b))
            rec.update(id=ident,t=[str(a),str(b)])
            if lower>0:
                rec['result']='POSITIVE'
            elif len(nodes)+len(stack)+2>max_nodes or depth>=max_depth:
                rec['result']='UNRESOLVED_BUDGET';unresolved.append(ident)
            else:
                m=(a+b)/2
                # This is only a subdivision heuristic, never a proof test.
                candidates=[((hi-lo)*(g.hi-g.lo),k) for k,((lo,hi),g) in enumerate(zip(face,grad)) if lo<hi]
                split_x=False
                if joint and candidates and b-a<=F(1,256):
                    _,_,at_mid,_=reduced_node(face,None,rival,t=(m,m))
                    split_x=at_mid<=0
                if split_x:
                    _,k=max(candidates);lo,hi=face[k];cut=(lo+hi)/2
                    left=list(face);right=list(face);left[k]=(lo,cut);right[k]=(cut,hi)
                    rec.update(result='X_SPLIT',coordinate=k,midpoint=str(cut))
                    stack.extend([(ident+'1',a,b,right,depth+1),(ident+'0',a,b,left,depth+1)])
                else:
                    rec.update(result='T_SPLIT',midpoint=str(m))
                    stack.extend([(ident+'1',m,b,face,depth+1),(ident+'0',a,m,face,depth+1)])
            nodes.append(rec)
        cases.append({'rival':rival,'node_count':len(nodes),'unresolved':unresolved,'nodes':nodes})
    return {'status':'PASS_COMPUTATIONAL' if all(not c['unresolved'] for c in cases) else 'UNRESOLVED',
            'model':'unchanged Seikh-Mandal V18','model_input_sha256':INPUT_SHA,'bits':PREC,
            'radius':str(h),'root_t':['0','1/4'],'scope':'every finite real q>=4; t=0 continuous limit',
            'metric':'40 normalized baseline-raw components, L-infinity, baseline rung 4',
            'max_nodes':max_nodes,'max_depth':max_depth,'joint_splits':joint,'cases':cases}


def main():
    p=argparse.ArgumentParser();p.add_argument('--radius',required=True)
    p.add_argument('--max-nodes',type=int,default=4095);p.add_argument('--max-depth',type=int,default=28)
    p.add_argument('--joint',action='store_true')
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    start=time.perf_counter();r=certify_uniform(a.radius,a.max_nodes,a.max_depth,a.joint)
    r['elapsed_seconds']=time.perf_counter()-start
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in r.items() if k!='cases'}))
    print(json.dumps([{'rival':c['rival'],'nodes':c['node_count'],'unresolved':len(c['unresolved'])} for c in r['cases']]))


if __name__=='__main__':main()
