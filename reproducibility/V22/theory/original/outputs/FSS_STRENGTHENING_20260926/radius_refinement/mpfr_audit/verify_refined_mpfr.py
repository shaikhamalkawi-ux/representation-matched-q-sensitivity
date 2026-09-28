"""Independent MPFR256 arithmetic audit of monotone-face Seikh records.

No import from the integer generator or its interval/model evaluators. The model
input is hash-pinned; bounded simplex extrema use exact rational vertex
enumeration, not the generator's dyadic greedy optimizer.
"""
from __future__ import annotations
import argparse
import copy
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import sys
import time

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'deps'))
import gmpy2 as g

BITS=256
g.get_context().precision=BITS
DOWN=g.context(precision=BITS,round=g.RoundDown)
UP=g.context(precision=BITS,round=g.RoundUp)
INPUT_SHA='71351416bb73b22f3781bad962a023b6ea7362e8020429aab97b47d3f36d1797'
INPUT_PATH=HERE/'model_inputs.json'

def require(v,msg):
    if not v: raise ValueError(msg)

def rat(x):
    if isinstance(x,F): return x
    if isinstance(x,g.mpfr):
        a,b=x.as_integer_ratio(); return F(int(a),int(b))
    return F(x)

def converted(x,ctx):
    x=rat(x)
    with ctx: return g.mpfr(g.mpq(x.numerator,x.denominator))

class I:
    __slots__=('lo','hi')
    def __init__(self,lo,hi=None):
        if hi is None: hi=lo
        self.lo=converted(lo,DOWN);self.hi=converted(hi,UP)
        require(self.lo<=self.hi,'reversed interval')
    @classmethod
    def direct(cls,lo,hi):
        obj=object.__new__(cls);obj.lo=lo;obj.hi=hi
        require(g.is_finite(lo) and g.is_finite(hi) and lo<=hi,'invalid interval')
        return obj
    def __add__(self,b):
        b=asiv(b)
        with DOWN: lo=self.lo+b.lo
        with UP: hi=self.hi+b.hi
        return I.direct(lo,hi)
    __radd__=__add__
    def __neg__(self):
        with DOWN:lo=-self.hi
        with UP:hi=-self.lo
        return I.direct(lo,hi)
    def __sub__(self,b): return self+-asiv(b)
    def __rsub__(self,b): return asiv(b)+-self
    def __mul__(self,b):
        b=asiv(b)
        with DOWN: lo=min(a*c for a in (self.lo,self.hi) for c in (b.lo,b.hi))
        with UP: hi=max(a*c for a in (self.lo,self.hi) for c in (b.lo,b.hi))
        return I.direct(lo,hi)
    __rmul__=__mul__
    def __truediv__(self,b):
        b=asiv(b); require(not b.lo<=0<=b.hi,'division through zero')
        with DOWN: lo=1/b.hi
        with UP: hi=1/b.lo
        return self*I.direct(lo,hi)
    def __rtruediv__(self,b):return asiv(b)/self
    def power(self,n):
        require(isinstance(n,int) and n>=0 and self.lo>=0,'nonnegative integer power domain')
        with DOWN:lo=self.lo**n
        with UP:hi=self.hi**n
        return I.direct(lo,hi)
    def log(self):
        require(self.lo>0,'log domain')
        with DOWN:lo=g.log(self.lo)
        with UP:hi=g.log(self.hi)
        return I.direct(lo,hi)
    def exp(self):
        with DOWN:lo=g.exp(self.lo)
        with UP:hi=g.exp(self.hi)
        return I.direct(lo,hi)
    def pair(self):return [str(rat(self.lo)),str(rat(self.hi))]

def asiv(x):return x if isinstance(x,I) else I(x)

def yexp(y):
    require(y.hi<=0,'positive entropy argument')
    candidates=[I.direct(v,v)*I.direct(v,v).exp() for v in (y.lo,y.hi)]
    lo=min(v.lo for v in candidates);hi=max(v.hi for v in candidates)
    if y.lo<=-1<=y.hi:lo=min(lo,(-I(-1).exp()).lo)
    return I.direct(lo,hi)

def weight_box(e):
    out=[]
    for j in range(5):
        l=I.direct(e[j].lo,e[j].lo)/(I.direct(e[j].lo,e[j].lo)+sum(I.direct(e[k].hi,e[k].hi) for k in range(5) if k!=j))
        u=I.direct(e[j].hi,e[j].hi)/(I.direct(e[j].hi,e[j].hi)+sum(I.direct(e[k].lo,e[k].lo) for k in range(5) if k!=j))
        out.append(I.direct(l.lo,u.hi))
    return out

def vertices(weights):
    bounds=[(rat(x.lo),rat(x.hi)) for x in weights]
    require(all(0<a<=b<1 for a,b in bounds),'weight bound domain')
    result=set()
    for free in range(5):
        rest=[j for j in range(5) if j!=free]
        for bits in itertools.product((0,1),repeat=4):
            w=[F(0)]*5
            for j,b in zip(rest,bits):w[j]=bounds[j][b]
            w[free]=1-sum(w)
            if bounds[free][0]<=w[free]<=bounds[free][1]:result.add(tuple(w))
    require(result,'empty simplex')
    return tuple(result)

def lp(vs,coeff):
    lo=[rat(x.lo) for x in coeff];hi=[rat(x.hi) for x in coeff]
    return I(min(sum(w*a for w,a in zip(v,lo)) for v in vs),
             max(sum(w*a for w,a in zip(v,hi)) for v in vs))

def evaluate(box,t,rival,gradients=True):
    x=[I(a,b) for a,b in box];p=4*I(*t)
    logs=[a.log() for a in x];z=[a.power(4) for a in x]
    e=[]
    for j in range(5):
        a=1+sum(yexp(p*logs[(i*5+j)*2+c]) for i in range(4) for c in range(2))/4
        e.append(I.direct(max(a.lo,I(F(1,4)).lo),min(a.hi,I(1).hi)))
    weights=weight_box(e);vs=vertices(weights)
    lp_rows=[[(1-z[(i*5+j)*2]).log() for j in range(5)] for i in range(4)]
    ln_rows=[[z[(i*5+j)*2+1].log() for j in range(5)] for i in range(4)]
    pp=[lp(vs,row).exp() for row in lp_rows];pn=[lp(vs,row).exp() for row in ln_rows]
    scores=[1-a-b for a,b in zip(pp,pn)];gap=scores[1]-scores[rival]
    if not gradients:return gap
    # derivative of gap with respect to each weight before normalization.
    dw=[-pp[1]*lp_rows[1][j]-pn[1]*ln_rows[1][j]
        +pp[rival]*lp_rows[rival][j]+pn[rival]*ln_rows[rival][j] for j in range(5)]
    common=lp(vs,dw);total=sum(e);derivatives=[]
    for k,a in enumerate(x):
        i,j,c=k//10,(k%10)//2,k%2
        de=p*((p-1)*logs[k]).exp()*(1+p*logs[k])/4
        partial=(dw[j]-common)*de/total
        sign=int(i==1)-int(i==rival)
        if sign:
            if c==0:partial+=sign*4*weights[j]*pp[i]*a.power(3)/(1-z[k])
            else:partial-=sign*4*weights[j]*pn[i]/a
        derivatives.append(partial)
    return gap,derivatives

def pair(value,label):
    require(isinstance(value,list) and len(value)==2,label+' not a pair')
    a,b=map(F,value);require(a<=b,label+' reversed');return a,b

def encloses(stored,computed,label):
    a,b=pair(stored,label)
    require(a<=rat(computed.lo) and rat(computed.hi)<=b,label+' fails MPFR enclosure')

def check_node(node,source_box,t,rival):
    box=list(source_box);reductions=node['reductions']
    require(isinstance(reductions,list),'reductions not a list')
    cursor=0;passes=0
    # Verify each stored reduction on a containing current box. Cached gradients
    # can justify several moves simultaneously; recompute when the old box's
    # bound cannot justify the next recorded reduction. Never trust its sign.
    while cursor<len(reductions):
        _,gradient=evaluate(box,t,rival);passes+=1;initial=cursor
        while cursor<len(reductions):
            red=reductions[cursor];k=red['coordinate']
            require(type(k) is int and 0<=k<40,'bad reduced coordinate')
            a,b=box[k];require(a<b,'duplicate or already collapsed reduction')
            v=F(red['endpoint']);require(v in (a,b),'reduction not at a face')
            gk=gradient[k]
            # A previously computed containing box stays valid after reductions.
            okay=(v==a and gk.lo>=0) or (v==b and gk.hi<=0)
            if not okay:break
            try:encloses(red['gradient'],gk,'stored gradient')
            except ValueError:break
            box[k]=(v,v);cursor+=1
        require(cursor>initial,'unproved face at '+node['id']+' reduction '+str(cursor))
    natural,gradient=evaluate(box,t,rival)
    midpoint=[(a+b)/2 for a,b in box]
    centered=evaluate([(x,x) for x in midpoint],t,rival,False)
    for (a,b),m,d in zip(box,midpoint,gradient):centered+=d*I(a-m,b-m)
    require(node['remaining_coordinates']==sum(a<b for a,b in box),'remaining count wrong')
    encloses(node['gap_natural'],natural,'natural gap')
    encloses(node['gap_centered'],centered,'centered gap')
    expected=max(F(node['gap_natural'][0]),F(node['gap_centered'][0]))
    require(F(node['lower'])==expected,'stored lower inconsistent')
    lower=max(rat(natural.lo),rat(centered.lo))
    if node['result']=='POSITIVE':require(lower>0 and expected>0,'nonpositive leaf')
    return box,lower,{'id':node['id'],'passes':passes,'reductions':len(reductions),
                     'lower_mpfr':str(lower),'natural_mpfr':natural.pair(),'centered_mpfr':centered.pair()}

def input_center():
    raw=INPUT_PATH.read_bytes();require(hashlib.sha256(raw).hexdigest()==INPUT_SHA,'input hash')
    d=json.loads(raw)['seikh'];require(d['r']==4 and d['shape']==[4,5,2],'model geometry')
    return [F(n,d['denominator']) for n in d['numerators']]

def verify(record):
    require(record['status']=='PASS_COMPUTATIONAL','incomplete record')
    require(record['model_input_sha256']==INPUT_SHA,'record model pin')
    require(record['bits']==224,'unexpected source record precision')
    require(record['metric']=='40 normalized baseline-raw components, L-infinity, baseline rung 4','wrong source metric')
    h=F(record['radius']);require(h>0,'nonpositive radius')
    root=[(a-h,a+h) for a in input_center()]
    require(all(0<a<=b<1 for a,b in root),'noninterior root')
    require(all(root[k][1]**4+root[k+1][1]**4<=1 for k in range(0,40,2)),'inadmissible root')
    uniform='root_t' in record
    if uniform:require(pair(record['root_t'],'root t')==(F(0),F(1,4)),'uniform root')
    else:require(F(record['q'])>=4,'q domain')
    root_t=(F(0),F(1,4)) if uniform else (1/F(record['q']),)*2
    require([c['rival'] for c in record['cases']]==[0,2,3],'rival set/order')
    output=[]
    for case in record['cases']:
        require(case['unresolved']==[],'unresolved case')
        nodes={n['id']:n for n in case['nodes']}
        require(len(nodes)==len(case['nodes'])==case['node_count'],'duplicate/count')
        stack=[('r',root,root_t)];seen=set();rows=[];leaves=0
        while stack:
            ident,box,t=stack.pop();require(ident in nodes and ident not in seen,'missing node')
            n=nodes[ident];seen.add(ident)
            if uniform:require(pair(n['t'],'node t')==t,'t cover mismatch')
            face,lower,row=check_node(n,box,t,case['rival']);rows.append(row)
            if n['result']=='POSITIVE':leaves+=1
            elif uniform and n['result']=='T_SPLIT':
                require(F(n['lower'])<=0,'split with stored positive lower')
                a,b=t;m=(a+b)/2;require(a<m<b and F(n['midpoint'])==m,'t midpoint')
                stack.extend([(ident+'1',face,(m,b)),(ident+'0',face,(a,m))])
            else:
                require(F(n['lower'])<=0,'split with stored positive lower')
                require(n['result']==('X_SPLIT' if uniform else 'SPLIT'),'source split kind')
                k=n['coordinate'];require(type(k)is int and 0<=k<40,'split coordinate')
                a,b=face[k];m=(a+b)/2;require(a<m<b and F(n['midpoint'])==m,'x midpoint')
                left=list(face);right=list(face);left[k]=(a,m);right[k]=(m,b)
                stack.extend([(ident+'1',right,t),(ident+'0',left,t)])
        require(seen==set(nodes) and len(nodes)==2*leaves-1,'unreachable or incomplete tree')
        output.append({'rival':case['rival'],'node_count':len(rows),'leaf_count':leaves,
                       'leaf_floor':str(min(F(r['lower_mpfr']) for r in rows if nodes[r['id']]['result']=='POSITIVE')),
                       'checks':rows})
    return {'status':'PASS_MPFR256','mpfr_version':g.mpfr_version(),'gmpy2_version':g.version(),
            'bits':BITS,'radius':str(h),'root_t':[str(x) for x in root_t],
            'model_input_sha256':INPUT_SHA,'cases':output,
            'scope':'Complete tree, stored interval consistency, all gradient face reductions, centered form and every positive leaf independently recomputed. Shared model assumptions; not human peer review.'}

def main():
    p=argparse.ArgumentParser();p.add_argument('records',type=Path,nargs='+');a=p.parse_args()
    for path in a.records:
        start=time.perf_counter();result=verify(json.loads(path.read_text()))
        result['record_sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
        result['checker_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result['elapsed_seconds']=time.perf_counter()-start
        out=HERE/(path.stem+'_MPFR256.json')
        out.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        print(json.dumps({'record':path.name,'status':result['status'],'seconds':result['elapsed_seconds'],
                          'cases':[{k:c[k] for k in ('rival','node_count','leaf_count','leaf_floor')} for c in result['cases']]}),flush=True)

if __name__=='__main__':main()
