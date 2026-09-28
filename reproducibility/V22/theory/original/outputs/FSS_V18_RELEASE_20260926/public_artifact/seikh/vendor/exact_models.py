"""Independent integer-interval implementation of two frozen mathematical models.
Uses only dyadic_interval.py plus Python standard library; no inherited .so/code.
Model formulas share mathematical definitions with the MPFR implementation.
"""
from fractions import Fraction as F
from pathlib import Path
import json,itertools
from dyadic_interval import IV,SCALE
HERE=Path(__file__).resolve().parent
INPUTS=json.loads((HERE/'model_inputs.json').read_text())

def domain(family,kind='exact',radius=F(0)):
    d=INPUTS[family]; raw=[F(n,d['denominator']) for n in d['numerators']]
    if radius<0:raise ValueError('negative radius')
    if kind=='exact':
        if radius!=0:raise ValueError('nonzero radius with exact domain')
        exact=[(x**d['r'],x**d['r']) for x in raw]
    elif kind=='source':
        if any(x-radius<=0 or x+radius>=1 for x in raw):raise ValueError('source domain not strictly inside (0,1)')
        exact=[((x-radius)**d['r'],(x+radius)**d['r']) for x in raw]
    elif kind=='canonical':exact=[(x**d['r']-radius,x**d['r']+radius) for x in raw]
    else:raise ValueError('unknown domain convention')
    if any(a<=0 or a>b or b>=1 for a,b in exact):raise ValueError('not an interior box')
    if any(exact[i][1]+exact[i+1][1]>1 for i in range(0,len(exact),2)):raise ValueError('inadmissible source box')
    return [IV.frac(a,b) for a,b in exact]

def tighter_weights(e):
    """Exact endpoint monotonicity of e_j/sum e, not denominator sampling."""
    out=[]
    for j,x in enumerate(e):
        lowden=x.lo+sum(y.hi for k,y in enumerate(e) if k!=j)
        highden=x.hi+sum(y.lo for k,y in enumerate(e) if k!=j)
        a=IV(x.lo,x.lo)/IV(lowden,lowden)
        b=IV(x.hi,x.hi)/IV(highden,highden)
        out.append(IV(a.lo,b.hi))
    return out

def simplex_linear_bounds(weights,values):
    """Exact LP of sum w_j*a_j with interval coefficients and sum w=1.
    Relaxation is valid even when actual weights depend on actual coefficients.
    Returns fixed dyadic endpoints via exact integer greedy allocation.
    """
    if any(w.lo<0 for w in weights):raise ValueError('negative weight')
    if sum(w.lo for w in weights)>SCALE or sum(w.hi for w in weights)<SCALE:raise ValueError('empty weight simplex')
    def extremum(upper):
        coeff=[v.hi if upper else v.lo for v in values]
        w=[v.lo for v in weights];left=SCALE-sum(w)
        for j in sorted(range(len(w)),key=lambda j:coeff[j],reverse=upper):
            delta=min(left,weights[j].hi-w[j]);w[j]+=delta;left-=delta
        if left:raise ValueError('LP allocation failure')
        val=sum(a*b for a,b in zip(w,coeff))
        return -((-val)//SCALE) if upper else val//SCALE
    return IV(extremum(False),extremum(True))

def seikh(t,inputs,tight=False,ratio=None,simplex=None):
    ratio=tight if ratio is None else ratio
    simplex=tight if simplex is None else simplex
    logs=[x.log() for x in inputs]
    e=[]
    for j in range(5):
        total=IV.integer(0)
        for i in range(4):
            for c in range(2):
                y=t*logs[(i*5+j)*2+c]
                total=total+y*y.exp()
        ej=1+total.div_int(4)
        # e_j >= 1-2/e > 1/4 and <=1, for all raw entries in [0,1].
        ej=IV(max(ej.lo,SCALE//4),min(ej.hi,SCALE))
        e.append(ej)
    w=tighter_weights(e) if ratio else [x/sum(e) for x in e]
    scores=[]
    for i in range(4):
        if simplex:
            lp=[(1-inputs[(i*5+j)*2]).log() for j in range(5)]
            ln=[logs[(i*5+j)*2+1] for j in range(5)]
            pp=simplex_linear_bounds(w,lp).exp()
            pn=simplex_linear_bounds(w,ln).exp()
        else:
            pp=pn=IV.integer(1)
            for j in range(5):
                pp=pp*(1-inputs[(i*5+j)*2]).pow_unit(w[j])
                pn=pn*inputs[(i*5+j)*2+1].pow_unit(w[j])
        scores.append(1-pp-pn)
    return scores,w

def mono_map(x,code):
    if x.lo<0 or x.hi>SCALE:raise ValueError('bounded transform input outside [0,1]')
    def at(n):
        y=IV(n,n)
        if code=='complement_cross':return (1-y)/(1+8*y)
        if code=='mu':return (1-y)/(1+2*y)
        if code=='nu':return y/(3-2*y)
        if code=='invnu':return 3*y/(1+2*y)
        raise ValueError('unknown map')
    decreasing=code in ('complement_cross','mu')
    l=at(x.hi if decreasing else x.lo).lo
    h=at(x.lo if decreasing else x.hi).hi
    return IV(max(0,l),min(SCALE,h))

def transform_pair(v):return [mono_map(v[0],'mu'),mono_map(v[1],'nu')]
def back_pair(v):return [mono_map(v[0],'mu'),mono_map(v[1],'invnu')]
def invol(x):return mono_map(x,'complement_cross')

def support_weights(vals,base,t):
    raws=[[(t*x.log()).exp() for x in v] for v in vals]
    num=[]
    for i,a in enumerate(raws):
        total=IV.integer(1)
        for j,b in enumerate(raws):
            if i==j:continue
            dm=(a[0]-b[0]).abs().pow_int(3)
            dn=(a[1]-b[1]).abs().pow_int(3)
            distance=(dm+dn).div_int(2).root(3).intersect(0,1)
            total=total+1-distance
        num.append(base[i]*total)
    return [z/sum(num) for z in num]

def expert_mean(vals,w):
    # Products in bounded coordinates, each power enclosed by exp(w log b).
    terms=[transform_pair(x) for x in vals]
    out=[]
    for c in range(2):
        p=IV.integer(1)
        for z,weight in zip(terms,w):p=p*z[c].pow_unit(weight)
        out.append(p)
    return back_pair(out)

def partition_mean(vals,w,d):
    size=len(vals)
    b=[transform_pair(x) for x in vals]
    a=[[invol(x.pow_unit(size*ww)) for x in v] for v,ww in zip(b,w)]
    prods=[IV.integer(1),IV.integer(1)];count=0
    for permutation in itertools.permutations(range(size)):
        count+=1
        for c in range(2):
            x=IV.integer(1)
            for j,k in enumerate(permutation):x=x*a[k][c].pow_int(d[j])
            prods[c]=prods[c]*invol(x)
    # 1/m! and 1/sum(d) use exact integer roots, a different expression path
    # than the floating-power implementation of the original kernel.
    out=[invol(invol(x.root(count)).root(sum(d))) for x in prods]
    return back_pair(out)

def qin(t,inputs,model):
    if model not in ('H','A'):raise ValueError('unknown source model')
    vals=[inputs[i:i+2] for i in range(0,200,2)]
    expert=[IV.frac(F(n,100)) for n in [30,22,28,20]]
    criterion=[IV.frac(F(n,100)) for n in [20,20,15,25,20]]
    scores=[]
    for alternative in range(5):
        coll=[]
        for c in range(5):
            v=[vals[e*25+alternative*5+c] for e in range(4)]
            coll.append(expert_mean(v,support_weights(v,expert,t)))
        w=support_weights(coll,criterion,t)
        first=partition_mean([coll[j] for j in [0,2,4]],[w[j] for j in [0,2,4]],[1,2,3])
        second=partition_mean([coll[j] for j in [1,3]],[w[j] for j in [1,3]],[1,2])
        if model=='A':
            u=1-((1-first[0])*(1-second[0])).intersect(0,1).root(2)
            v=(first[1]*second[1]).root(2)
        else:
            a,b=transform_pair(first),transform_pair(second)
            u,v=back_pair([(a[c]*b[c]).root(2) for c in range(2)])
        scores.append(u-v)
    return scores
