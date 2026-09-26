"""Self-contained source-form model snapshot for independent diagnostics.

Provenance: ../../FSS_REPAIR_ROUND2_20260925/mpfr_replay/analysis/
qin2019_source_form.py (logarithmic generators, input data and parameters),
and ../../FSS_REPAIR_ROUND2_20260925/qin_provenance/confirm_direct_formula.py
(java_direct_one_partition). Only the formula definitions are retained;
the archived experiment drivers, precision mutation, and local imports are not.
The literal one-decimal source inputs are recorded as integer tenths so that
no binary64 input conversion enters this independent evaluation.

No C++ bounded-coordinate transformations or interval-generator code imported.
Both models use local partition size m and globally normalized power weights;
they differ only in the final H (Hamacher) versus A (algebraic) partition merge.
"""
import itertools
import mpmath as mp

# Expert x alternative x criterion, alternating membership/non-membership.
RAW_TENTHS = [
    ["7 2 8 2 4 5 7 1 9 2", "8 6 7 6 4 5 5 3 7 2", "6 5 5 4 5 6 8 5 8 3", "7 2 6 5 6 5 6 2 6 5", "6 4 7 5 5 6 7 4 7 3"],
    ["7 1 7 3 3 5 6 1 8 2", "8 3 7 5 2 7 9 2 6 2", "9 2 8 3 5 3 8 4 9 2", "6 3 7 5 4 6 8 3 7 6", "7 3 8 4 3 4 7 2 6 3"],
    ["8 2 6 2 3 7 8 1 7 1", "8 6 7 4 4 8 7 4 6 4", "7 2 8 3 2 6 9 2 8 2", "4 7 9 2 2 5 5 5 9 6", "6 4 7 3 3 6 7 3 8 5"],
    ["9 3 6 2 3 6 7 3 7 2", "6 6 7 7 2 6 6 2 6 2", "7 1 9 4 4 7 8 3 9 4", "7 3 5 4 4 5 9 4 6 2", "8 4 8 5 2 5 7 4 8 3"],
]
PARTS = [[0, 2, 4], [1, 3]]
DELTAS = [[1, 2, 3], [1, 2]]


def base_weights():
    # Construct inside active precision; do not reuse low-precision constants.
    return ([mp.mpf(x)/100 for x in (30,22,28,20)],
            [mp.mpf(x)/100 for x in (20,20,15,25,20)])


def dist(a,b,p=3):
    return (mp.mpf('.5')*abs(a[0]-b[0])**p+mp.mpf('.5')*abs(a[1]-b[1])**p)**(mp.mpf(1)/p)


def pweights(vals,base,p=3):
    T=[sum(1-dist(x,y,p) for j,y in enumerate(vals) if j!=i) for i,x in enumerate(vals)]
    z=[base[i]*(1+T[i]) for i in range(len(vals))]
    d=sum(z)
    return [x/d for x in z]


class H:
    def __init__(self,q,lam=3): self.q=mp.mpf(q); self.lam=mp.mpf(lam)
    def f(self,t): return mp.log((self.lam+(1-self.lam)*t**self.q)/(t**self.q))
    def g(self,t): return mp.log((self.lam+(1-self.lam)*(1-t**self.q))/(1-t**self.q))
    def fi(self,s): return (self.lam/(mp.e**s+self.lam-1))**(1/self.q)
    def gi(self,s): return ((mp.e**s-1)/(mp.e**s+self.lam-1))**(1/self.q)


def add(a,b,G): return (G.gi(G.g(a[0])+G.g(b[0])),G.fi(G.f(a[1])+G.f(b[1])))
def mul(a,b,G): return (G.fi(G.f(a[0])+G.f(b[0])),G.gi(G.g(a[1])+G.g(b[1])))
def scalar(c,a,G): c=mp.mpf(c); return (G.gi(c*G.g(a[0])),G.fi(c*G.f(a[1])))
def power(a,c,G): c=mp.mpf(c); return (G.fi(c*G.f(a[0])),G.gi(c*G.g(a[1])))


def one_partition(vals,w,deltas,G):
    m=len(vals); acc=None
    for perm in itertools.permutations(range(m)):
        prod=None
        for pos,idx in enumerate(perm):
            if deltas[pos]==0: continue
            term=power(scalar(m*w[idx],vals[idx],G),deltas[pos],G)
            prod=term if prod is None else mul(prod,term,G)
        acc=prod if acc is None else add(acc,prod,G)
    return power(scalar(mp.mpf(1)/mp.factorial(m),acc,G),mp.mpf(1)/sum(deltas),G)


def java_direct_one_partition(vals,w,deltas,gen):
    n=len(vals); q=gen.q; lam=gen.lam
    products=[mp.mpf(1),mp.mpf(1)]
    for perm in itertools.permutations(range(n)):
        p0=p1=mp.mpf(1)
        for j,idx in enumerate(perm):
            mu,nu=vals[idx]
            ex=n*w[idx]
            am=(lam+(1-lam)*(1-mu**q))**ex
            bm=(1-mu**q)**ex
            an=(lam+(1-lam)*nu**q)**ex
            bn=nu**(q*ex)
            p0*=((am+(lam**2-1)*bm)/(am-bm))**deltas[j]
            p1*=((an+(lam**2-1)*bn)/(an-bn))**deltas[j]
        products[0]*=(p0+lam**2-1)/(p0-1)
        products[1]*=(p1+lam**2-1)/(p1-1)
    fact=1/mp.factorial(n); s=1/mp.mpf(sum(deltas))
    p0,p1=(v**fact for v in products)
    mu=(lam*(p0-1)**s/((p0+lam**2-1)**s+(lam-1)*(p0-1)**s))**(1/q)
    nu=(((p1+lam**2-1)**s-(p1-1)**s)/((p1+lam**2-1)**s+(lam-1)*(p1-1)**s))**(1/q)
    return mu,nu


def scores(canonical,q,partition,fixed_base_weights=False):
    """Both source-qualified merge scores with exact active-precision q."""
    q=mp.mpf(q)
    generator=H(q,3)
    ew,cw=base_weights()
    both=[[],[]]
    for alternative in range(5):
        collective=[]
        for criterion in range(5):
            values=[tuple(canonical[2*(expert*25+alternative*5+criterion)+component]**(mp.mpf(1)/q)
                          for component in range(2)) for expert in range(4)]
            weights=ew if fixed_base_weights else pweights(values,ew,3)
            collective.append(partition(values,weights,[1,0,0,0],generator))
        weights=cw if fixed_base_weights else pweights(collective,cw,3)
        parts=[partition([collective[i] for i in indices],[weights[i] for i in indices],deltas,generator)
               for indices,deltas in zip(PARTS,DELTAS)]
        hamacher=scalar(mp.mpf(1)/2,add(parts[0],parts[1],generator),generator)
        both[0].append(hamacher[0]**q-hamacher[1]**q)
        algebraic_u=1-mp.sqrt((1-parts[0][0]**q)*(1-parts[1][0]**q))
        algebraic_v=mp.sqrt(parts[0][1]**q*parts[1][1]**q)
        both[1].append(algebraic_u-algebraic_v)
    return both
