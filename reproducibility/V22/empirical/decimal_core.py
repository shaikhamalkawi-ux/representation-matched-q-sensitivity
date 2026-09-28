"""Unchanged R15 independent Decimal kernels; no primary evaluator import."""
from decimal import Decimal as D, localcontext
from fractions import Fraction as F

def decimal_fraction(value):
    r = F(value)
    return D(r.numerator) / D(r.denominator)


def decimal_pi():
    def atan(x):
        term, total, n = x, x, 0
        while True:
            n += 1
            term *= -x * x
            add = term / D(2*n+1)
            total += add
            if abs(add) < D("1e-89"):
                return total
    return 16*atan(D(1)/5) - 4*atan(D(1)/239)


def cosine(x):
    term, total, n = D(1), D(1), 0
    while True:
        n += 1
        term *= -x*x / D((2*n-1)*(2*n))
        total += term
        if abs(term) < D("1e-89"):
            return total


def power(x, w):
    if x == 0:
        assert w > 0
        return D(0)
    return (w*x.ln()).exp()


def product_power(values, weights):
    if any(x == 0 for x in values):
        return D(0)
    return sum((w*x.ln() for x, w in zip(values, weights)), D(0)).exp()


def reduced(u, v, w, pi):
    wa = 1-product_power([1-x for x in u], w)-product_power(v, w)
    wg = product_power(u, w)+product_power([1-x for x in v], w)-1
    plus = sum((a*cosine(pi*((1+y)**2-x*x)/8) for x,y,a in zip(u,v,w)), D(0))
    minus = sum((a*cosine(pi*((1+x)**2-y*y)/8) for x,y,a in zip(u,v,w)), D(0))
    assert plus+minus >= 1-D("1e-85")
    return wa, wg, plus/(plus+minus)


def literal_transport(u, v, w, q, pi):
    x, y = [power(a, 1/q) for a in u], [power(a, 1/q) for a in v]
    xp, yp = [power(a, q) for a in x], [power(a, q) for a in y]
    awm = power(1-product_power([1-a for a in xp],w),1/q)
    awn = product_power(y,w)
    gwm = product_power(x,w)
    gwn = power(1-product_power([1-a for a in yp],w),1/q)
    similarities = []
    for a,b in [(D(1),D(0)), (D(0),D(1))]:
        similarities.append(sum((weight*cosine(pi/4*(abs(s-a)+abs(t-b))
                             *(1-abs(1-s-t)/2))
                             for s,t,weight in zip(xp,yp,w)),D(0)))
    return power(awm,q)-power(awn,q), power(gwm,q)-power(gwn,q), \
        similarities[0]/sum(similarities)


def synthetic_gates():
    pi = decimal_pi()
    # Expected strings here are independent landmarks, not primary output values.
    assert abs(pi-D("3.14159265358979323846264338327950288419716939937510")) < D("1e-49")
    assert abs(cosine(pi/2)) < D("1e-85")
    assert abs(cosine(pi/3)-D("0.5")) < D("1e-85")
    weights = [D("0.3"),D("0.7")]
    examples = [([1,1],[0,0]),([0,0],[1,1]),([0,0],[0,0]),
                ([1,0],[0,1]),([D("0.4"),D("0.7")],[D("0.2"),D("0.1")]),
                ([0,D("0.7")],[D("0.2"),0])]
    max_error = D(0)
    comparisons = 0
    for us,vs in examples:
        u,v = list(map(D,us)),list(map(D,vs))
        expected = reduced(u,v,weights,pi)
        assert -1-D("1e-80") <= expected[0] <= 1+D("1e-80")
        assert -1-D("1e-80") <= expected[1] <= 1+D("1e-80")
        assert -D("1e-80") <= expected[2] <= 1+D("1e-80")
        for q in (D(1),D("1.5"),D(4),D(16)):
            actual = literal_transport(u,v,weights,q,pi)
            for a,b in zip(expected,actual):
                max_error=max(max_error,abs(a-b)); comparisons+=1
                assert abs(a-b)<D("1e-75")
        reverse = reduced(u[::-1],v[::-1],weights[::-1],pi)
        assert all(abs(a-b)<D("1e-80") for a,b in zip(expected,reverse))
    for u,v in [(D("0.6"),D("0.2")),(D(0),D(0)),(D(1),D(0)),(D(0),D(1))]:
        z = reduced([u],[v],[D(1)],pi)
        assert abs(z[0]-(u-v))<D("1e-80")
        assert abs(z[1]-(u-v))<D("1e-80")
        repeated=reduced([u,u],[v,v],weights,pi)
        assert all(abs(a-b)<D("1e-80") for a,b in zip(z,repeated))
    worst,perfect,empty,mixed = [reduced(list(map(D,a)),list(map(D,b)),weights,pi)
                                for a,b in [([0,0],[1,1]),([1,1],[0,0]),
                                            ([0,0],[0,0]),([1,0],[0,1])]]
    assert max(abs(a-b) for a,b in zip(worst,(-1,-1,0)))<D("1e-80")
    assert max(abs(a-b) for a,b in zip(perfect,(1,1,1)))<D("1e-80")
    assert abs(empty[2]-D("0.5"))<D("1e-80")
    assert abs(mixed[0]-1)<D("1e-80") and abs(mixed[1]+1)<D("1e-80")
    return {"status":"PASS_SYNTHETIC_ONLY", "literal_reduced_scalar_comparisons":comparisons,
            "max_absolute_error":str(max_error), "pi":str(pi),
            "gates":["Machin pi and cosine landmarks", "literal raw transport/reduction",
                     "criterion permutation", "one criterion", "identical criteria",
                     "perfect/worst/empty/mixed boundaries", "native score ranges"],
            "not_real_data_trials":True}


def caches_for(rows,pi):
    u=[[decimal_fraction(x) for x in r["mu"]] for r in rows]
    v=[[decimal_fraction(x) for x in r["nu"]] for r in rows]
    logu=[[None if x==0 else x.ln() for x in r] for r in u]
    logv=[[None if x==0 else x.ln() for x in r] for r in v]
    logs=[]; cells=[]
    for us,vs in zip(u,v):
        logs.append([[None if x==0 else x.ln() for x in factors]
                     for factors in ([1-x for x in us],vs,us,[1-x for x in vs])])
        cells.append([(cosine(pi*((1+y)**2-x*x)/8),cosine(pi*((1+x)**2-y*y)/8))
                      for x,y in zip(us,vs)])
    return logu,logv,logs,cells


def entropy_weights(logu,logv,q):
    count,dim=len(logu),len(logu[0])
    values=[]
    for j in range(dim):
        total=D(0)
        for ul,vl in zip(logu,logv):
            for s in (ul[j],vl[j]):
                if s is not None:
                    a=s/q; total+=a.exp()*a
        values.append(1+total/D(count))
    assert min(values)>0
    return [a/sum(values) for a in values]


def scores(logs,cells,w):
    result=[[],[],[]]
    for products,pairs in zip(logs,cells):
        pp=[D(0) if any(x is None for x in factors) else
            sum((weight*x for weight,x in zip(w,factors)),D(0)).exp() for factors in products]
        result[0].append(1-pp[0]-pp[1]); result[1].append(pp[2]+pp[3]-1)
        plus=sum((weight*pair[0] for weight,pair in zip(w,pairs)),D(0))
        minus=sum((weight*pair[1] for weight,pair in zip(w,pairs)),D(0))
        assert plus+minus>=1-D("1e-85")
        result[2].append(plus/(plus+minus))
    return result


def fmt(x):
    return format(x,".80g")
