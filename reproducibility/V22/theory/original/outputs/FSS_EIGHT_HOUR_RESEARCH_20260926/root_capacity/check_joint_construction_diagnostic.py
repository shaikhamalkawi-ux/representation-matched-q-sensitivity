"""Post-R4 full-table numerical diagnostic, NOT a certificate or a proof."""
from __future__ import annotations
import json
from pathlib import Path
from decimal import Decimal, getcontext

getcontext().prec = 90


class DecimalMath:
    """Small stdlib adapter; ordinary rounding, deliberately not interval proof."""
    mpf = Decimal
    inf = Decimal('Infinity')
    log = staticmethod(lambda x: Decimal(x).ln())
    exp = staticmethod(lambda x: Decimal(x).exp())
    sign = staticmethod(lambda x: (x > 0)-(x < 0))
    nstr = staticmethod(lambda x, n: format(x, f'.{n}g'))


mp = DecimalMath()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def check(m, k, n, r):
    t = (m-k)//n
    neutral = m-k-n*t
    require(k >= 2 and n >= 2 and t >= 1 and neutral >= 0, 'row budget')
    v = mp.mpf(1)/8
    delta = 11*mp.log(2)*t/(480*m)
    slopes = [mp.mpf(2*i-k-1)/(4*(k-1)) for i in range(1,k+1)]
    rows = []
    for slope in slopes:
        ell0 = -mp.log(2)-slope/2+delta*slope*slope
        ell1 = -mp.log(2)+slope/2+delta*slope*slope
        rows.append(((1-mp.exp(ell1))**(mp.mpf(1)/r),
                     (1-mp.exp(ell0))**(mp.mpf(1)/r)))
    for ell in range(n):
        a = mp.log(16)*16**ell
        rows.extend([(mp.exp(-a),mp.exp(-2*a))]*t)
    rows.extend([(mp.mpf(1)/16,mp.mpf(1)/16)]*neutral)
    require(len(rows) == m, 'm rows')
    require(all(0<x<1 and x**r+v**r<1 for row in rows for x in row), 'domain')
    logs = [[mp.log(x) for x in row] for row in rows]

    def weight(q):
        p = mp.mpf(r)/q
        vl = p*mp.log(v)
        e = [1+sum(mp.exp(p*row[j])*p*row[j]+mp.exp(vl)*vl for row in logs)/m
             for j in range(2)]
        return e[0]/sum(e)

    def scores(q):
        w = weight(q)
        return [1-mp.exp(w*mp.log(1-row[0]**r)+(1-w)*mp.log(1-row[1]**r))-v**r
                for row in rows]

    macro = [mp.mpf(r)*factor*16**ell for ell in range(n) for factor in (2,8)]
    endpoint_weights = [weight(q) for q in macro]
    require(all((w <= mp.mpf('.5')-delta if i%2 == 0 else w >= mp.mpf('.5')+delta)
                for i,w in enumerate(endpoint_weights)), 'entropy excursion')
    nodes = []
    labels = []
    min_gap = mp.inf
    for interval in range(len(macro)-1):
        increasing = endpoint_weights[interval] < endpoint_weights[interval+1]
        targets = sorted([(mp.mpf('.5')-2*delta*slope,i) for i,slope in enumerate(slopes)],
                         reverse=not increasing)
        left = macro[interval]
        right_end = macro[interval+1]
        for target, expected in targets:
            lo, hi = left, right_end
            sign_lo = mp.sign(weight(lo)-target)
            require(sign_lo*mp.sign(weight(hi)-target)<0, 'target bracket')
            for _ in range(180):
                mid = (lo+hi)/2
                if mp.sign(weight(mid)-target) == sign_lo:
                    lo = mid
                else:
                    hi = mid
            q = (lo+hi)/2
            ss = scores(q)
            gap = ss[expected]-max(x for j,x in enumerate(ss) if j != expected)
            require(gap > 0, 'unique global candidate winner')
            require(not nodes or q > nodes[-1], 'strict node chronology')
            nodes.append(q)
            labels.append(expected)
            min_gap = min(min_gap,gap)
            left = q
    analytic_gap = 77*mp.log(2)*t/(30720*m*(k-1)**2)
    require(min_gap >= analytic_gap, 'analytic conservative margin diagnostic')
    transitions = sum(a != b for a,b in zip(labels,labels[1:]))
    require(transitions == (2*n-1)*(k-1), 'certified-pattern transition target')
    return dict(m=m,k=k,n=n,r=r,multiplicity=t,neutral_rows=neutral,
                sample_count=len(nodes),transitions=transitions,
                minimum_sample_gap=mp.nstr(min_gap,35),
                analytic_lower_gap=mp.nstr(analytic_gap,35),
                ratio=mp.nstr(min_gap/analytic_gap,20),passed=True)


def main():
    cases=[check(m,k,n,r) for m,k,n in [(4,2,2),(7,3,3),(8,3,2),(12,5,3),(20,2,2),(11,2,3)]
           for r in (1,4,16)]
    result=dict(status='PASS_NUMERICAL_DIAGNOSTIC_ONLY',precision_decimal_digits=90,
                not_a_continuum_or_outward_certificate=True,
                implementation='Stdlib Decimal90; full shared table and all row scores; no injected weight path',cases=cases)
    path=Path(__file__).with_name('POST_R4_JOINT_CONSTRUCTION_DIAGNOSTIC.json')
    if path.exists():
        require(json.loads(path.read_text(encoding='utf-8')) == result,
                'Existing diagnostic differs; preserve it and investigate')
        mode='READ_COMPARE_EXISTING_NO_OVERWRITE'
    else:
        path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
        mode='CREATE_NEW'
    print(json.dumps({'status':result['status'],'cases':len(cases),'output':str(path),'mode':mode}))


if __name__ == '__main__':
    main()
