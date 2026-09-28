"""V17 containing enclosures of the unchanged Seikh model, integer arithmetic.

Only the interval evaluation is strengthened: h(y)=y*exp(y) is enclosed by
its endpoint values and its unique critical minimum at -1. All input units,
entropy normalization and score formulas remain the V16 definitions.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import sys
import time

HERE = Path(__file__).resolve().parent
EXPECTED_VENDOR_HASHES = {
    'model_inputs.json': '71351416bb73b22f3781bad962a023b6ea7362e8020429aab97b47d3f36d1797',
    'exact_models.py': '8d5941c32cd7cf7140c3273ac6805a370f73bcea627d321cd2ba3c2bfc62f427',
    'dyadic_interval.py': '38cae911d1feb49562dcfd3693bb56a8e1bb8219c8afa51f89d44aff5998852e',
}
for _name, _expected in EXPECTED_VENDOR_HASHES.items():
    if hashlib.sha256((HERE / 'vendor' / _name).read_bytes()).hexdigest() != _expected:
        raise ValueError('Immutable model/arithmetic source hash mismatch: ' + _name)
sys.path.insert(0, str(HERE / "vendor"))
from dyadic_interval import IV, SCALE, PREC
from exact_models import domain, tighter_weights, simplex_linear_bounds

INPUT_SHA = EXPECTED_VENDOR_HASHES['model_inputs.json']


def entropy_term(y):
    if y.hi > 0:
        raise ValueError("entropy argument must be nonpositive")
    left, right = IV(y.lo, y.lo), IV(y.hi, y.hi)
    a, b = left * left.exp(), right * right.exp()
    lows = [a.lo, b.lo]
    if y.lo <= -SCALE <= y.hi:
        lows.append((-IV.integer(-1).exp()).lo)
    return IV(min(lows), max(a.hi, b.hi))


def scores(t, inputs):
    if not (0 <= t.lo <= t.hi <= SCALE // 4):
        raise ValueError("t outside [0,1/4]")
    logs = [x.log() for x in inputs]
    factors = []
    for j in range(5):
        total = sum((entropy_term(t * logs[(i*5+j)*2+c])
                     for i in range(4) for c in range(2)), IV.integer(0))
        e = 1 + total.div_int(4)
        factors.append(IV(max(e.lo, SCALE//4), min(e.hi, SCALE)))
    weights = tighter_weights(factors)
    result = []
    for i in range(4):
        pos = [(1-inputs[(i*5+j)*2]).log() for j in range(5)]
        neg = [logs[(i*5+j)*2+1] for j in range(5)]
        p = simplex_linear_bounds(weights, pos).exp()
        n = simplex_linear_bounds(weights, neg).exp()
        result.append(1-p-n)
    return result, weights


def certify(radius, q=None, max_nodes=2047, max_depth=18):
    radius = F(radius)
    if q is not None and F(q) < 4:
        raise ValueError("fixed q must be >=4")
    root = (F(0), F(1,4)) if q is None else (F(1,F(q)), F(1,F(q)))
    inputs = domain('seikh', 'source', radius)
    pending = [(root[0], root[1], 'r', 0)]
    nodes, leaves, unresolved = [], [], []
    start = time.perf_counter()
    while pending:
        a,b,ident,depth = pending.pop()
        sc,w = scores(IV.frac(a,b), inputs)
        margins = [F(sc[1].lo-sc[j].hi,SCALE) for j in (0,2,3)]
        good = min(margins)>0
        rec = {'id':ident, 't':[str(a),str(b)],
               'scores':[x.endpoints() for x in sc],
               'weights':[x.endpoints() for x in w],
               'margins':[str(x) for x in margins], 'positive':good}
        nodes.append(rec)
        if good:
            leaves.append(rec)
        elif a==b or depth>=max_depth or len(nodes)+len(pending)+2>max_nodes:
            unresolved.append(rec)
        else:
            midpoint=(a+b)/2
            pending.extend([(midpoint,b,ident+'1',depth+1),(a,midpoint,ident+'0',depth+1)])
    return {'status':'PASS_COMPUTATIONAL' if not unresolved else 'UNRESOLVED',
            'version':17, 'family':'Seikh-Mandal',
            'model':'formula entropy and product-generator qROFAWA, unchanged',
            'enclosure':'scalar_yexp_endpoint_critical+weight_simplex_LP',
            'arithmetic':'Python integer dyadic intervals with rational remainder bounds',
            'bits':PREC, 'input_model_sha256':INPUT_SHA,
            'normalized_baseline_raw_halfwidth':str(radius), 'baseline_rung':4,
            'root_t':[str(x) for x in root], 'q_exact':None if q is None else str(F(q)),
            'scope':'every finite real q>=4; t=0 limiting endpoint only' if q is None else 'one fixed finite rung',
            'winner':'Y2', 'node_count':len(nodes), 'leaf_count':len(leaves),
            'unresolved_count':len(unresolved),
            'margin_floor':str(min((F(v) for n in leaves for v in n['margins']),default=F(0))),
            'max_nodes':max_nodes, 'max_depth':max_depth,
            'elapsed_seconds':time.perf_counter()-start,
            'nodes':nodes, 'leaves':leaves, 'unresolved':unresolved}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--radius', required=True)
    p.add_argument('--q')
    p.add_argument('--max-nodes',type=int,default=2047)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():
        raise FileExistsError('Refusing to overwrite a research result')
    result=certify(F(args.radius),args.q,args.max_nodes)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('nodes','leaves','unresolved')},indent=2))

if __name__=='__main__':
    main()
