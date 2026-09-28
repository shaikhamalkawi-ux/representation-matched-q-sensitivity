"""Independent literal Decimal100 replay of a finite-node raw source-box claim.

No imports from the construction, discovery, or MPFR model evaluator. The
separate Decimal arithmetic kernel is pinned; all formulas are restated here.
The JSON input supplies exact source rationals, rational q nodes, expected
winner labels and one all-coordinate baseline-raw L-infinity halfwidth.
Every coordinate can vary independently; no common-nu/symmetry constraint.
This is a finite-sequence preservation check, not a continuum winner/radius
optimum claim, an empirical uncertainty calibration, or a new general theorem.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import copy
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
KERNEL = HERE.parent / 'entropy_family/pinar_independent/replay_decimal.py'
KERNEL_SHA = '82fc5a524072a0da7dac34aafa53e3b96db39df727d42c3bc8158b3adfd36235'
if hashlib.sha256(KERNEL.read_bytes()).hexdigest() != KERNEL_SHA:
    raise ValueError('independent Decimal kernel changed')
SPEC = importlib.util.spec_from_file_location('compact_independent_decimal', KERNEL)
kernel = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(kernel)
Box = kernel.Box
require = kernel.require


def rational(value):
    require(type(value) in (str, int), 'rational must be a string or integer, never float/bool')
    return F(value)


def parse(data):
    require(type(data) is dict, 'source spec must be an object')
    r = rational(data['r'])
    require(r >= 1, 'baseline below1')
    h = rational(data['halfwidth'])
    require(h > 0, 'strictly positive source halfwidth required')
    labels = data['labels']
    require(type(labels) is list and len(labels) >= 2, 'need at least2 labels')
    require(all(type(x) is str and x for x in labels), 'nonempty string labels required')
    require(len(set(labels)) == len(labels), 'duplicate labels')
    source = data['source']
    require(type(source) is list and len(source) == len(labels), 'row count mismatch')
    source_boxes = []
    exact = []
    for row in source:
        require(type(row) is list and len(row) == 2, 'exactly2 criteria required')
        boxes = []
        exact_row = []
        for pair in row:
            require(type(pair) is list and len(pair) == 2, 'orthopair length')
            mu, nu = map(rational, pair)
            require(0 < mu-h < mu+h < 1 and 0 < nu-h < nu+h < 1, 'raw source box outside interior')
            mb, nb = Box(mu-h, mu+h), Box(nu-h, nu+h)
            require((mb.power(r)+nb.power(r)).hi < 1, 'whole raw source box not strictly r-admissible')
            boxes.append((mb, nb))
            exact_row.append((mu, nu))
        source_boxes.append(boxes)
        exact.append(exact_row)
    nodes = data['nodes']
    require(type(nodes) is list and len(nodes) >= 2, 'need ordered nodes')
    checked = []
    previous = None
    for node in nodes:
        require(type(node) is dict, 'node must be object')
        q = rational(node['q'])
        winner = node['winner']
        require(q >= r and (previous is None or q > previous), 'q below baseline or not strictly ordered')
        require(winner in labels, 'unknown winner')
        previous = q
        checked.append((q, labels.index(winner)))
    require(any(a[1] != b[1] for a, b in zip(checked, checked[1:])), 'no claimed winner changes')
    return r, h, labels, exact, source_boxes, checked


def product(values):
    out = Box(1)
    for value in values:
        out *= value
    return out


def literal_pipeline(matrix, r, q):
    m = len(matrix)
    transported = [[(mu.power(r/q), nu.power(r/q)) for mu, nu in row] for row in matrix]
    factors = [1 + sum(pair[j][0]*pair[j][0].log()+pair[j][1]*pair[j][1].log()
                       for pair in transported)/m for j in range(2)]
    require(all(e.lo > 0 for e in factors), 'nonpositive entropy factor')
    total = sum(factors)
    weights = [e/total for e in factors]
    require(all(w.lo > 0 and w.hi < 1 for w in weights), 'weight interval not interior')
    scores = []
    for row in transported:
        membership = (1-product((1-mu.power(q)).power(w) for (mu, nu), w in zip(row, weights))).power(1/q)
        nonmembership = product(nu.power(w) for (mu, nu), w in zip(row, weights))
        scores.append(membership.power(q)-nonmembership.power(q))
    return factors, weights, scores


def verify(data):
    r, h, labels, exact, matrix, nodes = parse(data)
    records = []
    center = [[(Box(mu), Box(nu)) for mu, nu in row] for row in exact]
    for q, winner in nodes:
        factors, weights, scores = literal_pipeline(matrix, r, q)
        _, _, center_scores = literal_pipeline(center, r, q)
        gaps = [scores[winner]-other for i, other in enumerate(scores) if i != winner]
        require(all(g.lo > 0 for g in gaps), 'claimed whole-box winner not certified')
        for interval, c in zip(scores, center_scores):
            require(interval.lo <= c.lo <= c.hi <= interval.hi, 'center not contained in box enclosure')
        records.append({'q':str(q), 'winner':labels[winner],
                        'minimum_box_gap_lower':str(min(g.lo for g in gaps)),
                        'entropy_factors':[e.pack() for e in factors],
                        'weights':[w.pack() for w in weights],
                        'scores':[s.pack() for s in scores],
                        'center_scores':[s.pack() for s in center_scores]})
    changes = sum(a[1] != b[1] for a, b in zip(nodes, nodes[1:]))
    return {'status':'PASS_INDEPENDENT_DECIMAL100_RAW_SOURCE_BOX',
            'r':str(r), 'm':len(labels), 'independently_variable_raw_coordinates':4*len(labels),
            'halfwidth':str(h), 'q_window':[str(nodes[0][0]),str(nodes[-1][0])],
            'winner_sequence':[labels[i] for q,i in nodes], 'at_least_changes':changes,
            'minimum_box_gap_lower':min(records,key=lambda x:F(x['minimum_box_gap_lower']))['minimum_box_gap_lower'],
            'precision_decimal_digits':100, 'records':records,
            'scope':'Every fixed raw source matrix in the all-coordinate L-infinity box preserves these finite-node winners. By continuity this forces at least the reported number of winner changes. No common nonmembership/symmetry restriction on box members; no practical calibration, optimal radius, exact event count, or common-nu global theorem generalization.'}


def malformed_tests(data):
    mutations = []
    def add(label, change):
        altered = copy.deepcopy(data)
        change(altered)
        mutations.append((label, altered))
    add('float baseline', lambda x:x.update(r=4.0))
    add('boolean baseline', lambda x:x.update(r=True))
    add('zero radius', lambda x:x.update(halfwidth='0'))
    add('negative radius', lambda x:x.update(halfwidth='-1/100'))
    add('outside raw domain', lambda x:x['source'][0][0].__setitem__(0,'2'))
    add('nonadmissible pair', lambda x:x['source'][0].__setitem__(0,['999/1000','999/1000']))
    add('third criterion', lambda x:x['source'][0].append(['1/2','1/8']))
    add('duplicate label', lambda x:x['labels'].__setitem__(0,x['labels'][1]))
    add('unknown winner', lambda x:x['nodes'][0].update(winner='ABSENT'))
    add('unordered nodes', lambda x:x['nodes'].__setitem__(1,copy.deepcopy(x['nodes'][0])))
    add('wrong winner', lambda x:x['nodes'][0].update(winner=next(t for t in x['labels'] if t!=x['nodes'][0]['winner'])))
    rejected = []
    for label, altered in mutations:
        try:
            verify(altered)
        except (ValueError, KeyError, TypeError, ZeroDivisionError):
            rejected.append(label)
        else:
            raise ValueError('malformed input accepted: '+label)
    return rejected


def main(source, output):
    raw = source.read_bytes()
    original = json.loads(raw)
    require(original.get('schema') == 'compact-multiwinner-raw-box-v1', 'unknown frozen source schema')
    data = {'r':original['baseline_r'], 'halfwidth':original['uniform_raw_halfwidth'],
            'source':original['source'], 'labels':original['labels'], 'nodes':original['q_nodes']}
    result = verify(data)
    result['malformed_input_rejections'] = malformed_tests(data)
    result['primitive_selftest_checks'] = kernel.selftests()
    result['source_spec_sha256'] = hashlib.sha256(raw).hexdigest()
    result['decimal_kernel_sha256'] = KERNEL_SHA
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    with output.open('x',encoding='utf8') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({key:result[key] for key in ['status','m','halfwidth','at_least_changes','minimum_box_gap_lower','source_spec_sha256']},indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    main(args.source,args.output)
