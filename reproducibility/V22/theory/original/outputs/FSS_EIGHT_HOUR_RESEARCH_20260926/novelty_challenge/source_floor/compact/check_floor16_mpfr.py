"""Research-only fixed-floor m16 witness; finite-node directed MPFR256 checks.

Discovery uses ordinary binary64 only to choose integer q nodes. All source
rounding and final center/full64-coordinate box checks use directed MPFR256.
Both literal transported and canonical two-product score formulas are checked.
No common-nu or mirrored-coordinate constraint is imposed on box members.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import copy
import hashlib
import importlib.util
import json
import math

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
KERNEL = ROOT / 'outputs/FSS_STRENGTHENING_20260926/radius_refinement/mpfr_audit/verify_refined_mpfr.py'
KERNEL_SHA = '867b871950cf41960ee1274411b2581cda88377ecdc512617991c1471487592a'
if hashlib.sha256(KERNEL.read_bytes()).hexdigest() != KERNEL_SHA:
    raise ValueError('pinned directed kernel changed')
spec = importlib.util.spec_from_file_location('floor16_mpfr', KERNEL)
kernel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(kernel)
I, require = kernel.I, kernel.require
LABELS = [f'C{i}' for i in range(1, 9)] + [f'D{i}' for i in range(1, 9)]
PATTERN = ['C8', 'C7', 'C6', 'C5']


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def power(x, exponent):
    return (x.log() * (exponent if isinstance(exponent, I) else I(exponent))).exp()


def product(values):
    ans = I(1)
    for value in values:
        ans *= value
    return ans


def pack(x):
    return {'lower': str(x.lo), 'upper': str(x.hi), 'exact_dyadic_endpoints': x.pair()}


def save(path, value):
    with path.open('x', encoding='utf8') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')


def freeze_source():
    delta = (I(F(1, 4)) * I(F(-1, 4)).exp()
             - I(F(1, 6)) * I(F(-1, 6)).exp()) / 8
    rows = []
    rounding = []
    for i in range(1, 9):
        c = F(2*i-9, 28)
        row = []
        for w in [F(1), F(0)]:
            ell = -I(2).log() + c*(w-F(1, 2)) + delta*c*c
            mu = power(1-ell.exp(), F(1, 4))
            lo = round(kernel.rat(mu.lo)*10**6)
            hi = round(kernel.rat(mu.hi)*10**6)
            require(lo == hi, 'six-place rounding unresolved')
            value = F(lo, 10**6)
            row.append([str(value), '1/8'])
            rounding.append({'candidate': f'C{i}', 'criterion': int(w == 0)+1,
                             'unrounded_raw_membership': pack(mu),
                             'frozen_six_decimal_membership': f'{float(value):.6f}'})
        rows.append(row)
    require(all(rows[i] == rows[7-i][::-1] for i in range(8)), 'mirror symmetry lost')
    rows += [[['1/8', '1/8'], ['1/4', '1/8']] for _ in range(8)]
    return rows, pack(delta), rounding


def diagnostic(source, q):
    p = 4.0/q
    raw = [[(float(F(mu))**p, float(F(nu))**p) for mu, nu in row] for row in source]
    entropy = [1+sum(row[j][0]*math.log(row[j][0])+row[j][1]*math.log(row[j][1])
                     for row in raw)/16 for j in range(2)]
    weights = [x/sum(entropy) for x in entropy]
    scores = [1-math.prod((1-float(F(mu))**4)**w for (mu, nu), w in zip(row, weights))
                -math.prod((float(F(nu))**4)**w for (mu, nu), w in zip(row, weights))
              for row in source]
    order = sorted(range(16), key=lambda i: scores[i], reverse=True)
    return LABELS[order[0]], scores[order[0]]-scores[order[1]], weights[0]


def discover(output):
    source, delta, rounding = freeze_source()
    best = {}
    ranges = {}
    for q in range(34, 5001):
        winner, margin, w = diagnostic(source, q)
        ranges.setdefault(winner, [q, q, 0])
        ranges[winner][1] = q
        ranges[winner][2] += 1
        if winner in PATTERN and (winner not in best or margin > best[winner]['binary64_margin']):
            best[winner] = {'q': str(q), 'winner': winner, 'binary64_margin': margin,
                            'binary64_weight1': w}
    require(set(best) == set(PATTERN), 'not all four candidates found at integer nodes')
    nodes = [{'q': best[label]['q'], 'winner': label} for label in PATTERN]
    require(all(int(a['q']) < int(b['q']) for a, b in zip(nodes, nodes[1:])), 'wrong chronology')
    record = {'schema': 'floor16-integer-node-source-v1', 'baseline_r': '4',
              'arithmetic_bits': 256, 'labels': LABELS, 'source': source,
              'source_sha256': hashlib.sha256(canonical_bytes(source)).hexdigest(),
              'candidate_decimal_places': 6, 'candidate_rounding_rule': 'nearest, ties to even',
              'exact_transcendental_design_delta_enclosure': delta,
              'source_rounding_evidence': rounding, 'q_nodes': nodes,
              'requested_raw_halfwidth': '1/2000000',
              'binary64_discovery_only': {'integer_q_range': [34, 5000], 'best_nodes': best,
                                        'winning_ranges': ranges},
              'interval_kernel_sha256': KERNEL_SHA,
              'scope': 'Synthetic fixed m16/r4 exact six-decimal candidate source; integer-node sequence only. Directed verification is separate. No published data, whole-q-path certificate, exact switch count, calibrated uncertainty or maximal source radius.'}
    save(output, record)
    print(json.dumps({'nodes': best, 'source_sha256': record['source_sha256']}, indent=2))


def source_box(source, h):
    require(len(source) == 16 and all(len(row) == 2 for row in source), '16 by 2 source required')
    require(h >= 0, 'negative halfwidth')
    result = []
    for row in source:
        target = []
        for pair in row:
            require(len(pair) == 2, 'orthopair shape')
            mu, nu = map(F, pair)
            require(mu >= F(1, 8) and nu >= F(1, 8), 'center source floor')
            require(0 < mu-h <= mu+h < 1 and 0 < nu-h <= nu+h < 1, 'whole box interior')
            require((mu+h)**4 + (nu+h)**4 < 1, 'whole box strict r4 admissibility')
            target.append((I(mu-h, mu+h), I(nu-h, nu+h)))
        result.append(target)
    return result


def evaluate(box, q, mode):
    q = F(q)
    require(q >= 4, 'q below baseline')
    raw = [[(power(mu, 4/q), power(nu, 4/q)) for mu, nu in row] for row in box]
    entropy = [1+sum(row[j][0]*row[j][0].log()+row[j][1]*row[j][1].log()
                     for row in raw)/16 for j in range(2)]
    require(all(x.lo > 0 for x in entropy), 'entropy factors not positive')
    weights = [x/sum(entropy) for x in entropy]
    require(all(0 < w.lo <= w.hi < 1 for w in weights), 'weights not interior')
    scores = []
    if mode == 'literal':
        for row in raw:
            agg_mu = power(1-product(power(1-power(mu, q), w)
                                    for (mu, nu), w in zip(row, weights)), 1/q)
            agg_nu = product(power(nu, w) for (mu, nu), w in zip(row, weights))
            scores.append(power(agg_mu, q)-power(agg_nu, q))
    elif mode == 'canonical':
        for row in box:
            scores.append(1-product(power(1-power(mu, 4), w) for (mu, nu), w in zip(row, weights))
                          -product(power(power(nu, 4), w) for (mu, nu), w in zip(row, weights)))
    else:
        raise ValueError('unknown model mode')
    return scores, entropy, weights


def checked(record, h, mode):
    box = source_box(record['source'], h)
    rows = []
    for node in record['q_nodes']:
        scores, entropy, weights = evaluate(box, F(node['q']), mode)
        winner = LABELS.index(node['winner'])
        gaps = [scores[winner]-s for i, s in enumerate(scores) if i != winner]
        rows.append({'q': node['q'], 'winner': node['winner'],
                     'winner_verified': all(g.lo > 0 for g in gaps),
                     'minimum_margin_lower': str(min(g.lo for g in gaps)),
                     'all_scores': [pack(s) for s in scores],
                     'all_rival_gaps': [pack(g) for g in gaps],
                     'entropy_factors': [pack(e) for e in entropy],
                     'weights': [pack(w) for w in weights]})
    return {'model_mode': mode, 'halfwidth': str(h),
            'all_winners_verified': all(row['winner_verified'] for row in rows), 'records': rows}


def validate(record):
    require(record['schema'] == 'floor16-integer-node-source-v1', 'schema')
    require(record['baseline_r'] == '4' and record['arithmetic_bits'] == 256, 'baseline/precision')
    require(record['labels'] == LABELS, 'labels')
    require(record['candidate_decimal_places'] == 6, 'candidate precision')
    require(record['interval_kernel_sha256'] == KERNEL_SHA, 'kernel pin')
    require(record['source_sha256'] == hashlib.sha256(canonical_bytes(record['source'])).hexdigest(), 'source hash')
    for row in record['source']:
        for pair in row:
            require(all((F(x)*10**6).denominator == 1 for x in pair), 'source not on six-place grid')
    nodes = record['q_nodes']
    require([n['winner'] for n in nodes] == PATTERN, 'winner sequence')
    qs = [F(n['q']) for n in nodes]
    require(all(q.denominator == 1 and q >= 4 for q in qs), 'integer q required')
    require(all(a < b for a, b in zip(qs, qs[1:])), 'chronology')
    require(F(record['requested_raw_halfwidth']) == F(1, 2000000), 'requested six-place halfunit')


def verify(source, output):
    original = source.read_bytes()
    own_before = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    record = json.loads(original)
    validate(record)
    attempts = []
    for h in [F(0), F(1, 2000000), F(1, 10000000)]:
        for mode in ['literal', 'canonical']:
            attempts.append(checked(record, h, mode))
    centers = [x for x in attempts if x['halfwidth'] == '0']
    required = [x for x in attempts if F(x['halfwidth']) == F(1, 2000000)]
    require(all(x['all_winners_verified'] for x in centers), 'center sequence fails')
    success = all(x['all_winners_verified'] for x in required)
    # Failed desired boxes remain in the receipt; smaller-box success does not
    # silently replace the requested claim.
    negatives = []
    for label, mutate in [
        ('changed unpinned source', lambda x: x['source'][0][0].__setitem__(0, '1/2')),
        ('wrong precision', lambda x: x.__setitem__('arithmetic_bits', 53)),
        ('wrong winner sequence', lambda x: x['q_nodes'][0].__setitem__('winner', 'C1')),
        ('unordered nodes', lambda x: x['q_nodes'][1].__setitem__('q', x['q_nodes'][0]['q'])),
        ('wrong halfwidth', lambda x: x.__setitem__('requested_raw_halfwidth', '0')),
    ]:
        altered = copy.deepcopy(record)
        mutate(altered)
        try:
            validate(altered)
        except (ValueError, AssertionError):
            negatives.append({'mutation': label, 'rejected': True})
        else:
            raise ValueError('malformed input accepted: '+label)
    require(source.read_bytes() == original, 'source changed during verification')
    require(hashlib.sha256(KERNEL.read_bytes()).hexdigest() == KERNEL_SHA, 'kernel changed during verification')
    require(hashlib.sha256(Path(__file__).read_bytes()).hexdigest() == own_before, 'checker changed during verification')
    result = {'status': 'PASS_DIRECTED_MPFR256_FULL64_HALFUNIT_BOX' if success else 'CENTER_PASS_REQUESTED_BOX_NOT_CERTIFIED',
              'baseline_r': '4', 'm': 16, 'criteria': 2, 'independent_raw_coordinates': 64,
              'source_record_sha256': hashlib.sha256(original).hexdigest(),
              'checker_sha256': own_before, 'kernel_sha256': KERNEL_SHA,
              'source_center_floor': '1/8',
              'requested_box_floor': str(F(1, 8)-F(1, 2000000)),
              'requested_halfwidth': '1/2000000', 'common_nu_or_mirror_constraints_on_box': False,
              'all_source_box_points_strictly_r4_admissible': True,
              'winner_sequence': PATTERN, 'integer_q_nodes': record['q_nodes'],
              'at_least_genuine_real_q_changes': 3 if success else None,
              'attempts': attempts, 'negative_tests': negatives,
              'scope': 'Finite integer-node winners and all64 independent baseline-raw coordinates. Box nu values may differ by row/column. Box floor is 0.1249995, not 1/8. No continuum-q certificate, exact event count, published data, optimum radius or empirically calibrated uncertainty.'}
    save(output, result)
    print(json.dumps({'status': result['status'], 'attempts': [
        {'mode': x['model_mode'], 'halfwidth': x['halfwidth'], 'pass': x['all_winners_verified'],
         'margins': [r['minimum_margin_lower'] for r in x['records']]} for x in attempts]}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['discover', 'verify'])
    parser.add_argument('--source', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'discover':
        discover(args.output)
    else:
        require(args.source is not None, 'source is required')
        verify(args.source, args.output)
