"""Compact exact-rational multiwinner witness and full 28-coordinate box.

Discovery freezes a source, ordered rational rung nodes and halfwidth. Verify
then evaluates the literal transported qROFWA pipeline with directed MPFR256.
No equality constraints on the perturbed memberships/nonmemberships are used.
Only a finite strict-winner pattern is certified, not exact event count/radius.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import copy
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
KERNEL = ROOT / 'outputs/FSS_STRENGTHENING_20260926/radius_refinement/mpfr_audit/verify_refined_mpfr.py'
KERNEL_SHA = '867b871950cf41960ee1274411b2581cda88377ecdc512617991c1471487592a'
if hashlib.sha256(KERNEL.read_bytes()).hexdigest() != KERNEL_SHA:
    raise ValueError('MPFR arithmetic kernel hash mismatch')
spec = importlib.util.spec_from_file_location('compact_multiwinner_mpfr', KERNEL)
kernel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(kernel)
I, require = kernel.I, kernel.require
R = F(4)
LABELS = ['C1', 'C2', 'C3', 'D1', 'D2', 'D3', 'D4']
AUXILIARIES = [(F(16157, 250000), F(8, 125)),
               (F(19313, 250000), F(2, 25)),
               (F(2099, 20000), F(1, 10)),
               (F(24221, 200000), F(31, 250))]


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':')).encode()


def power(x, exponent):
    exponent = exponent if isinstance(exponent, I) else I(exponent)
    return (x.log() * exponent).exp()


def product(values):
    ans = I(1)
    for value in values:
        ans *= value
    return ans


def pack(x):
    return {'lower_decimal': str(x.lo), 'upper_decimal': str(x.hi),
            'exact_dyadic_endpoints': x.pair()}


def source_box(source, halfwidth):
    require(len(source) == 7 and all(len(row) == 2 for row in source), '7 by 2 source')
    require(halfwidth >= 0, 'negative halfwidth')
    box = []
    for row in source:
        b = []
        for pair in row:
            require(len(pair) == 2, 'orthopair length')
            mu, nu = map(F, pair)
            require(F(1, 16) < mu - halfwidth <= mu + halfwidth < F(7, 8), 'compact membership box')
            require(F(1, 16) < nu - halfwidth <= nu + halfwidth < F(7, 8), 'compact nonmembership box')
            require((mu + halfwidth)**4 + (nu + halfwidth)**4 < 1, 'entire source box r4 admissible')
            b.append((I(mu - halfwidth, mu + halfwidth), I(nu - halfwidth, nu + halfwidth)))
        box.append(b)
    return box


def raw_and_weights(box, q):
    q = F(q)
    require(q >= R, 'rung below baseline')
    raw = [[(power(mu, R/q), power(nu, R/q)) for mu, nu in row] for row in box]
    entropy = [1 + sum(mu * mu.log() + nu * nu.log()
                       for mu, nu in [row[j] for row in raw]) / 7 for j in range(2)]
    require(all(e.lo > 0 for e in entropy), 'positive entropy factors')
    weights = [e/sum(entropy) for e in entropy]
    require(all(w.lo > 0 and w.hi < 1 for w in weights), 'interior weights')
    return raw, entropy, weights


def evaluate(box, q):
    q = F(q)
    raw, entropy, weights = raw_and_weights(box, q)
    scores = []
    for row in raw:
        # Literal target-q raw aggregator, followed by its q-powered score.
        complement = product(power(1 - power(mu, q), w) for (mu, nu), w in zip(row, weights))
        agg_mu = power(1 - complement, 1/q)
        agg_nu = product(power(nu, w) for (mu, nu), w in zip(row, weights))
        scores.append(power(agg_mu, q) - power(agg_nu, q))
    return scores, entropy, weights


def checked_nodes(source, nodes, halfwidth):
    box = source_box(source, halfwidth)
    rows = []
    for node in nodes:
        q = F(node['q'])
        winner = LABELS.index(node['winner'])
        scores, entropy, weights = evaluate(box, q)
        gaps = [scores[winner] - value for i, value in enumerate(scores) if i != winner]
        rows.append({'q': str(q), 'winner': node['winner'],
                     'winner_verified': all(g.lo > 0 for g in gaps),
                     'minimum_margin_lower': str(min(g.lo for g in gaps)),
                     'all_scores': [pack(s) for s in scores],
                     'all_rival_gaps': [pack(g) for g in gaps],
                     'entropy_factors': [pack(e) for e in entropy],
                     'weights': [pack(w) for w in weights]})
    return rows


def candidates(d, digits):
    rows = []
    for c in [F(-1, 4), F(0), F(1, 4)]:
        row = []
        for w in [F(1), F(0)]:
            ell = -I(2).log() + c*(w-F(1, 2)) + d*c*c
            mu = power(1 - ell.exp(), F(1, 4))
            lo = round(kernel.rat(mu.lo)*10**digits)
            hi = round(kernel.rat(mu.hi)*10**digits)
            require(lo == hi, 'rational freezing unresolved')
            row.append([str(F(lo, 10**digits)), '1/8'])
        rows.append(row)
    require(rows[0] == rows[2][::-1] and rows[1][0] == rows[1][1], 'frozen seed symmetry')
    return rows + [[[str(a), '1/8'], [str(b), '1/8']] for a, b in AUXILIARIES]


def central_node(source, left, right):
    box = source_box(source, F(0))
    def direction(q):
        _, _, weights = raw_and_weights(box, q)
        displacement = weights[0] - F(1, 2)
        return 1 if displacement.lo > 0 else -1 if displacement.hi < 0 else 0
    lo, hi = F(left), F(right)
    sign = direction(lo)
    require(sign*direction(hi) == -1, 'initial entropy crossing bracket')
    for _ in range(40):
        mid = (lo+hi)/2
        mid_sign = direction(mid)
        if not mid_sign:
            lo = hi = mid
            break
        if mid_sign == sign:
            lo = mid
        else:
            hi = mid
    # Discovery only: this exact rational point is subsequently fully verified.
    return F(round(((lo+hi)/2)*10**10), 10**10)


def discover(output):
    endpoints = [{'q': str(q), 'winner': label} for q, label in
                 [(4, 'C1'), (5, 'C3'), (6, 'C1'), (8, 'C3')]]
    attempts = []
    chosen = None
    for exponent in range(6, 13):
        d = F(1, 10**exponent)
        source = candidates(d, 16)
        rows = checked_nodes(source, endpoints, F(0))
        ok = all(row['winner_verified'] for row in rows)
        attempts.append({'d': str(d), 'pass': ok,
                         'endpoint_margin_lowers': [row['minimum_margin_lower'] for row in rows]})
        if ok:
            chosen = (d, source)
            break
    require(chosen is not None, 'no finite source design admitted')
    d, source = chosen
    nodes = []
    for index, node in enumerate(endpoints):
        nodes.append(node)
        if index < 3:
            nodes.append({'q': str(central_node(source, F(node['q']), F(endpoints[index+1]['q']))), 'winner': 'C2'})
    centre_rows = checked_nodes(source, nodes, F(0))
    require(all(row['winner_verified'] for row in centre_rows), 'centre pattern rejected')
    radius_attempts = []
    selected_radius = None
    for exponent in range(7, 25):
        radius = F(1, 10**exponent)
        rows = checked_nodes(source, nodes, radius)
        ok = all(row['winner_verified'] for row in rows)
        radius_attempts.append({'halfwidth': str(radius), 'pass': ok,
                               'minimum_margin_lowers': [row['minimum_margin_lower'] for row in rows]})
        if ok:
            selected_radius = radius
            break
    require(selected_radius is not None, 'no positive full-coordinate box admitted')
    record = {'schema': 'compact-multiwinner-raw-box-v1', 'baseline_r': '4',
              'arithmetic_bits': 256, 'source': source, 'labels': LABELS,
              'source_sha256': hashlib.sha256(canonical_bytes(source)).hexdigest(),
              'source_design_d': str(d), 'candidate_decimal_places': 16,
              'q_nodes': nodes, 'uniform_raw_halfwidth': str(selected_radius),
              'metric': 'baseline raw-coordinate L-infinity, all 28 mu/nu coordinates independent',
              'rung_window': ['4', '8'], 'minimum_genuine_changes': 6,
              'source_design_attempts': attempts, 'box_search_attempts': radius_attempts,
              'interval_kernel_sha256': KERNEL_SHA,
              'scope': 'Constructed fixed exact-rational seven-row example; positive box preserves finite strict real-q winner sequence. Not a published centre, exact event count, optimum radius, useful error calibration, or universal unequal-nu complexity theorem.'}
    with output.open('x', encoding='utf-8') as stream:
        json.dump(record, stream, indent=2)
    print(json.dumps({key: record[key] for key in ['source_sha256', 'source_design_d', 'q_nodes', 'uniform_raw_halfwidth']}, indent=2))


def verify_record(record):
    require(record['schema'] == 'compact-multiwinner-raw-box-v1', 'schema')
    require(record['baseline_r'] == '4' and record['arithmetic_bits'] == 256, 'baseline/precision')
    require(record['labels'] == LABELS and record['rung_window'] == ['4', '8'], 'labels/window')
    require(record['minimum_genuine_changes'] == 6, 'event count')
    require(record['interval_kernel_sha256'] == KERNEL_SHA, 'kernel pin')
    require(record['source_sha256'] == hashlib.sha256(canonical_bytes(record['source'])).hexdigest(), 'source hash')
    radius = F(record['uniform_raw_halfwidth'])
    require(radius > 0, 'positive source-box width')
    nodes = record['q_nodes']
    require(len(nodes) == 7 and [node['winner'] for node in nodes] == ['C1', 'C2', 'C3', 'C2', 'C1', 'C2', 'C3'], 'seven-node winner pattern')
    qs = [F(node['q']) for node in nodes]
    require(qs[0] == 4 and qs[-1] == 8 and all(a < b for a, b in zip(qs, qs[1:])), 'chronological window')
    centre = checked_nodes(record['source'], nodes, F(0))
    box = checked_nodes(record['source'], nodes, radius)
    require(all(row['winner_verified'] for row in centre+box), 'strict literal-pipeline winner not enclosed')
    return {'status': 'PASS_DIRECTED_MPFR256_FULL_28_COORDINATE_BOX',
            'centre_records': centre, 'box_records': box,
            'uniform_halfwidth': str(radius), 'minimum_genuine_real_q_changes': 6,
            'common_nu_or_symmetry_constraints_on_box': False,
            'all_source_box_points_strictly_r4_admissible': True,
            'raw_coordinate_bounds_for_entire_box': ['1/16', '7/8']}


def verify(input_path, output):
    raw = input_path.read_bytes()
    record = json.loads(raw)
    result = verify_record(record)
    negatives = []
    for label, mutate in [
        ('wrong winner', lambda x: x['q_nodes'][1].__setitem__('winner', 'C1')),
        ('wrong precision', lambda x: x.__setitem__('arithmetic_bits', 53)),
        ('zero radius', lambda x: x.__setitem__('uniform_raw_halfwidth', '0')),
        ('changed unpinned source', lambda x: x['source'][0][0].__setitem__(0, '1/2')),
        ('unordered q nodes', lambda x: x['q_nodes'][1].__setitem__('q', '4')),
        ('invalid large box', lambda x: x.__setitem__('uniform_raw_halfwidth', '1/10')),
    ]:
        bad = copy.deepcopy(record)
        mutate(bad)
        try:
            verify_record(bad)
        except (ValueError, AssertionError):
            negatives.append({'mutation': label, 'rejected': True})
        else:
            raise ValueError('malformed record accepted: ' + label)
    result.update({'input_record_sha256': hashlib.sha256(raw).hexdigest(),
                   'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   'interval_kernel_sha256': KERNEL_SHA, 'negative_tests': negatives,
                   'source': record['source'], 'q_nodes': record['q_nodes'],
                   'scope': record['scope']})
    with output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2)
    print(json.dumps({'status': result['status'], 'halfwidth': result['uniform_halfwidth'],
                      'box_margin_lowers': [row['minimum_margin_lower'] for row in result['box_records']],
                      'negative_tests_rejected': len(negatives)}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['discover', 'verify'])
    parser.add_argument('--input', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.action == 'discover':
        discover(args.output)
    else:
        require(args.input is not None, '--input needed')
        verify(args.input, args.output)
