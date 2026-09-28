"""Independent Decimal100 enclosure replay of the exact floor16 witness.

Only the separately audited Decimal primitive is imported. No MPFR evaluator,
discovery model, entropy-gap reduction or common-nu cancellation is imported.
All 64 raw coordinates vary independently in the nonzero box.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SOURCE = ROOT/'novelty_challenge/source_floor/compact/FLOOR16_SOURCE_AND_INTEGER_NODES_R1.json'
SOURCE_SHA = '059d8456e553f2e0e8ccb894ac75b58ca99635045ab42b1017d30ad251e1d79f'
RECORD = ROOT/'novelty_challenge/source_floor/compact/FLOOR16_FULL64_MPFR256_R1.json'
RECORD_SHA = 'd72356e5a37dc057d622e04d05147fbbfa317d79efc0f1fb06463e389412327c'
KERNEL = ROOT/'entropy_family/pinar_independent/replay_decimal.py'
KERNEL_SHA = '82fc5a524072a0da7dac34aafa53e3b96db39df727d42c3bc8158b3adfd36235'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


if digest(KERNEL) != KERNEL_SHA:
    raise ValueError('Decimal primitive identity')
spec = importlib.util.spec_from_file_location('floor16_independent_decimal', KERNEL)
kernel = importlib.util.module_from_spec(spec)
spec.loader.exec_module(kernel)
Box, require = kernel.Box, kernel.require
LABELS = [f'C{i}' for i in range(1, 9)]+[f'D{i}' for i in range(1, 9)]
NODES = [('34', 'C8'), ('142', 'C7'), ('238', 'C6'), ('721', 'C5')]


def product(values):
    result = Box(1)
    for value in values:
        result = result*value
    return result


def validate(data):
    require(data['schema'] == 'floor16-integer-node-source-v1', 'schema')
    require(data['baseline_r'] == '4' and data['labels'] == LABELS, 'model')
    require([(n['q'], n['winner']) for n in data['q_nodes']] == NODES, 'exact nodes')
    require(data['candidate_decimal_places'] == 6, 'decimal precision')
    require(F(data['requested_raw_halfwidth']) == F(1, 2000000), 'halfunit')
    source = [[[F(x) for x in pair] for pair in row] for row in data['source']]
    require(len(source) == 16 and all(len(row) == 2 for row in source), 'matrix shape')
    require(all(len(pair) == 2 for row in source for pair in row), 'pair shape')
    require(all(F(1, 8) <= x < 1 and (x*10**6).denominator == 1
                for row in source for pair in row for x in pair), 'center floor and grid')
    require(all(source[i] == source[7-i][::-1] for i in range(8)), 'center mirrors')
    require(all(row == [[F(1, 8), F(1, 8)], [F(1, 4), F(1, 8)]]
                for row in source[8:]), 'auxiliary center')
    require(all(nu == F(1, 8) for row in source for mu, nu in row), 'center nu')
    canonical = json.dumps(data['source'], sort_keys=True, separators=(',', ':')).encode()
    require(hashlib.sha256(canonical).hexdigest() == data['source_sha256'], 'embedded source digest')
    return source


def check_rounding(source):
    delta = (Box(F(1, 4))*Box(F(-1, 4)).exp()
             -Box(F(1, 6))*Box(F(-1, 6)).exp())/8
    checked = []
    for i in range(8):
        slope = F(2*(i+1)-9, 28)
        for j, weight in enumerate([F(1), F(0)]):
            ell = -Box(2).log()+slope*(weight-F(1, 2))+delta*slope*slope
            mu = (1-ell.exp()).power(F(1, 4))
            rounded_lo = F(round(F(mu.lo)*10**6), 10**6)
            rounded_hi = F(round(F(mu.hi)*10**6), 10**6)
            require(rounded_lo == rounded_hi == source[i][j][0], 'independent nearest rounding')
            h = F(1, 2000000)
            require(source[i][j][0]-h < F(mu.lo) <= F(mu.hi) < source[i][j][0]+h,
                    'unrounded design strictly inside halfunit cell')
            checked.append({'candidate': LABELS[i], 'criterion': j+1,
                            'design_membership': mu.pack(),
                            'six_decimal_exact': str(source[i][j][0])})
    return checked


def boxed_source(source, halfwidth):
    boxes = []
    maxima = []
    count = 0
    for row in source:
        boxed_row = []
        for mu, nu in row:
            ml, mh = mu-halfwidth, mu+halfwidth
            nl, nh = nu-halfwidth, nu+halfwidth
            require(0 < ml <= mh < 1 and 0 < nl <= nh < 1, 'strict interior box')
            require(mh**4+nh**4 < 1, 'every box pair strictly admissible')
            maxima.append(mh**4+nh**4)
            boxed_row.append((Box(ml, mh), Box(nl, nh)))
            count += 2
        boxes.append(boxed_row)
    require(count == 64, 'independent coordinate count')
    return boxes, {'independent_raw_coordinates': count,
                   'minimum_raw_box_endpoint': str(min(x-halfwidth for row in source for pair in row for x in pair)),
                   'maximum_upper_corner_powered_pair_sum': str(max(maxima)),
                   'common_nu_or_mirror_equalities_imposed': False}


def evaluate(boxes, q, mode):
    p = F(4)/q
    # Reconstruct transported raw pairs and raw-entropy weights literally.
    raw = [[(mu.power(p), nu.power(p)) for mu, nu in row] for row in boxes]
    factors = []
    for j in range(2):
        total = Box(0)
        for row in raw:
            mu, nu = row[j]
            total = total+mu*mu.log()+nu*nu.log()
        factors.append(1+total/16)
    require(all(e.lo > 0 for e in factors), 'positive entropy factors')
    denominator = factors[0]+factors[1]
    weights = [e/denominator for e in factors]
    require(all(0 < w.lo <= w.hi < 1 for w in weights), 'interior weights')
    scores = []
    for index in range(16):
        if mode == 'literal':
            transported = raw[index]
            complement = product((1-mu.power(q)).power(w)
                                 for (mu, nu), w in zip(transported, weights))
            agg_mu = (1-complement).power(1/q)
            agg_nu = product(nu.power(w) for (mu, nu), w in zip(transported, weights))
            scores.append(agg_mu.power(q)-agg_nu.power(q))
        elif mode == 'canonical':
            original = boxes[index]
            member_product = product((1-mu.power(4)).power(w)
                                     for (mu, nu), w in zip(original, weights))
            nonmember_product = product(nu.power(4).power(w)
                                        for (mu, nu), w in zip(original, weights))
            scores.append(1-member_product-nonmember_product)
        else:
            raise ValueError('unrecognized formula')
    return scores, factors, weights


def main(output):
    own_hash = digest(Path(__file__))
    require(digest(SOURCE) == SOURCE_SHA and digest(RECORD) == RECORD_SHA, 'source/record pins')
    data = json.loads(SOURCE.read_text(encoding='utf8'))
    mpfr = json.loads(RECORD.read_text(encoding='utf8'))
    require(mpfr['status'] == 'PASS_DIRECTED_MPFR256_FULL64_HALFUNIT_BOX', 'MPFR claimed status')
    source = validate(data)
    rounding = check_rounding(source)
    primitive_checks = kernel.selftests()
    attempts = []
    score_containment = []
    for h in [F(0), F(1, 2000000)]:
        boxes, domain = boxed_source(source, h)
        for mode in ['literal', 'canonical']:
            reference = next(x for x in mpfr['attempts']
                             if F(x['halfwidth']) == h and x['model_mode'] == mode)
            records = []
            for node_index, (qs, label) in enumerate(NODES):
                scores, factors, weights = evaluate(boxes, F(qs), mode)
                winner = LABELS.index(label)
                gaps = [scores[winner]-score for i, score in enumerate(scores) if i != winner]
                require(all(g.lo > 0 for g in gaps), 'winner not certified by independent enclosure')
                contained = []
                for i, (score, expected) in enumerate(zip(scores, reference['records'][node_index]['all_scores'])):
                    lo, hi = map(F, expected['exact_dyadic_endpoints'])
                    contained.append(lo <= F(score.lo) <= F(score.hi) <= hi)
                    score_containment.append(contained[-1])
                records.append({'q': qs, 'winner': label,
                                'minimum_global_margin_lower': str(min(g.lo for g in gaps)),
                                'all_scores': [s.pack() for s in scores],
                                'all_rival_gaps': [g.pack() for g in gaps],
                                'entropy_factors': [e.pack() for e in factors],
                                'weights': [w.pack() for w in weights],
                                'score_enclosures_inside_MPFR_record': contained})
            attempts.append({'halfwidth': str(h), 'model_mode': mode, 'domain': domain,
                             'all_four_global_winners_verified': True, 'records': records})
            print(json.dumps({'halfwidth': str(h), 'model': mode,
                              'minimum_margin_lower': min(F(x['minimum_global_margin_lower']) for x in records).__str__(),
                              'scores_inside_MPFR': sum(sum(x['score_enclosures_inside_MPFR_record']) for x in records)}), flush=True)
    # Containment is a cross-record diagnostic, not assumed by the positivity proof.
    require(digest(SOURCE) == SOURCE_SHA and digest(RECORD) == RECORD_SHA
            and digest(KERNEL) == KERNEL_SHA and digest(Path(__file__)) == own_hash,
            'input or checker changed during replay')
    result = {'status': 'PASS_INDEPENDENT_DECIMAL100_LITERAL_AND_CANONICAL_FULL64_BOX',
              'source_sha256': SOURCE_SHA, 'MPFR_record_sha256': RECORD_SHA,
              'primitive_sha256': KERNEL_SHA, 'checker_sha256': own_hash,
              'python': sys.version, 'decimal_precision': 100,
              'primitive_selftest_count': primitive_checks,
              'source_candidate_rounding_checks': rounding,
              'source_center_floor': '1/8', 'box_floor': str(F(1, 8)-F(1, 2000000)),
              'halfwidth': '1/2000000', 'integer_nodes': NODES,
              'scores_compared_to_MPFR': len(score_containment),
              'scores_contained_in_MPFR': sum(score_containment),
              'attempts': attempts,
              'scope': 'One exact six-decimal synthetic table, four integer nodes, at least three forced real-rung winner changes for every point in the full64 raw box. The box floor is 0.1249995, not 1/8. No whole-q certificate, exact switch count, maximal radius, empirical calibration, generic six-decimal guarantee, or public benchmark claim.'}
    with output.open('x', encoding='utf8') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'status': result['status'], 'scores_compared': len(score_containment),
                      'scores_contained': sum(score_containment)}, indent=2), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    main(parser.parse_args().output)
