"""Exact rational domain and all-rung candidate/auxiliary separation checks.

No imports from either MPFR or Decimal evaluator. The source bound deductions
use monotonicity of qROFWA in the matched source powers, with positive normalized
weights. This checks that the full source box need not preserve common nu.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json


def main(input_path, output):
    raw = input_path.read_bytes()
    record = json.loads(raw)
    radius = F(record['uniform_raw_halfwidth'])
    assert record['baseline_r'] == '4' and radius > 0
    source = record['source']
    assert len(source) == 7 and all(len(row) == 2 for row in source)
    rows = []
    coordinate_bits = []
    for index, row in enumerate(source):
        assert all(len(pair) == 2 and all(isinstance(x, str) for x in pair) for pair in row)
        pairs = [[F(x) for x in pair] for pair in row]
        for mu, nu in pairs:
            assert F(1, 16) < mu-radius <= mu+radius < F(7, 8)
            assert F(1, 16) < nu-radius <= nu+radius < F(7, 8)
            assert (mu+radius)**4 + (nu+radius)**4 < 1
        lower_score = min((mu-radius)**4 for mu, nu in pairs)-max((nu+radius)**4 for mu, nu in pairs)
        upper_score = max((mu+radius)**4 for mu, nu in pairs)-min((nu-radius)**4 for mu, nu in pairs)
        if index < 3:
            assert lower_score > F(1, 4)
        else:
            assert upper_score < 0
        rows.append({'label': record['labels'][index], 'all_finite_q_ge_4_score_lower': str(lower_score),
                     'all_finite_q_ge_4_score_upper': str(upper_score)})
        coordinate_bits.extend(max(x.numerator.bit_length(), x.denominator.bit_length()) for pair in pairs for x in pair)
    result = {'status': 'PASS_EXACT_RATIONAL_DOMAIN_AND_GLOBAL_AUXILIARY_EXCLUSION',
              'source_record_sha256': hashlib.sha256(raw).hexdigest(),
              'checker_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'baseline_r': '4', 'halfwidth': str(radius), 'independent_raw_coordinates': 28,
              'r_admissible_orthopairs_checked': 14, 'maximum_source_rational_bit_length': max(coordinate_bits),
              'candidate_score_uniform_lower': '1/4', 'auxiliary_score_uniform_upper': '0',
              'row_bounds': rows,
              'scope': 'All points in the 28-coordinate source box, all finite real q>=4, any positive normalized criterion weights. Candidate pair switching itself is only certified at the seven recorded q nodes.'}
    with output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2)
    print(json.dumps({k: v for k, v in result.items() if k != 'row_bounds'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    main(args.input, args.output)
