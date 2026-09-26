"""Targeted independent regression checks of the V17 scalar entropy enclosure.

These finite checks supplement the analytic derivative/range argument; they do
not prove the implementation correct. Optional mpmath checks are sanity checks,
not directed-rounding certificates. Run with the project bundled Python -B.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import random
import sys
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from seikh_certificates import certify, domain, entropy_term, scores
from dyadic_interval import IV, PREC, SCALE

try:
    import mpmath as mp
except ImportError:
    mp = None


def point_value(n: int) -> IV:
    """Evaluate a dyadic point without invoking the new range function."""
    point = IV(n, n)
    return point * point.exp()


def endpoint_pair(value: IV) -> tuple[int, int]:
    return value.lo, value.hi


class EntropyRangeTests(unittest.TestCase):
    def assertContains(self, outer: IV, inner: IV) -> None:
        self.assertLessEqual(outer.lo, inner.lo)
        self.assertGreaterEqual(outer.hi, inner.hi)

    def test_exact_zero(self):
        self.assertEqual(endpoint_pair(entropy_term(IV.integer(0))), (0, 0))

    def test_dyadic_singletons(self):
        for value in (F(-3), F(-2), F(-1), F(-1, 2), F(-1, 8)):
            with self.subTest(value=str(value)):
                interval = IV.frac(value)
                self.assertEqual(interval.lo, interval.hi)
                self.assertEqual(endpoint_pair(entropy_term(interval)),
                                 endpoint_pair(point_value(interval.lo)))

    def test_singleton_critical_value(self):
        critical = -IV.integer(-1).exp()
        self.assertEqual(endpoint_pair(entropy_term(IV.integer(-1))),
                         endpoint_pair(critical))

    def test_wholly_decreasing_interval(self):
        y = IV.frac(-3, -2)
        left, right = point_value(y.lo), point_value(y.hi)
        # Both points are below -1, hence h(left)>h(right).
        self.assertGreater(left.lo, right.hi)
        self.assertEqual(endpoint_pair(entropy_term(y)), (right.lo, left.hi))

    def test_wholly_increasing_interval(self):
        y = IV.frac(F(-3, 4), F(-1, 4))
        left, right = point_value(y.lo), point_value(y.hi)
        self.assertLess(left.hi, right.lo)
        self.assertEqual(endpoint_pair(entropy_term(y)), (left.lo, right.hi))

    def test_strict_critical_straddle(self):
        y = IV.frac(-2, F(-1, 2))
        critical = -IV.integer(-1).exp()
        left, right = point_value(y.lo), point_value(y.hi)
        self.assertLess(critical.hi, min(left.lo, right.lo))
        self.assertEqual(entropy_term(y).lo, critical.lo)
        self.assertEqual(entropy_term(y).hi, max(left.hi, right.hi))

    def test_touch_critical_from_left(self):
        y = IV.frac(-2, -1)
        self.assertEqual(endpoint_pair(entropy_term(y)),
                         ((-IV.integer(-1).exp()).lo, point_value(y.lo).hi))

    def test_touch_critical_from_right(self):
        y = IV.frac(-1, F(-1, 2))
        self.assertEqual(endpoint_pair(entropy_term(y)),
                         ((-IV.integer(-1).exp()).lo, point_value(y.hi).hi))

    def test_one_ulp_critical_straddles_and_neighbors(self):
        critical = -IV.integer(-1).exp()
        for a, b in ((-SCALE-1, -SCALE+1), (-SCALE-1, -SCALE),
                     (-SCALE, -SCALE+1), (-SCALE-2, -SCALE-1),
                     (-SCALE+1, -SCALE+2)):
            with self.subTest(a=a, b=b):
                y = IV(a, b)
                result = entropy_term(y)
                self.assertContains(result, point_value(a))
                self.assertContains(result, point_value(b))
                if a <= -SCALE <= b:
                    self.assertLessEqual(result.lo, critical.lo)
                    self.assertGreaterEqual(result.hi, critical.hi)

    def test_zero_as_right_endpoint(self):
        result = entropy_term(IV.frac(-2, 0))
        self.assertEqual(result.hi, 0)
        self.assertEqual(result.lo, (-IV.integer(-1).exp()).lo)

    def test_seeded_exact_dyadic_interior_points(self):
        rng = random.Random(170926)
        count = 0
        for _ in range(64):
            a_units = rng.randrange(-256, -1)
            b_units = rng.randrange(a_units + 1, 1)
            y = IV.frac(F(a_units, 16), F(b_units, 16))
            enclosure = entropy_term(y)
            for numerator in (0, 1, 2, 3, 4):
                point = y.lo + (y.hi - y.lo) * numerator // 4
                self.assertContains(enclosure, point_value(point))
                count += 1
        self.assertEqual(count, 320)

    def test_nonpositive_range(self):
        for a, b in ((-10, -5), (-3, -2), (-2, 0), (-1, 0)):
            self.assertLessEqual(entropy_term(IV.frac(a, b)).hi, 0)

    def test_reject_positive_singleton(self):
        with self.assertRaises(ValueError):
            entropy_term(IV(1, 1))

    def test_reject_interval_crossing_positive(self):
        with self.assertRaises(ValueError):
            entropy_term(IV(-SCALE, 1))

    def test_reject_outside_exp_implementation_domain(self):
        with self.assertRaises(ValueError):
            entropy_term(IV.integer(-1025))

    def test_reject_reversed_interval_at_constructor(self):
        with self.assertRaises(ValueError):
            IV(0, -1)

    def test_compactification_endpoint_uniform_weights(self):
        inputs = domain('seikh', 'source', F(1, 100))
        result, weights = scores(IV.integer(0), inputs)
        fifth = F(1, 5)
        self.assertEqual(len(result), 4)
        self.assertEqual(len(weights), 5)
        self.assertEqual(len(set(endpoint_pair(w) for w in weights)), 1)
        for w in weights:
            self.assertLessEqual(F(w.lo, SCALE), fifth)
            self.assertGreaterEqual(F(w.hi, SCALE), fifth)
            self.assertGreater(w.lo, 0)

    def test_reject_t_outside_declared_parameter_domain(self):
        inputs = domain('seikh', 'source', F(0))
        for t in (IV(-1, 0), IV(0, SCALE // 4 + 1)):
            with self.subTest(t=endpoint_pair(t)):
                with self.assertRaises(ValueError):
                    scores(t, inputs)

    def test_reject_invalid_source_radius(self):
        for radius in (F(-1, 100), F(1, 20), F(1, 5)):
            with self.subTest(radius=str(radius)):
                with self.assertRaises(ValueError):
                    domain('seikh', 'source', radius)

    def test_reject_fixed_q_below_four(self):
        for q in (F(-1), F(1), F(3), F(399, 100)):
            with self.subTest(q=str(q)):
                with self.assertRaises(ValueError):
                    certify(F(1, 100), q=q)

    @unittest.skipIf(mp is None, 'mpmath not installed; optional sanity check')
    def test_high_precision_finite_point_sanity(self):
        count = 0
        with mp.workdps(210):
            scale = mp.mpf(SCALE)
            rng = random.Random(170927)
            intervals = [IV(-SCALE-1, -SCALE+1), IV(-SCALE, -SCALE),
                         IV.frac(-2, 0), IV.frac(-3, -2),
                         IV.frac(F(-3, 4), F(-1, 4))]
            for _ in range(32):
                a_units = rng.randrange(-256, -1)
                b_units = rng.randrange(a_units + 1, 1)
                intervals.append(IV.frac(F(a_units, 16), F(b_units, 16)))
            for y in intervals:
                outer = entropy_term(y)
                lower, upper = mp.mpf(outer.lo)/scale, mp.mpf(outer.hi)/scale
                points = {y.lo, y.hi, (y.lo+y.hi)//2}
                if y.lo <= -SCALE <= y.hi:
                    points.add(-SCALE)
                for n in points:
                    x = mp.mpf(n)/scale
                    value = x*mp.exp(x)
                    self.assertLessEqual(lower, value)
                    self.assertGreaterEqual(upper, value)
                    count += 1
        self.assertGreaterEqual(count, 100)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(EntropyRangeTests)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    report = {
        'status': 'PASS' if result.wasSuccessful() else 'FAIL',
        'tests_run': result.testsRun,
        'failure_count': len(result.failures),
        'error_count': len(result.errors),
        'skipped_count': len(result.skipped),
        'failures': [{'test': str(test), 'traceback': trace}
                     for test, trace in result.failures + result.errors],
        'precision_bits': PREC,
        'mpmath_version': None if mp is None else mp.__version__,
        'mpmath_role': 'finite 210-decimal sanity checks only, not a proof',
        'analytic_reference': 'design_review.md: exact entropy summand range',
        'implementation_sha256': hashlib.sha256(
            (HERE/'seikh_certificates.py').read_bytes()).hexdigest(),
        'test_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'review_scope': 'AI-assisted internal review; not a proof-assistant verification',
    }
    if args.output:
        with args.output.open('x', encoding='utf-8') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')
    print(json.dumps(report, indent=2))
    raise SystemExit(0 if result.wasSuccessful() else 1)


if __name__ == '__main__':
    main()
