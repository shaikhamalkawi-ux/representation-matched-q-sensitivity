"""Ordinary Decimal100 diagnostic R2 of the separately transcribed 2022 case.

R1 is preserved. R2 additionally pins the saved operator transcription and
reports post-hoc truncation compatibility without attributing that convention.

No interval certificate or source-PDF verification is claimed. The literal
Eq7/Eq8 Frank formulas were supplied by root after direct authorized HTML read.
All input data are from this 5x5 case, not the different 2023 benchmark.
Output creation is exclusive; previous diagnostics are never overwritten.
"""
from decimal import Decimal as D, localcontext, ROUND_HALF_EVEN, ROUND_DOWN
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

HERE = Path(__file__).resolve().parent
INPUT = HERE / "SEIKH2022_HTML_TRANSCRIPTION_R1.json"
OPERATOR = HERE / "SEIKH2022_OPERATOR_TRANSCRIPTION_R1.md"
OUT = HERE / "SEIKH2022_NUMERICAL_DIAGNOSTIC_R2.json"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def stringify(x):
    if isinstance(x, D):
        return str(x)
    if isinstance(x, dict):
        return {k: stringify(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [stringify(v) for v in x]
    return x


def compare(values, printed):
    residuals = [x - y for x, y in zip(values, printed)]
    unit = D("0.0001")
    return {
        "residuals_computed_minus_printed": residuals,
        "max_absolute_residual": max(abs(x) for x in residuals),
        "rounding_convention_for_diagnostic_only": "nearest four decimals, ties to even",
        "computed_rounded_4dp": [x.quantize(unit, rounding=ROUND_HALF_EVEN) for x in values],
        "all_equal_when_rounded_to_printed_4dp": all(
            x.quantize(unit, rounding=ROUND_HALF_EVEN) == y
            for x, y in zip(values, printed)
        ),
        "each_within_4dp_half_unit": [abs(x) <= D("0.00005") for x in residuals],
        "posthoc_truncation_diagnostic_not_author_attribution": {
            "computed_truncated_toward_zero_4dp": [x.quantize(unit, rounding=ROUND_DOWN) for x in values],
            "all_equal_to_printed": all(x.quantize(unit, rounding=ROUND_DOWN) == y for x, y in zip(values, printed)),
        },
    }


def main():
    require(not OUT.exists(), "Refuse to overwrite existing R2 diagnostic")
    source_bytes = INPUT.read_bytes()
    operator_bytes = OPERATOR.read_bytes()
    source = json.loads(source_bytes)
    script_bytes = Path(__file__).read_bytes()
    with localcontext() as ctx:
        ctx.prec = 100
        raw = [[[D(x) for x in cell] for cell in row] for row in source["raw_cells"]]
        require(len(raw) == 5 and all(len(row) == 5 for row in raw), "Expected literal 5x5 case")
        require(all(len(cell) == 2 for row in raw for cell in row), "Expected 50 coordinates")
        require(source["baseline"] == {"q": 4, "tau": 2}, "Unexpected source baseline")
        q, tau = D(4), D(2)
        require(all(0 < x < 1 for row in raw for cell in row for x in cell), "Raw grade outside (0,1)")
        require(all(mu**4 + nu**4 < 1 for row in raw for mu, nu in row), "Source is not strictly q4 admissible")
        costs = source["cost_criteria_one_based"]
        require(costs == [1, 4, 5], "Unexpected source cost convention")
        normalized = [[([nu, mu] if j + 1 in costs else [mu, nu])
                       for j, (mu, nu) in enumerate(row)] for row in raw]
        mean_raw_log_terms = [
            sum(x * x.ln() for row in raw for x in row[j]) / D(5)
            for j in range(5)
        ]
        normalized_terms = [
            sum(x * x.ln() for row in normalized for x in row[j]) / D(5)
            for j in range(5)
        ]
        require(mean_raw_log_terms == normalized_terms, "Symmetric entropy must be swap invariant")
        printed_weights = [D(x) for x in source["printed_entropy_weights_eq10"]]
        require(sum(printed_weights) == 1, "Printed Eq10 weights do not sum to one; do not silently renormalize")
        variants = {"printed_eq10_as_written": {"weights": printed_weights}}
        log_denominators = {"natural_log": D(1), "base10_log": D(10).ln(), "base2_log": D(2).ln()}
        entropy = {}
        for name, logbase in log_denominators.items():
            factors = [D(1) + s / logbase for s in mean_raw_log_terms]
            require(all(x > 0 for x in factors), "A tested factor is nonpositive")
            weights = [e / sum(factors) for e in factors]
            entropy[name] = {
                "mean_raw_log_natural_terms": mean_raw_log_terms,
                "log_base_conversion_denominator": logbase,
                "unnormalized_raw_factors": factors,
                "normalized_weights": weights,
                "comparison_to_eq10": compare(weights, printed_weights),
            }
            variants[name] = {"weights": weights}

        logtau = tau.ln()

        def log_frank_product(canonical, weights):
            # log_tau(1 + product_j (tau^x_j - 1)^w_j).
            inner = sum(w * ((logtau * x).exp() - D(1)).ln()
                        for x, w in zip(canonical, weights))
            return (D(1) + inner.exp()).ln() / logtau

        def aggregate(row, weights, operator):
            mup = [mu**4 for mu, nu in row]
            nup = [nu**4 for mu, nu in row]
            if operator == "WA":
                muq = D(1) - log_frank_product([D(1) - x for x in mup], weights)
                nuq = log_frank_product(nup, weights)
            elif operator == "WG":
                muq = log_frank_product(mup, weights)
                nuq = D(1) - log_frank_product([D(1) - x for x in nup], weights)
            else:
                raise ValueError("Unknown declared Frank operator")
            require(0 < muq < 1 and 0 < nuq < 1 and muq + nuq < 1, "Aggregate outside strict qROF domain")
            return {
                "mu_power_q": muq,
                "nu_power_q": nuq,
                "aggregate_mu": (muq.ln() / q).exp(),
                "aggregate_nu": (nuq.ln() / q).exp(),
                "printed_scale_score": (D(1) + muq - nuq) / D(2),
                "unshifted_powered_difference": muq - nuq,
            }

        baseline = {}
        for name, variant in variants.items():
            weights = variant["weights"]
            require(abs(sum(weights) - 1) < D("1e-95"), "Unexpected weight normalization drift")
            result = {"weights_used": weights, "weight_sum": sum(weights)}
            for operator, key in [
                ("WA", "printed_frank_wa_entropy_scores_table4"),
                ("WG", "printed_frank_wg_entropy_scores_table4"),
            ]:
                aggregates = [aggregate(row, weights, operator) for row in normalized]
                scores = [row["printed_scale_score"] for row in aggregates]
                result[operator] = {
                    "rows": aggregates,
                    "scores": scores,
                    "descending_order": sorted(range(1, 6), key=lambda i: scores[i-1], reverse=True),
                    "comparison_to_table4": compare(scores, [D(x) for x in source[key]]),
                }
            baseline[name] = result
        record = {
            "status": "ORDINARY_DECIMAL100_HTML_SOURCE_DIAGNOSTIC",
            "utc": datetime.now(timezone.utc).isoformat(),
            "precision_decimal_digits": 100,
            "input_name": INPUT.name,
            "input_sha256": sha(source_bytes),
            "operator_provenance_name": OPERATOR.name,
            "operator_provenance_sha256": sha(operator_bytes),
            "checker_sha256": sha(script_bytes),
            "source_status": source["status"],
            "doi": source["doi"],
            "baseline": {"q": "4", "tau": "2", "alternatives": 5, "criteria": 5, "raw_coordinates": 50},
            "normalization": {"operation": "swap membership and nonmembership only", "cost_criteria_one_based": costs},
            "normalized_cells": normalized,
            "entropy_before_after_cost_swap_identical": True,
            "source_eq7_WA": "mu^q=1-log_tau(1+prod(tau^(1-mu_j^q)-1)^w_j); nu^q=log_tau(1+prod(tau^(nu_j^q)-1)^w_j)",
            "source_eq8_WG": "mu^q=log_tau(1+prod(tau^(mu_j^q)-1)^w_j); nu^q=1-log_tau(1+prod(tau^(1-nu_j^q)-1)^w_j)",
            "formula_provenance": "Exact Eq7/Eq8 supplied by root from direct authorized ProQuest full-text HTML read; PDF visual verification pending",
            "printed_eq10_weight_sum": sum(printed_weights),
            "entropy_log_base_diagnostics": entropy,
            "frank_baseline_diagnostics": baseline,
            "limitations": [
                "Ordinary rounded Decimal100 arithmetic, not directed intervals or a continuum certificate.",
                "The transcribed source is HTML, not yet visually verified against the PDF.",
                "The inspected text does not explicitly specify a log base; numerical agreement is evidence, not an author declaration.",
                "No source weights or formulas changed to force score agreement.",
                "No AHP renormalization, 2023 benchmark import, matched transport or parameter fitting is performed.",
                "Four-decimal nearest-rounding and post-hoc truncation comparisons are diagnostic conventions, not a claim about the authors' internal rounding.",
            ],
        }
        require(INPUT.read_bytes() == source_bytes, "Input changed while computing")
        require(OPERATOR.read_bytes() == operator_bytes, "Operator provenance changed while computing")
        require(Path(__file__).read_bytes() == script_bytes, "Checker changed while computing")
        output = json.dumps(stringify(record), indent=2) + "\n"
        with OUT.open("x", encoding="utf-8", newline="\n") as f:
            f.write(output)
        print(json.dumps({
            "output": str(OUT),
            "sha256": sha(OUT.read_bytes()),
            "entropy_4dp_matches": {k: v["comparison_to_eq10"]["all_equal_when_rounded_to_printed_4dp"] for k, v in entropy.items()},
            "baseline_4dp_matches": {k: {op: v[op]["comparison_to_table4"]["all_equal_when_rounded_to_printed_4dp"] for op in ("WA", "WG")} for k, v in baseline.items()},
            "baseline_truncation_compatibility": {k: {op: v[op]["comparison_to_table4"]["posthoc_truncation_diagnostic_not_author_attribution"]["all_equal_to_printed"] for op in ("WA", "WG")} for k, v in baseline.items()},
        }, indent=2))


if __name__ == "__main__":
    main()

