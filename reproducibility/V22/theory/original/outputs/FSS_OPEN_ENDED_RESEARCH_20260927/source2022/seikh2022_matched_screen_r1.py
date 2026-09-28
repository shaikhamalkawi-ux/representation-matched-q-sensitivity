"""Bounded ordinary Decimal100 matched-path screen of the 2022 HTML case.

Fixed source rung r=4, Frank tau=2. Exact rational requested grid plus analytic
equal-weight p=0 limit; no continuum certificate or PDF verification.
"""
from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "SEIKH2022_HTML_TRANSCRIPTION_R1.json"
OPERATOR = HERE / "SEIKH2022_OPERATOR_TRANSCRIPTION_R1.md"
BASELINE = HERE / "SEIKH2022_NUMERICAL_DIAGNOSTIC_R2.json"
OUTPUT = HERE / "SEIKH2022_MATCHED_SCREEN_R1.json"


def sha(x):
    return hashlib.sha256(x).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def strings(x):
    if isinstance(x, D):
        return str(x)
    if isinstance(x, dict):
        return {k: strings(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [strings(v) for v in x]
    return x


def main():
    require(not OUTPUT.exists(), "Refuse to overwrite existing matched screen")
    pinned = {path: path.read_bytes() for path in [SOURCE, OPERATOR, BASELINE, Path(__file__)]}
    source = json.loads(pinned[SOURCE])
    baseline = json.loads(pinned[BASELINE])
    with localcontext() as ctx:
        ctx.prec = 100
        raw = [[[D(x) for x in cell] for cell in row] for row in source["raw_cells"]]
        require(len(raw) == 5 and all(len(row) == 5 for row in raw), "Wrong case dimensions")
        require(source["cost_criteria_one_based"] == [1, 4, 5], "Wrong cost convention")
        normalized = [[([nu, mu] if j in (0, 3, 4) else [mu, nu])
                       for j, (mu, nu) in enumerate(row)] for row in raw]
        logs = [[[-x.ln() for x in cell] for cell in row] for row in normalized]
        canonical = [[[mu**4, nu**4] for mu, nu in row] for row in normalized]
        logtau = D(2).ln()

        def weights(p):
            if p == 0:
                return [D(1) / D(5)] * 5
            factors = [D(1) - sum((a*p) * (-a*p).exp() for row in logs for a in row[j]) / D(5)
                       for j in range(5)]
            require(all(x > 0 for x in factors), "Nonpositive raw entropy factor")
            return [x / sum(factors) for x in factors]

        def frank(values, w):
            inner = sum(wj * ((logtau*x).exp() - 1).ln() for x, wj in zip(values, w))
            return (1 + inner.exp()).ln() / logtau

        def scores(w, kind):
            result = []
            for row in canonical:
                mu = [x[0] for x in row]
                nu = [x[1] for x in row]
                if kind == "WA":
                    muq = 1 - frank([1-x for x in mu], w)
                    nuq = frank(nu, w)
                elif kind == "WG":
                    muq = frank(mu, w)
                    nuq = 1 - frank([1-x for x in nu], w)
                else:
                    raise ValueError("Unknown Frank operator")
                require(0 < muq < 1 and 0 < nuq < 1 and muq + nuq < 1, "Invalid aggregate")
                result.append((1 + muq - nuq) / 2)
            return result

        def evaluation(p):
            w = weights(p)
            out = {"p": p, "weights": w}
            for kind in ["WA", "WG"]:
                s = scores(w, kind)
                order = sorted(range(5), key=lambda i: s[i], reverse=True)
                out[kind] = {
                    "scores_printed_scale": s,
                    "descending_order_one_based": [i+1 for i in order],
                    "winner_one_based": order[0]+1,
                    "all_rival_winner_margin_printed_scale": s[order[0]]-s[order[1]],
                }
            return out

        grid = sorted(set([F(4)*F(3, 2)**k for k in range(65)] +
                          [F(q) for q in [4,5,6,8,10,12,16,20,32,64,128]]))
        evaluated = []
        for q in grid:
            qd = D(q.numerator) / D(q.denominator)
            node = evaluation(D(4) / qd)
            node["requested_q_exact_rational"] = str(q)
            node["q_decimal"] = qd
            evaluated.append(node)
        limit = evaluation(D(0))
        baseline_residuals = {}
        for kind in ["WA", "WG"]:
            old = [D(x) for x in baseline["frank_baseline_diagnostics"]["natural_log"][kind]["scores"]]
            residual = max(abs(x-y) for x,y in zip(evaluated[0][kind]["scores_printed_scale"], old))
            require(residual < D("1e-94"), "q4 did not reproduce current exact-formula baseline")
            baseline_residuals[kind] = residual
        summary = {}
        for kind in ["WA", "WG"]:
            summary[kind] = {
                "observed_finite_grid_winners": sorted(set(x[kind]["winner_one_based"] for x in evaluated)),
                "observed_adjacent_winner_changes": [
                    [a["requested_q_exact_rational"], b["requested_q_exact_rational"]]
                    for a,b in zip(evaluated, evaluated[1:])
                    if a[kind]["winner_one_based"] != b[kind]["winner_one_based"]
                ],
                "minimum_sampled_winner_margin_printed_scale": min(x[kind]["all_rival_winner_margin_printed_scale"] for x in evaluated),
                "analytic_equal_weight_limit_winner_one_based": limit[kind]["winner_one_based"],
                "limit_winner_margin_printed_scale_ordinary_evaluation": limit[kind]["all_rival_winner_margin_printed_scale"],
                "observed_finite_grid_orders": sorted(set(tuple(x[kind]["descending_order_one_based"]) for x in evaluated)),
            }
        record = {
            "status": "NUMERICAL_ONLY_BOUNDED_MATCHED_SCREEN_HTML_INPUT",
            "utc": datetime.now(timezone.utc).isoformat(),
            "precision_decimal_digits": 100,
            "input_pins": {p.name: {"sha256": sha(b), "bytes": len(b)} for p,b in pinned.items()},
            "model": {
                "source_rung": "4",
                "frank_tau": "2",
                "transport": "(mu,nu) -> (mu^(4/q),nu^(4/q))",
                "entropy": "Full recomputation of Eq5 raw natural-log factors after matching",
                "aggregation": "Literal Eq7/Eq8 with fixed canonical mu^4 and nu^4, variable recomputed weights",
                "score": "(1+mu_out^q-nu_out^q)/2; unshifted margins are twice these",
                "all_rows_retained": 5,
            },
            "grid": {
                "definition": "Union of q=4*(3/2)^k, k=0,...,64, and the listed small integer nodes",
                "additional_integer_nodes": [4,5,6,8,10,12,16,20,32,64,128],
                "finite_node_count": len(grid),
                "minimum_q": str(grid[0]),
                "maximum_q_exact_rational": str(grid[-1]),
            },
            "baseline_q4_max_abs_score_residual": baseline_residuals,
            "finite_grid": evaluated,
            "equal_weight_limit": limit,
            "limit_justification": "For each strictly positive source grade a=-ln(x) is finite; phi(a*p)->0 as p->0, so e_j->1 and w_j->1/5. Fixed canonical Frank losses are continuous in w.",
            "summary": summary,
            "limitations": [
                "No conclusion about unobserved q intervals or all-q winner stability follows from the finite grid.",
                "The limiting weights are analytic; numerical limit scores are ordinary rounded Decimal100, not an outward certificate.",
                "This is a representation-matched r4/tau2 control, not the source Section6 fixed-raw tau4 sensitivity protocol.",
                "Natural log is a declared model convention supported by Eq10 numerical agreement, not an explicit source log-base declaration.",
                "The source input is authorized HTML transcription; PDF visual verification remains pending.",
            ],
        }
        for path, data in pinned.items():
            require(path.read_bytes() == data, "Input drift: "+path.name)
        with OUTPUT.open("x", encoding="utf-8", newline="\n") as f:
            json.dump(strings(record), f, indent=2)
            f.write("\n")
        print(json.dumps(strings({"output": str(OUTPUT), "sha256": sha(OUTPUT.read_bytes()), "finite_nodes": len(grid), "summary": summary}), indent=2))


if __name__ == "__main__":
    main()
