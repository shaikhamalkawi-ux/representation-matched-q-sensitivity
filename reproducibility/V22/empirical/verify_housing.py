"""Independent ordinary Decimal60 replay of all frozen ACS housing nodes.

No import of any primary evaluator. Decimal arithmetic is not outward rounded.
Only frozen exact inputs/results are read; public path/hash routing was ported.
"""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
from decimal import Decimal, getcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time
import traceback

import numpy as np  # Read the frozen NPZ container only, not arithmetic.

D = Decimal
getcontext().prec = 60
RANK_TOL = D("1e-12")
FLOAT_TOL = D("2e-15")
COUNTS = Counter()
ERRORS = defaultdict(lambda: D(0))


def pin(path):
    data = path.read_bytes()
    return {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def dump_new(path, obj):
    with path.open("x", encoding="utf-8") as stream:
        json.dump(obj, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def check(flag, label):
    COUNTS[label] += 1
    if not bool(flag):
        raise AssertionError(label)


def close(actual, expected, label, tolerance=FLOAT_TOL):
    expected = D(str(expected))
    error = abs(actual - expected)
    ERRORS[label] = max(ERRORS[label], error)
    check(error <= tolerance, label)


def compare_tree(actual, expected, label):
    if isinstance(actual, dict):
        check(set(actual) == set(expected), label + ".keys")
        for key, value in actual.items():
            compare_tree(value, expected[key], label + "." + key)
    elif isinstance(actual, list):
        check(len(actual) == len(expected), label + ".length")
        for value, old in zip(actual, expected, strict=True):
            compare_tree(value, old, label)
    elif isinstance(actual, D):
        close(actual, expected, label,
              D("2e-13") if label.endswith("range_percentage_points") else FLOAT_TOL)
    else:
        check(actual == expected, label)


def decimalize(frac):
    return D(frac.numerator) / D(frac.denominator)


def serializable(item):
    if isinstance(item, D):
        return str(item)
    if isinstance(item, list):
        return [serializable(x) for x in item]
    if isinstance(item, dict):
        return {key: serializable(value) for key, value in item.items()}
    return item


def run(workspace, out):
    start = time.perf_counter()
    out.mkdir(parents=True, exist_ok=False)
    phase = Path(__file__).resolve().parent
    primary = phase / "expected/housing"
    source = phase / "data/housing"
    paths = [source / "EXACT_INPUTS_R1.json", source / "root_results_R1/ALL_NODE_SCORES_R1.npz"]
    paths += sorted(p for p in primary.iterdir() if p.is_file())
    workspace = phase
    guards = {str(p.relative_to(workspace)): pin(p) for p in paths}
    import portable_io
    portable_io.verify_module()
    dump_new(out / "PRE_RUN_PIN_R1.json", {
        "utc": datetime.now(timezone.utc).isoformat(), "script": pin(Path(__file__)),
        "guards": guards, "precision": 60, "rounding": str(getcontext().rounding),
        "tolerance": str(FLOAT_TOL), "rank_tolerance": str(RANK_TOL),
        "scope": "Post-hoc independent finite-grid replay; not outward or continuum certification",
        "python": sys.version,
    })
    try:
        inp = load(source / "EXACT_INPUTS_R1.json")
        group = inp["groups"]["primary"]
        rows = group["rows"]
        check(len(rows) == 51 and group["criteria"] == ["renters", "owners"], "cohort")
        copied = load(primary / "PRIMARY_EXACT_INPUTS_R1.json")
        check(copied["rows"] == rows and copied["criteria"] == group["criteria"], "copied_input")
        names = [r["name"] for r in rows]
        ids = [r["id"] for r in rows]
        nd, wy = names.index("North Dakota"), names.index("Wyoming")
        exact_mu = [[Fraction(a) for a in r["mu"]] for r in rows]
        exact_nu = [[Fraction(a) for a in r["nu"]] for r in rows]
        for aa, bb in zip(exact_mu, exact_nu, strict=True):
            for a, b in zip(aa, bb, strict=True):
                check(0 < a < 1 and 0 < b < 1 and a + b < 1, "exact_admissibility")
        mu = [[decimalize(x) for x in row] for row in exact_mu]
        nu = [[decimalize(x) for x in row] for row in exact_nu]
        log_mu = [[x.ln() for x in row] for row in mu]
        log_nu = [[x.ln() for x in row] for row in nu]
        log_complement = [[(D(1) - x).ln() for x in row] for row in mu]
        old = np.load(source / "root_results_R1/ALL_NODE_SCORES_R1.npz", allow_pickle=False)
        old_weights = old["primary__weights"]
        old_scores = old["primary__matched_raw_entropy"]
        old_q = old["q"]
        weights, scores, ranks, nodes = [], [], [], []
        qgrid = [D(64+k)/64 for k in range(961)]
        runner_ids = []
        smallest_nonzero_separation = None
        frozen_rank_differences = 0
        all_near_ties = 0
        for k, q in enumerate(qgrid):
            # Direct exact-fraction q, natural logarithm, and Decimal exponential.
            factors = []
            for j in range(2):
                total = D(0)
                for i in range(51):
                    lm, ln = log_mu[i][j]/q, log_nu[i][j]/q
                    total -= lm.exp()*lm + ln.exp()*ln
                factors.append(D(1)-total/51)
            check(all(x > 0 for x in factors), "positive_factors")
            weight = [x/sum(factors) for x in factors]
            score = []
            for i in range(51):
                membership_loss = sum(weight[j]*log_complement[i][j] for j in range(2)).exp()
                nonmembership_loss = sum(weight[j]*log_nu[i][j] for j in range(2)).exp()
                score.append(D(1)-membership_loss-nonmembership_loss)
            rank = [1+sum(other > value+RANK_TOL for other in score) for value in score]
            if k == 0:
                baseline_rank = rank[:]
            for j in range(2):
                close(weight[j], old_weights[k,j], "all_frozen_weights")
            for i in range(51):
                close(score[i], old_scores[k,i], "all_frozen_scores")
                old_rank = 1 + sum(float(other) > float(old_scores[k,i])+1e-12
                                   for other in old_scores[k])
                check(rank[i] == old_rank, "all_frozen_ranks")
            close(q, old_q[k], "frozen_q")
            runner = max((i for i in range(51) if i != nd), key=lambda i: score[i])
            runner_ids.append(runner)
            gap = score[nd]-score[runner]
            check(rank[nd] == 1 and gap > RANK_TOL, "north_dakota_strict_grid_winner")
            ordered = sorted(range(51), key=lambda i: -score[i])
            position = ordered.index(wy)
            above, below = ordered[position-1], ordered[position+1]
            displacement = max(abs(a-b) for a,b in zip(rank, baseline_rank, strict=True))
            for i in range(51):
                for j in range(i+1,51):
                    distance = abs(score[i]-score[j])
                    if distance <= RANK_TOL:
                        all_near_ties += 1
                    if distance and (smallest_nonzero_separation is None or distance < smallest_nonzero_separation):
                        smallest_nonzero_separation = distance
            nodes.append({
                "node_index": k, "q": q,
                "weight_renters": weight[0], "weight_owners": weight[1],
                "north_dakota_score": score[nd], "north_dakota_gap": gap,
                "runner_up_name": names[runner], "wyoming_score": score[wy],
                "wyoming_rank": rank[wy], "above_wyoming_name": names[above],
                "above_minus_wyoming_score": score[above]-score[wy],
                "below_wyoming_name": names[below],
                "wyoming_minus_below_score": score[wy]-score[below],
                "maximum_displacement_vs_q1": displacement,
                "normalized_maximum_displacement": D(displacement)/50,
            })
            weights.append(weight); scores.append(score); ranks.append(rank)
            if k % 160 == 0:
                print(f"Independent Decimal60 node {k}/960", flush=True)
        plotted = load(primary / "HOUSING_PLOTTED_NODES_R1.json")
        close(RANK_TOL, plotted["rank_tolerance"], "plotted_rank_tolerance", D(0))
        compare_tree(nodes, plotted["nodes"], "every_plotted_field")
        with (primary / "HOUSING_PLOTTED_NODES_R1.csv").open(encoding="utf-8", newline="") as f:
            csv_rows = list(csv.DictReader(f))
        check(len(csv_rows) == 961, "csv_rows")
        for a,b in zip(nodes,csv_rows,strict=True):
            check(set(a) == set(b), "csv_fields")
            for key,value in a.items():
                if isinstance(value, D):
                    close(value,b[key],"csv_numeric")
                elif isinstance(value,int):
                    check(value == int(b[key]),"csv_integer")
                else:
                    check(value == b[key],"csv_text")
        ranges = {}
        for j,label in enumerate(group["criteria"]):
            low = min(range(961),key=lambda k:weights[k][j])
            high = max(range(961),key=lambda k:weights[k][j])
            ranges[label] = {"minimum":weights[low][j],"first_q_at_minimum":qgrid[low],
                             "maximum":weights[high][j],"first_q_at_maximum":qgrid[high],
                             "q1":weights[0][j],"q16":weights[-1][j],
                             "range_percentage_points":100*(weights[high][j]-weights[low][j])}
        low = min(range(961),key=lambda k:nodes[k]["north_dakota_gap"])
        high = max(range(961),key=lambda k:nodes[k]["north_dakota_gap"])
        nd_summary = {"unique_winner_all961_nodes":True,
                      "runner_up_names":sorted({names[i] for i in runner_ids}),
                      "minimum_sampled_all_rival_gap":nodes[low]["north_dakota_gap"],
                      "q_at_gap_minimum":qgrid[low],
                      "maximum_sampled_all_rival_gap":nodes[high]["north_dakota_gap"],
                      "q_at_gap_maximum":qgrid[high],
                      "minimum_sampled_score":min(s[nd] for s in scores),
                      "maximum_sampled_score":max(s[nd] for s in scores)}
        def sign(x):
            return 1 if x > RANK_TOL else -1 if x < -RANK_TOL else 0
        brackets=[]
        wy_ties=0
        for rival in range(51):
            if rival == wy:
                continue
            differences=[s[wy]-s[rival] for s in scores]
            wy_ties+=sum(abs(x)<=RANK_TOL for x in differences)
            for k in range(960):
                if sign(differences[k])*sign(differences[k+1]) == -1:
                    brackets.append({"rival_id":ids[rival],"rival_name":names[rival],
                        "lower_node_q":qgrid[k],"upper_node_q":qgrid[k+1],
                        "wyoming_minus_rival_lower":differences[k],
                        "wyoming_minus_rival_upper":differences[k+1],
                        "ranks_lower":[ranks[k][wy],ranks[k][rival]],
                        "ranks_upper":[ranks[k+1][wy],ranks[k+1][rival]],
                        "wyoming_minus_rival_q1":differences[0],
                        "wyoming_minus_rival_q16":differences[-1]})
        brackets.sort(key=lambda b:(b["lower_node_q"],b["rival_name"]))
        primary_brackets=load(primary / "WYOMING_SAMPLED_CROSSING_BRACKETS_R1.json")
        compare_tree(brackets,primary_brackets["brackets"],"all_crossing_brackets")
        participating={b["rival_id"] for b in brackets}|{ids[wy],ids[nd]}
        check(primary_brackets["rival_coordinate_rows"]==[row for row in rows if row["id"] in participating],"bracket_exact_rows")
        wy_summary={"rank_q1":ranks[0][wy],"rank_q16":ranks[-1][wy],
                    "minimum_sampled_rank":min(row[wy] for row in ranks),
                    "maximum_sampled_rank":max(row[wy] for row in ranks),
                    "sampled_opposite_sign_brackets":len(brackets),"near_tied_rival_nodes":wy_ties}
        summary=load(primary / "HOUSING_MECHANISM_SUMMARY_R1.json")
        compare_tree(ranges,summary["weight_ranges_sampled"],"summary_ranges")
        compare_tree(nd_summary,summary["north_dakota"],"summary_north_dakota")
        compare_tree(wy_summary,summary["wyoming"],"summary_wyoming")
        selected=[nodes[k] for k in [0,64,192,448,960]]
        compare_tree(selected,summary["selected_nodes"],"selected_nodes")
        max_displacement=max(n["maximum_displacement_vs_q1"] for n in nodes)
        check(max_displacement==summary["maximum_displacement_vs_q1_all_rows_nodes"],"summary_max_displacement")
        close(D(max_displacement)/50,summary["normalized_maximum_displacement_vs_q1_all_rows_nodes"],"summary_normalized")
        check(guards=={str(p.relative_to(workspace)):pin(p) for p in paths},"all_guarded_inputs_unchanged")
        dump_new(out / "DECIMAL60_ALL_NODES_R1.json",serializable({"precision":60,"nodes":nodes,
                  "weights":weights,"scores":scores,"ranks":ranks,"brackets":brackets}))
        result={"status":"PASS","precision":60,"outward_rounded":False,"continuum_certificate":False,
                "model":"Unchanged primary ACS2024 r=1 matched raw-natural-log algebraic WA",
                "nodes":961,"rows":51,"independently_recomputed_weights":1922,
                "independently_recomputed_scores":49011,"rank_checks":49011,
                "weight_ranges":serializable(ranges),"north_dakota":serializable(nd_summary),
                "wyoming":wy_summary,"normalized_maximum_displacement":str(D(max_displacement)/50),
                "smallest_nonzero_pair_score_separation_at_nodes":str(smallest_nonzero_separation),
                "all_near_tied_pairs_at_nodes":all_near_ties,
                "check_counts":dict(COUNTS),"total_assertions":sum(COUNTS.values()),
                "maximum_numeric_errors":{k:str(v) for k,v in ERRORS.items()},
                "guarded_inputs_unchanged":True,"elapsed_seconds":time.perf_counter()-start,
                "limitations":["Independent arithmetic from existing exact table, not a second official-source extraction.",
                    "Ordinary Decimal60 finite-grid numerical confirmation, not rigorous outward intervals or continuum claim.",
                    "Post-hoc selected illustration, not held-out validation, clinical/housing advice, or prevalence evidence."]}
        dump_new(out / "INDEPENDENT_RESULT_R1.json",result)
        dump_new(out / "OUTPUT_MANIFEST_R1.json",{p.name:pin(p) for p in sorted(out.iterdir()) if p.is_file()})
        print(json.dumps({"status":"PASS","total_assertions":result["total_assertions"],
                          "max_score_error":str(ERRORS["all_frozen_scores"]),
                          "elapsed_seconds":result["elapsed_seconds"]}),flush=True)
    except Exception:
        dump_new(out / "FAILURE_R1.json",{"status":"FAIL","traceback":traceback.format_exc(),
                                          "check_counts":dict(COUNTS),"elapsed_seconds":time.perf_counter()-start})
        raise


if __name__ == "__main__":
    sys.dont_write_bytecode=True
    parser=argparse.ArgumentParser(description="Independent Decimal60 full-grid post-hoc housing replay")
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    assert not args.out.resolve().is_relative_to(Path(__file__).resolve().parent), "Output must be outside release module"
    run(Path(__file__).resolve().parent,args.out.resolve())
