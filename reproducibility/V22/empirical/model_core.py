"""Mechanically extracted unchanged R15 numerical/metric kernels; portable wrapper is separate."""
import hashlib
import json
import math
import platform
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
import numpy as np
METHODS=("wa","wg","cosine_topsis")
ARMS=("dynamic","frozen_q1","equal")
Q=1+np.arange(961,dtype=np.float64)/64
KEY_Q=(1,2,3,4,8,12,16)
KEY_INDEX=tuple(int((q-1)*64) for q in KEY_Q)
TOL=1e-12

def utc():
    return datetime.now(timezone.utc).isoformat()


def canonical(obj):
    return json.dumps(obj, sort_keys=True, ensure_ascii=False,
                      separators=(",", ":"), allow_nan=False).encode("utf-8")


def digest(obj):
    return hashlib.sha256(canonical(obj)).hexdigest()


def file_pin(path):
    raw = Path(path).read_bytes()
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def write_json(path, obj):
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(obj, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def write_npz(path, arrays):
    with Path(path).open("xb") as stream:
        np.savez_compressed(stream, **arrays)


def require(condition, label):
    if not condition:
        raise ValueError(label)


def entropy_weights(u, v, q):
    p = 1.0 / q
    total = np.zeros_like(u)
    for values in (u, v):
        positive = values > 0
        x = np.power(values[positive], p)
        total[positive] += x * np.log(x)
    e = 1 + total.mean(axis=0)
    require(np.all(np.isfinite(e)) and np.all(e > 0), "finite positive entropy factors")
    return e / e.sum()


def log_products(values, weights, complement=False):
    """Rows by criteria times nodes by criteria, with exact zero-factor branches."""
    good = values < 1 if complement else values > 0
    logs = np.full(values.shape, -np.inf)
    if complement:
        logs[good] = np.log1p(-values[good])
    else:
        logs[good] = np.log(values[good])
    # Every weight is strictly positive. No 0 * infinity term is possible.
    return np.exp(logs @ weights.T).T


def scores_for_weights(u, v, weights):
    weights = np.atleast_2d(weights)
    require(weights.shape[1] == u.shape[1] and np.all(weights > 0), "positive common weights")
    require(np.allclose(weights.sum(axis=1), 1, atol=2e-15, rtol=0), "normalized weights")
    wa = 1 - log_products(u, weights, complement=True) - log_products(v, weights)
    wg = log_products(u, weights) + log_products(v, weights, complement=True) - 1
    plus_cells = np.cos((math.pi / 8) * ((1 + v) ** 2 - u ** 2))
    minus_cells = np.cos((math.pi / 8) * ((1 + u) ** 2 - v ** 2))
    plus, minus = (plus_cells @ weights.T).T, (minus_cells @ weights.T).T
    denominator = plus + minus
    require(np.all(denominator > 0) and np.all(denominator >= 1 - 2e-15), "positive fixed-ideal cosine denominator")
    result = {"wa": wa, "wg": wg, "cosine_topsis": plus / denominator}
    for method, values in result.items():
        require(np.all(np.isfinite(values)), method + " finite scores")
        lower = 0 if method == "cosine_topsis" else -1
        require(np.all(values >= lower - 2e-15) and np.all(values <= 1 + 2e-15), method + " range")
    return result


def ranks(scores):
    ordered = np.sort(scores)
    return 1 + len(scores) - np.searchsorted(ordered, scores + TOL, side="right")


def top_indices(scores):
    return np.flatnonzero(scores >= scores.max() - TOL)


def pair_categories(a, b, chunk=128):
    n = len(a)
    count = {"strict_concordance": 0, "strict_discordance": 0,
             "exactly_one_tied": 0, "both_tied": 0}
    js = np.arange(n)[None, :]
    for start in range(0, n, chunk):
        end = min(start + chunk, n)
        mask = np.arange(start, end)[:, None] < js
        da = (a[start:end, None] - a[None, :])[mask]
        db = (b[start:end, None] - b[None, :])[mask]
        sa = (da > TOL).astype(np.int8) - (da < -TOL).astype(np.int8)
        sb = (db > TOL).astype(np.int8) - (db < -TOL).astype(np.int8)
        both = (sa != 0) & (sb != 0)
        count["strict_concordance"] += int(np.count_nonzero(both & (sa == sb)))
        count["strict_discordance"] += int(np.count_nonzero(both & (sa != sb)))
        count["exactly_one_tied"] += int(np.count_nonzero((sa == 0) ^ (sb == 0)))
        count["both_tied"] += int(np.count_nonzero((sa == 0) & (sb == 0)))
    total = n * (n - 1) // 2
    require(sum(count.values()) == total, "complete unordered-pair partition")
    return {"total_unordered_pairs": total, **count}


def trajectory_summary(values, ids, names, invariant=False):
    baseline_rank = ranks(values[0])
    baseline_top = tuple(map(int, top_indices(values[0])))
    runs, changes, changed_nodes = [], [], []
    max_rank_displacement = 0
    unique_transitions = 0
    previous = None
    iterable = range(1) if invariant else range(len(Q))
    for k in iterable:
        row = values[0] if invariant else values[k]
        current = tuple(map(int, top_indices(row)))
        max_rank_displacement = max(max_rank_displacement, int(np.abs(ranks(row) - baseline_rank).max()))
        if current != baseline_top:
            changed_nodes.append(k)
        if current != previous:
            if previous is not None:
                changes.append(k)
                unique_transitions += len(previous) == len(current) == 1
                runs[-1]["last_index"] = k - 1
            runs.append({"first_index": k, "last_index": k,
                         "indices": list(current), "ids": [ids[i] for i in current],
                         "names": [names[i] for i in current]})
            previous = current
    runs[-1]["last_index"] = len(Q) - 1
    adjacent = sorted({i for k in changes for i in (k - 1, k)})
    return {"q_nodes": len(Q), "invariant_broadcast": invariant,
            "q1_top_ids": [ids[i] for i in baseline_top],
            "q1_top_names": [names[i] for i in baseline_top],
            "q1_top_cardinality": len(baseline_top), "top_set_runs": runs,
            "nodes_different_top_set_from_q1": len(changed_nodes),
            "different_top_set_node_indices": changed_nodes,
            "consecutive_top_set_changes": len(changes), "top_change_right_indices": changes,
            "consecutive_unique_winner_changes": int(unique_transitions),
            "independent_change_adjacent_indices_required": adjacent,
            "max_competition_rank_displacement_from_q1": max_rank_displacement}


def exact_dominance(cohort, chunk=128):
    """Fraction-sorted ordinal coordinates give exactly equivalent comparisons."""
    fu, fv = cohort["fractions_u"], cohort["fractions_v"]
    n, d = cohort["u"].shape
    ru, rv = np.empty((n, d), dtype=np.int32), np.empty((n, d), dtype=np.int32)
    for values, ordinal in ((fu, ru), (fv, rv)):
        for j in range(d):
            distinct = sorted({row[j] for row in values})
            lookup = {x: k for k, x in enumerate(distinct)}
            ordinal[:, j] = [lookup[row[j]] for row in values]
    all_i, all_j = [], []
    for start in range(0, n, chunk):
        end = min(start + chunk, n)
        mask = np.ones((end - start, n), dtype=bool)
        strict = np.zeros_like(mask)
        for j in range(d):
            du = ru[start:end, j, None] - ru[None, :, j]
            dv = rv[start:end, j, None] - rv[None, :, j]
            mask &= (du >= 0) & (dv <= 0)
            strict |= (du > 0) | (dv < 0)
        i, j = np.nonzero(mask & strict)
        all_i.append((i + start).astype(np.int32)); all_j.append(j.astype(np.int32))
    ii, jj = np.concatenate(all_i), np.concatenate(all_j)
    same_u = np.all(ru[ii] == ru[jj], axis=1)
    same_v = np.all(rv[ii] == rv[jj], axis=1)
    masks = {}
    for key, exact, endpoint in (("u0", fu, 0), ("u1", fu, 1), ("v0", fv, 0), ("v1", fv, 1)):
        masks[key] = np.array([any(x == endpoint for x in row) for row in exact])
    wa_equal = (same_u | (masks["u1"][ii] & masks["u1"][jj])) & (same_v | (masks["v0"][ii] & masks["v0"][jj]))
    wg_equal = (same_u | (masks["u0"][ii] & masks["u0"][jj])) & (same_v | (masks["v1"][ii] & masks["v1"][jj]))
    return {"dominator_indices": ii, "dominated_indices": jj,
            "wa_exact_equal_positive_weight_products": wa_equal,
            "wg_exact_equal_positive_weight_products": wg_equal}


def dominance_summary(score, method, dominance):
    ii, jj = dominance["dominator_indices"], dominance["dominated_indices"]
    delta = score[ii] - score[jj]
    losing = delta < -TOL
    tied = np.abs(delta) <= TOL
    exact = dominance.get(method + "_exact_equal_positive_weight_products", np.zeros(len(ii), dtype=bool))
    unresolved = tied & ~exact
    return {"exact_strict_componentwise_pairs": len(ii),
            "numeric_dominator_losing_count": int(np.count_nonzero(losing)),
            "numeric_dominator_tied_within_tolerance_count": int(np.count_nonzero(tied)),
            "algebraically_confirmed_boundary_equal_count": int(np.count_nonzero(tied & exact)),
            "unresolved_near_tie_count": int(np.count_nonzero(unresolved)),
            "dominator_strictly_better_beyond_tolerance_count": int(np.count_nonzero(delta > TOL)),
            "losing_pair_offsets": np.flatnonzero(losing).tolist(),
            "unresolved_near_tie_pair_offsets": np.flatnonzero(unresolved).tolist(),
            "scope": "q1 only for numerical counts; exact WA/WG product equalities hold for every positive common weight vector. Near ties are not asserted exact."}


def evaluate_cohort(cohort, out):
    out.mkdir(parents=True, exist_ok=False)
    started = utc()
    u, v = cohort["u"], cohort["v"]
    weights = np.array([entropy_weights(u, v, q) for q in Q])
    equal = np.full(u.shape[1], 1 / u.shape[1])
    dynamic = scores_for_weights(u, v, weights)
    equally = scores_for_weights(u, v, equal)
    arrays = {"q": Q, "weights_dynamic": weights, "weights_frozen_q1": weights[0],
              "weights_equal": equal, "ids": np.array(cohort["ids"]), "names": np.array(cohort["names"])}
    report = {"cohort": cohort["key"], "started_utc": started,
              "inventory": cohort["inventory"], "names_in_order": cohort["names"],
              "ids_in_order": cohort["ids"], "methods": {}, "q1_cross_method": {},
              "q1_dominance": {}, "validation_nodes": list(KEY_Q), "tie_tolerance": TOL}
    dominance = exact_dominance(cohort)
    write_npz(out / "EXACT_DOMINANCE_R1.npz", dominance)
    for method in METHODS:
        frozen = dynamic[method][0:1].copy()
        values_by_arm = {"dynamic": dynamic[method], "frozen_q1": frozen, "equal": equally[method]}
        report["methods"][method] = {}
        report["q1_dominance"][method] = {}
        for arm, values in values_by_arm.items():
            arrays[method + "__" + arm] = values if arm == "dynamic" else values[0]
            summary = trajectory_summary(values, cohort["ids"], cohort["names"], invariant=arm != "dynamic")
            if arm == "dynamic":
                summary["seven_node_pair_comparisons_to_q1"] = [
                    {"q": q, "index": k, **pair_categories(values[0], values[k])}
                    for q, k in zip(KEY_Q, KEY_INDEX)]
                summary["max_strict_pair_discordance_SEVEN_NODES_ONLY"] = max(
                    x["strict_discordance"] for x in summary["seven_node_pair_comparisons_to_q1"])
            report["methods"][method][arm] = summary
            if arm != "dynamic":
                report["q1_dominance"][method][arm] = dominance_summary(values[0], method, dominance)
    for arm in ("frozen_q1", "equal"):
        baseline = arrays["wa__" + arm]
        report["q1_cross_method"][arm] = {}
        for method in ("wg", "cosine_topsis"):
            target = arrays[method + "__" + arm]
            report["q1_cross_method"][arm][method + "_versus_wa"] = {
                **pair_categories(baseline, target),
                "max_competition_rank_displacement": int(np.abs(ranks(baseline) - ranks(target)).max())}
    write_npz(out / "ALL_METHOD_SCORES_R1.npz", arrays)
    report["array_layout"] = {name: {"shape": list(a.shape), "dtype": str(a.dtype)} for name, a in arrays.items()}
    report["invariant_array_semantics"] = "Frozen/equal one-dimensional score arrays broadcast identically to all 961 q nodes by proved matching identities; these are not 961 independent evaluations."
    report["exact_dominance_array_pin"] = file_pin(out / "EXACT_DOMINANCE_R1.npz")
    report["score_array_pin"] = file_pin(out / "ALL_METHOD_SCORES_R1.npz")
    report["completed_utc"] = utc()
    write_json(out / "RESULTS_R1.json", report)
    print(json.dumps({"cohort": cohort["key"], "status": "COMPUTED_PENDING_INDEPENDENT_REVIEW",
                      "shape": list(u.shape), "result_sha256": file_pin(out / "RESULTS_R1.json")["sha256"]}), flush=True)
    return {"cohort": cohort["key"], "shape": list(u.shape),
            "result": file_pin(out / "RESULTS_R1.json"), "scores": report["score_array_pin"],
            "exact_dominance": report["exact_dominance_array_pin"]}


def synthetic_tests(out):
    out.mkdir(parents=True, exist_ok=False)
    checks = []
    def test(ok, label):
        require(bool(ok), "synthetic: " + label); checks.append(label)
    u = np.array([[1, 1], [0, 0], [0, 0], [.5, .25], [1, .25], [.25, 1], [.25, .25]])
    v = np.array([[0, 0], [1, 1], [0, 0], [.25, .5], [0, .5], [.5, 0], [.25, .25]])
    w = np.array([.375, .625])
    reduced = scores_for_weights(u, v, w)
    for q in KEY_Q:
        x, y = u ** (1 / q), v ** (1 / q)
        pu, pv = x ** q, y ** q
        ma = (1 - np.prod((1 - pu) ** w, axis=1)) ** (1 / q)
        na = np.prod(y ** w, axis=1)
        mg = np.prod(x ** w, axis=1)
        ng = (1 - np.prod((1 - pv) ** w, axis=1)) ** (1 / q)
        literal = {"wa": ma ** q - na ** q, "wg": mg ** q - ng ** q}
        h = 1 - pu - pv
        dp = np.sum(w * np.cos(math.pi / 4 * (np.abs(pu - 1) + np.abs(pv)) * (1 - np.abs(h) / 2)), axis=1)
        dm = np.sum(w * np.cos(math.pi / 4 * (np.abs(pu) + np.abs(pv - 1)) * (1 - np.abs(h) / 2)), axis=1)
        literal["cosine_topsis"] = dp / (dp + dm)
        for method in METHODS:
            test(np.max(np.abs(literal[method] - reduced[method][0])) < 2e-12, f"literal raw transport {method} q={q}")
        weights = entropy_weights(u, v, q)
        test(np.all(weights > 0) and abs(weights.sum() - 1) < 2e-15, f"closed entropy q={q}")
        revw = entropy_weights(u[::-1], v[::-1], q)
        test(np.max(np.abs(weights - revw)) < 2e-15, f"row permutation entropy q={q}")
    for method, values in reduced.items():
        permuted = scores_for_weights(u[:, ::-1], v[:, ::-1], w[::-1])[method]
        test(np.max(np.abs(permuted - values)) < 2e-15, method + " criterion permutation")
        test(abs(values[0, 0] - 1) < 2e-15, method + " perfect")
        test(abs(values[0, 1] - (0 if method == "cosine_topsis" else -1)) < 2e-15, method + " worst")
        test(abs(values[0, 2] - (.5 if method == "cosine_topsis" else 0)) < 2e-15, method + " zero pair")
        test(abs(values[0, 6] - (.5 if method == "cosine_topsis" else 0)) < 2e-15, method + " equal positive negative")
        repeat = scores_for_weights(np.tile(u[3], (4, 1)), np.tile(v[3], (4, 1)), w)[method]
        test(np.all(repeat == repeat[0, 0]), method + " identical rows")
    single = scores_for_weights(u[:, :1], v[:, :1], [1])
    for method in ("wa", "wg"):
        test(np.max(np.abs(single[method][0] - (u[:, 0] - v[:, 0]))) < 2e-15, method + " single criterion")
    test(np.array_equal(ranks(np.array([.2, .2, .1, .3])), [2, 2, 4, 1]), "competition ranks")
    test(pair_categories(np.array([0., 0., 1.]), np.array([0., 1., 0.])) == {
        "total_unordered_pairs": 3, "strict_concordance": 0, "strict_discordance": 1,
        "exactly_one_tied": 2, "both_tied": 0}, "pair partition with ties")
    fixture = {"u": u, "fractions_u": [[Fraction(float(x)) for x in row] for row in u],
               "fractions_v": [[Fraction(float(x)) for x in row] for row in v]}
    dom = exact_dominance(fixture)
    brute = {(i, j) for i in range(len(u)) for j in range(len(u))
             if np.all(u[i] >= u[j]) and np.all(v[i] <= v[j])
             and (np.any(u[i] > u[j]) or np.any(v[i] < v[j]))}
    test(set(zip(dom["dominator_indices"].tolist(), dom["dominated_indices"].tolist())) == brute,
         "exact Fraction-ordinal dominance")
    for method in METHODS:
        summary = dominance_summary(reduced[method][0], method, dom)
        test(summary["numeric_dominator_losing_count"] == 0, method + " synthetic Pareto")
        if method == "cosine_topsis":
            test(summary["numeric_dominator_tied_within_tolerance_count"] == 0, "cosine strict synthetic Pareto")
    report = {"status": "PASS_SYNTHETIC_ONLY", "created_utc": utc(), "checks": checks,
              "count": len(checks), "real_inputs_loaded": False, "real_models_evaluated": False,
              "code": file_pin(__file__), "python": platform.python_version(), "numpy": np.__version__,
              "scope": "Closed-boundary/formula/identity synthetic checks only; not primary-source example reproduction or real-data admission."}
    write_json(out / "SYNTHETIC_TESTS_R1.json", report)
    print(json.dumps({"status": report["status"], "checks": len(checks), "code": report["code"]}))
