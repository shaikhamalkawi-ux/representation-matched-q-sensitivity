"""Read-only replay optimization of the frozen Decimal90 numerical diagnostic.

This is not the original frozen-code replay, outward arithmetic, a continuum
certificate, or a new scientific experiment. The frozen check body, including
180 bisections and all result formatting, is reused without edits; only its
nested weight/scores functions are redirected to identical-argument caches.
"""
from __future__ import annotations

import argparse
import ast
import copy
from datetime import datetime, timezone
from decimal import Decimal, getcontext
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys
import time


ORIGINAL_SHA = "8f37e4b1c777d497f75342d3745d463cf8f8e050de0a9a5c297ac68b24757d7d"
EXPECTED_SHA = "4f5b4ac7fd18f5829b8b369951f7b7ca71ce68cc4b3590933f0475d57b9cd3b0"
CASE_SPECS = [
    (m, k, n, r)
    for m, k, n in [(4, 2, 2), (7, 3, 3), (8, 3, 2),
                       (12, 5, 3), (20, 2, 2), (11, 2, 3)]
    for r in (1, 4, 16)
]
OUT = Path(__file__).resolve().parent
WORKSPACE = OUT.parents[2]
SOURCE_DIR = WORKSPACE / "outputs/FSS_EIGHT_HOUR_RESEARCH_20260926/root_capacity"
ORIGINAL = SOURCE_DIR / "check_joint_construction_diagnostic.py"
EXPECTED = SOURCE_DIR / "POST_R4_JOINT_CONSTRUCTION_DIAGNOSTIC.json"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require_integrity(base):
    require(sha(ORIGINAL) == ORIGINAL_SHA,
            "Frozen original script changed during execution")
    require(sha(EXPECTED) == EXPECTED_SHA,
            "Frozen expected JSON changed during execution")
    require(sha(Path(__file__)) == base["checker_sha256"],
            "Executing checker source changed during execution")


def utc():
    return datetime.now(timezone.utc).isoformat()


def emit(value):
    print(json.dumps(value, separators=(",", ":")), flush=True)


def write_new(path, value):
    # Exclusive create: never overwrite an earlier receipt or frozen evidence.
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2)
        stream.write("\n")


class BudgetExpired(RuntimeError):
    pass


def check_time(deadline):
    if deadline is not None and time.perf_counter() >= deadline:
        raise BudgetExpired("The bounded phase time budget expired")


class CachedModel:
    """Cache constant subexpressions and identical per-q raw-log atom results."""

    def __init__(self, rows, logs, v, r, deadline=None):
        self.rows, self.logs, self.v, self.r = rows, logs, v, r
        self.deadline = deadline
        self.log_v = Decimal(v).ln()
        self.v_power = v ** r
        # Each expression below is exactly the original expression, evaluated
        # once instead of on every score call. No algebraic rewriting is used.
        self.log_u = [
            (Decimal(1-row[0]**r).ln(), Decimal(1-row[1]**r).ln())
            for row in rows
        ]
        # as_tuple() avoids treating different decimal representations as one
        # cache key even when Decimal numeric equality happens to identify them.
        self.unique_logs = []
        self.row_indices = []
        lookup = {}
        for row in logs:
            indices = []
            for log_x in row:
                key = log_x.as_tuple()
                if key not in lookup:
                    lookup[key] = len(self.unique_logs)
                    self.unique_logs.append(log_x)
                indices.append(lookup[key])
            self.row_indices.append(tuple(indices))
        self.weight_calls = 0
        self.score_calls = 0
        self.cached_weight_exp_calls = 0
        self.literal_equivalent_weight_exp_calls = 0

    def weight(self, q):
        check_time(self.deadline)
        p = Decimal(self.r) / q
        vl = p * self.log_v
        nu_atom = Decimal(vl).exp() * vl
        # Deliberately preserve LEFT ASSOCIATION:
        # original is (exp(p*log_x) * p) * log_x, NOT exp(y)*y.
        atoms = [Decimal(p*log_x).exp()*p*log_x
                 for log_x in self.unique_logs]
        # Preserve original j order, row order, per-row addition, builtin
        # sum's integer-zero start, division, and final normalization.
        e = [
            1 + sum(atoms[row[j]] + nu_atom for row in self.row_indices)
            / len(self.rows)
            for j in range(2)
        ]
        self.weight_calls += 1
        self.cached_weight_exp_calls += 1 + len(self.unique_logs)
        self.literal_equivalent_weight_exp_calls += 4 * len(self.rows)
        return e[0] / sum(e)

    def scores(self, q):
        self.score_calls += 1
        w = self.weight(q)
        return [
            1-Decimal(w*row[0]+(1-w)*row[1]).exp()-self.v_power
            for row in self.log_u
        ]

    def stats(self):
        return dict(
            unique_raw_log_representations=len(self.unique_logs),
            weight_calls=self.weight_calls,
            score_calls=self.score_calls,
            cached_weight_exp_calls=self.cached_weight_exp_calls,
            literal_equivalent_weight_exp_calls=self.literal_equivalent_weight_exp_calls,
        )


def load_original():
    require(sha(ORIGINAL) == ORIGINAL_SHA, "Frozen original script hash mismatch")
    require(sha(EXPECTED) == EXPECTED_SHA, "Frozen expected JSON hash mismatch")
    spec = importlib.util.spec_from_file_location("frozen_joint_oracle", ORIGINAL)
    module = importlib.util.module_from_spec(spec)
    # Import executes definitions, NOT main/check or any JSON write.
    old_dont_write = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = old_dont_write
    require(getcontext().prec == 90, "Original did not establish Decimal90")
    require(getcontext().rounding == "ROUND_HALF_EVEN", "Unexpected Decimal rounding")
    tree = ast.parse(ORIGINAL.read_text(encoding="utf-8"))
    check = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                 and n.name == "check")
    require(
        sum(isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
            and n.func.id == "range" and len(n.args) == 1
            and isinstance(n.args[0], ast.Constant) and n.args[0].value == 180
            for n in ast.walk(check)) == 1,
        "Expected exactly one unchanged 180-bisection loop",
    )
    return module, check


def build_functions(module, original_check):
    """AST extraction avoids manually recoding the literal oracle/check loop."""
    probe = copy.deepcopy(original_check)
    probe.name = "build_literal_probe"
    macro_index = next(
        i for i, node in enumerate(probe.body)
        if isinstance(node, ast.Assign)
        and any(isinstance(target, ast.Name) and target.id == "macro"
                for target in node.targets)
    )
    probe.body = probe.body[:macro_index] + ast.parse(
        "return rows, logs, v, delta, slopes, weight, scores"
    ).body
    # AST return is the only addition to the exact original initialization
    # and nested literal weight/scores definitions.

    optimized = copy.deepcopy(original_check)
    optimized.name = "optimized_check"
    transformed = []
    replaced = []
    for node in optimized.body:
        if isinstance(node, ast.FunctionDef) and node.name == "weight":
            transformed += ast.parse(
                "_cached = CachedModel(rows, logs, v, r, FULL_DEADLINE)\n"
                "MODELS.append(_cached)\n"
                "def weight(q):\n"
                "    return _cached.weight(q)\n"
            ).body
            replaced.append("weight")
        elif isinstance(node, ast.FunctionDef) and node.name == "scores":
            transformed += ast.parse(
                "def scores(q):\n"
                "    return _cached.scores(q)\n"
            ).body
            replaced.append("scores")
        else:
            transformed.append(node)
    require(replaced == ["weight", "scores"], "Unexpected nested function layout")
    optimized.body = transformed
    original_core = [
        ast.dump(node, include_attributes=False) for node in original_check.body
        if not (isinstance(node, ast.FunctionDef)
                and node.name in ("weight", "scores"))
    ]
    optimized_core = [
        ast.dump(node, include_attributes=False) for node in transformed
        if not (isinstance(node, ast.FunctionDef)
                and node.name in ("weight", "scores"))
        and not (isinstance(node, ast.Assign)
                 and isinstance(node.targets[0], ast.Name)
                 and node.targets[0].id == "_cached")
        and not (isinstance(node, ast.Expr)
                 and isinstance(node.value, ast.Call)
                 and isinstance(node.value.func, ast.Attribute)
                 and isinstance(node.value.func.value, ast.Name)
                 and node.value.func.value.id == "MODELS")
    ]
    require(original_core == optimized_core, "A non-cache check statement changed")
    namespace = vars(module).copy()
    namespace.update(CachedModel=CachedModel, FULL_DEADLINE=None, MODELS=[])
    derived = ast.fix_missing_locations(ast.Module(
        body=[probe, optimized], type_ignores=[]
    ))
    exec(compile(derived, "<read-only-frozen-AST-derived-check>", "exec"), namespace)
    return namespace, hashlib.sha256(
        "\n".join(original_core).encode("utf-8")
    ).hexdigest()


def same_decimal(a, b):
    return a == b and a.as_tuple() == b.as_tuple()


def equivalence_probes(namespace, budget_seconds, base):
    started = time.perf_counter()
    deadline = started + budget_seconds
    cases = []
    for case_index, spec in enumerate(CASE_SPECS, 1):
        check_time(deadline)
        m, k, n, r = spec
        build_start = time.perf_counter()
        rows, logs, v, delta, slopes, literal_w, literal_s = (
            namespace["build_literal_probe"](*spec)
        )
        literal_build_seconds = time.perf_counter() - build_start
        cache_start = time.perf_counter()
        cached = CachedModel(rows, logs, v, r, deadline)
        cache_build_seconds = time.perf_counter() - cache_start
        macro = [Decimal(r)*factor*16**ell
                 for ell in range(n) for factor in (2, 8)]
        # Exact same Decimal construction on both paths; broad deterministic
        # coverage, not an exhaustive equivalence proof.
        points = [Decimal(r), *macro]
        points += [(a+b)/2 for a, b in zip(macro, macro[1:])]
        points += [(a+2*b)/3 for a, b in zip(macro, macro[1:])]
        literal_seconds = optimized_seconds = 0.0
        outputs = []
        for q in points:
            check_time(deadline)
            start = time.perf_counter()
            w0, s0 = literal_w(q), literal_s(q)
            literal_seconds += time.perf_counter() - start
            start = time.perf_counter()
            w1, s1 = cached.weight(q), cached.scores(q)
            optimized_seconds += time.perf_counter() - start
            require(same_decimal(w0, w1), f"Weight mismatch {spec}, q={q}")
            require(len(s0) == len(s1) == m, "Score length mismatch")
            require(all(same_decimal(a, b) for a, b in zip(s0, s1)),
                    f"Score mismatch {spec}, q={q}")
            outputs.append(dict(q=str(q), weight_numeric_and_tuple_equal=True,
                                all_scores_numeric_and_tuple_equal=True))
        record = dict(
            m=m, k=k, n=n, r=r, probe_count=len(points),
            literal_build_seconds=literal_build_seconds,
            cache_build_seconds=cache_build_seconds,
            literal_evaluation_seconds=literal_seconds,
            cached_evaluation_seconds=optimized_seconds,
            measured_evaluation_speed_ratio=literal_seconds/optimized_seconds,
            cache_stats=cached.stats(), probes=outputs, passed=True,
        )
        cases.append(record)
        emit(dict(phase="literal_equivalence", case=case_index, total=18,
                  spec=spec, probes=len(points), status="PASS_EXACT_DECIMAL_AND_TUPLE"))
    require_integrity(base)
    return dict(
        **base, status="PASS_FINITE_LITERAL_EQUIVALENCE",
        phase="deterministic_original_literal_vs_cached_probes",
        elapsed_seconds=time.perf_counter()-started,
        cases=cases, total_cases=len(cases),
        total_probes=sum(c["probe_count"] for c in cases),
        assertion="Decimal numeric equality AND as_tuple equality for weights and every row score",
        limitation="Finite probes plus operation-preserving cache audit, not exhaustive numerical sampling or a continuum certificate",
    )


def full_compare(namespace, budget_seconds, expected, base):
    started = time.perf_counter()
    namespace["FULL_DEADLINE"] = started + budget_seconds
    cases = []
    timings = []
    status = "NOT_STARTED"
    failure = None
    try:
        for index, spec in enumerate(CASE_SPECS):
            check_time(namespace["FULL_DEADLINE"])
            start = time.perf_counter()
            result = namespace["optimized_check"](*spec)
            seconds = time.perf_counter()-start
            require(result == expected["cases"][index],
                    f"Retained JSON case differs at index {index}: {spec}")
            cases.append(result)
            timings.append(dict(spec=spec, elapsed_seconds=seconds,
                                cache_stats=namespace["MODELS"][-1].stats()))
            emit(dict(phase="optimized_full_compare", case=index+1, total=18,
                      spec=spec, elapsed_seconds=seconds, status="EXACT_JSON_MATCH"))
        result_payload = dict(
            status="PASS_NUMERICAL_DIAGNOSTIC_ONLY", precision_decimal_digits=90,
            not_a_continuum_or_outward_certificate=True,
            implementation="Stdlib Decimal90; full shared table and all row scores; no injected weight path",
            cases=cases,
        )
        require(result_payload == expected, "Whole retained JSON differs")
        require_integrity(base)
        status = "PASS_OPTIMIZED_FULL_READ_COMPARE"
    except BudgetExpired as exc:
        status, failure = "TIME_LIMIT_WITH_PARTIAL_PROGRESS", str(exc)
        result_payload = None
    except Exception as exc:
        status, failure = "FAIL_OPTIMIZED_COMPARE", f"{type(exc).__name__}: {exc}"
        result_payload = None
    finally:
        namespace["FULL_DEADLINE"] = None
    receipt = dict(
        **base, status=status,
        finished_utc=utc(), budget_seconds=budget_seconds,
        elapsed_seconds=time.perf_counter()-started, completed_cases=len(cases),
        case_timings=timings, failure=failure,
        expected_mode="READ_COMPARE_EXISTING_NO_OVERWRITE",
        original_full_script_executed=False,
        original_180_bisection_body_preserved=True,
        original_and_expected_hashes_unchanged=(
            sha(ORIGINAL) == ORIGINAL_SHA and sha(EXPECTED) == EXPECTED_SHA
        ),
        executing_checker_hash_unchanged=(
            sha(Path(__file__)) == base["checker_sha256"]
        ),
        integrity_checks_required_before_pass=True,
        limitation="Optimized diagnostic read/compare only; NOT original frozen-code replay, interval arithmetic, or continuum certification",
    )
    return receipt, result_payload


def main():
    start_checker_sha = sha(Path(__file__))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suffix", default="R1")
    parser.add_argument("--equivalence-budget-seconds", type=float, default=180)
    parser.add_argument("--full-budget-seconds", type=float, default=600)
    args = parser.parse_args()
    require(re.fullmatch(r"R[1-9][0-9]*", args.suffix), "Invalid receipt suffix")
    paths = {name: OUT/f"{name}_{args.suffix}.json" for name in
             ("LITERAL_EQUIVALENCE", "OPTIMIZED_FULL_COMPARE", "OPTIMIZED_RESULT")}
    require(all(not p.exists() for p in paths.values()), "Receipt exists; choose new suffix")
    module, original_check = load_original()
    namespace, body_sha = build_functions(module, original_check)
    expected = json.loads(EXPECTED.read_text(encoding="utf-8"))
    require([tuple(c[k] for k in ("m", "k", "n", "r")) for c in expected["cases"]]
            == CASE_SPECS, "Expected record case order mismatch")
    base = dict(
        started_utc=utc(), checker=str(Path(__file__).resolve()),
        checker_sha256=start_checker_sha, original_script=str(ORIGINAL),
        original_script_sha256=ORIGINAL_SHA, expected_json=str(EXPECTED),
        expected_json_sha256=EXPECTED_SHA, unchanged_check_core_ast_sha256=body_sha,
        python=sys.version, decimal_context=str(getcontext()),
        precision_decimal_digits=90, bisections_per_target=180,
        optimizations=[
            "hoist identical log(v)",
            "cache identical raw-log atom results per q with Decimal-tuple keys",
            "preserve exp(p*logx)*p*logx left association",
            "preserve row/column summation order and normalization",
            "hoist identical per-row log(1-mu**r) and v**r score constants",
        ],
        protected_inputs_are_read_only=True,
    )
    require_integrity(base)
    try:
        equivalence = equivalence_probes(namespace, args.equivalence_budget_seconds, base)
    except Exception as exc:
        failure = dict(**base, status="FAIL_OR_TIMEOUT_LITERAL_EQUIVALENCE",
                       failure=f"{type(exc).__name__}: {exc}", finished_utc=utc(),
                       full_phase_started=False)
        write_new(paths["LITERAL_EQUIVALENCE"], failure)
        emit(failure)
        return 2
    write_new(paths["LITERAL_EQUIVALENCE"], equivalence)
    emit(dict(phase="equivalence_complete", status=equivalence["status"],
              probes=equivalence["total_probes"],
              receipt=str(paths["LITERAL_EQUIVALENCE"])))
    receipt, result = full_compare(namespace, args.full_budget_seconds, expected, base)
    if result is not None:
        write_new(paths["OPTIMIZED_RESULT"], result)
        receipt["optimized_result_path"] = str(paths["OPTIMIZED_RESULT"])
        receipt["optimized_result_sha256"] = sha(paths["OPTIMIZED_RESULT"])
    receipt["literal_equivalence_receipt"] = str(paths["LITERAL_EQUIVALENCE"])
    receipt["literal_equivalence_receipt_sha256"] = sha(paths["LITERAL_EQUIVALENCE"])
    write_new(paths["OPTIMIZED_FULL_COMPARE"], receipt)
    emit(dict(status=receipt["status"], completed_cases=receipt["completed_cases"],
              elapsed_seconds=receipt["elapsed_seconds"],
              receipt=str(paths["OPTIMIZED_FULL_COMPARE"]),
              receipt_sha256=sha(paths["OPTIMIZED_FULL_COMPARE"])))
    return 0 if receipt["status"] == "PASS_OPTIMIZED_FULL_READ_COMPARE" else 2


if __name__ == "__main__":
    raise SystemExit(main())
