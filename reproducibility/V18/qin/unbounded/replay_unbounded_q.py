"""Independent 256-bit replay of the t=1/q closed cover of [0,1/3].

This file imports neither the generator nor its helpers. It reconstructs the
claimed input cube and rational tree, and re-evaluates every retained leaf.
The shared native interval kernel and MPFR/GMP remain explicit trust boundaries.
"""
from __future__ import annotations

import argparse
import copy
import ctypes
from fractions import Fraction as Q
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time

HERE = Path(__file__).resolve().parent
SOURCES = ("mpfr_interval_kernel.cpp", "continuous_q_kernel.cpp",
           "compact_q_kernel.cpp", "certify_unbounded_q.py")


class InvalidCertificate(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise InvalidCertificate(message)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_binary(value):
    require(type(value) in (float, int) and math.isfinite(value),
            "finite numeric endpoint required")
    return Q.from_float(float(value))


def exact_text(value):
    require(isinstance(value, str), "rational value must be a string")
    number = Q(value)
    require(str(number) == value, "rational string is not canonical")
    return number


def outward(number, upper):
    value = float(number)
    represented = Q.from_float(value)
    if (upper and represented < number) or (not upper and represented > number):
        value = math.nextafter(value, math.inf if upper else -math.inf)
    return value


def floor12(number):
    scaled = number.numerator * 10**12 // number.denominator
    return ("-" if scaled < 0 else "") + f"{abs(scaled)//10**12}.{abs(scaled)%10**12:012d}"


def scores_and_margins(scores):
    require(isinstance(scores, list) and len(scores) == 5, "five score intervals required")
    rational = []
    for pair in scores:
        require(isinstance(pair, list) and len(pair) == 2, "score interval has wrong shape")
        low, high = map(exact_binary, pair)
        require(low <= high, "reversed score interval")
        rational.append((low, high))
    return [rational[0][0] - rational[j][1] for j in range(1, 5)]


class Native:
    def __init__(self, library, runtime_dirs):
        self.handles = []
        if os.name == "nt":
            for directory in (library.parent, *runtime_dirs):
                self.handles.append(os.add_dll_directory(str(directory.resolve())))
        self.lib = ctypes.CDLL(str(library.resolve()))
        pointer = ctypes.POINTER(ctypes.c_double)
        self.lib.kernel_precision.argtypes = [ctypes.c_uint]
        self.lib.kernel_precision.restype = None
        self.lib.kernel_precision(256)
        self.lib.kernel_scores_compact.argtypes = [pointer, pointer, ctypes.c_double,
                                                   ctypes.c_double, ctypes.c_int, pointer]
        self.lib.kernel_scores_compact.restype = ctypes.c_int
        self.lib.kernel_scores_limit.argtypes = [pointer, pointer, ctypes.c_int, pointer]
        self.lib.kernel_scores_limit.restype = ctypes.c_int
        self.lib.kernel_error.restype = ctypes.c_char_p
        self.lib.kernel_version.restype = ctypes.c_char_p
        self.calls = 0
        self.version = self.lib.kernel_version().decode("ascii")

    def evaluate(self, lo, hi, tlo, thi, merge, limit=False):
        lower = None if lo is None else (ctypes.c_double * 200)(*lo)
        upper = None if hi is None else (ctypes.c_double * 200)(*hi)
        output = (ctypes.c_double * 10)()
        if limit:
            status = self.lib.kernel_scores_limit(lower, upper, merge, output)
        else:
            status = self.lib.kernel_scores_compact(
                lower, upper, outward(tlo, False), outward(thi, True), merge, output)
        self.calls += 1
        require(status == 0, "native evaluation rejected domain: " + self.lib.kernel_error().decode())
        intervals = [[output[2*i], output[2*i+1]] for i in range(5)]
        scores_and_margins(intervals)
        return intervals


def interval_from_identifier(identifier):
    require(isinstance(identifier, str) and re.fullmatch(r"r[01]{0,64}", identifier),
            "invalid tree identifier")
    lower, upper = Q(0), Q(1, 3)
    for digit in identifier[1:]:
        midpoint = (lower + upper) / 2
        if digit == "0":
            upper = midpoint
        else:
            lower = midpoint
    return lower, upper


def check_row(row):
    require(isinstance(row, dict), "node must be an object")
    lower, upper = interval_from_identifier(row["id"])
    require(exact_text(row["t_lower_exact"]) == lower and
            exact_text(row["t_upper_exact"]) == upper,
            "node interval does not match exact rational subdivision")
    require(type(row["depth"]) is int and row["depth"] == len(row["id"])-1,
            "node depth does not match identifier")
    require(type(row["positive"]) is bool, "positive flag is not boolean")
    if "error" in row:
        require(row["positive"] is False and isinstance(row["error"], str),
                "error node falsely marked positive")
        require("score_intervals" not in row and "margin_lowers_exact" not in row,
                "error node has incompatible score payload")
        return None
    margins = scores_and_margins(row["score_intervals"])
    require(row["margin_lowers_exact"] == list(map(str, margins)),
            "stored margin differs from exact endpoint subtraction")
    require(row["positive"] == (min(margins) > 0), "stored positivity inconsistent")
    return margins


def structure(document, library, merge, radius, library_policy):
    require(document["status"] == "PASS" and document["cover_complete"] is True,
            "certificate not marked complete PASS")
    require(document["winner"] == "A1" and document["all_four_challengers"] is True,
            "wrong winner/challenger scope")
    require(type(document["merge"]) is int and document["merge"] == merge and
            document["model"] == ("H" if merge == 0 else "A"), "model differs from expected model")
    require(document["q_min_exact"] == "3" and document["t_domain_exact"] == ["0", "1/3"],
            "wrong q/t domain")
    require(type(document["bits"]) is int and 64 <= document["bits"] <= 256,
            "unsupported generating precision")
    require(document["unresolved_leaves"] == 0 and document["unresolved"] == [],
            "unresolved leaves cannot be certified")
    require(exact_text(document["canonical_radius_exact"]) == radius and radius > 0,
            "radius differs from independently requested positive radius")
    require(set(document["source_hashes"]) == set(SOURCES), "wrong source manifest")
    for filename in SOURCES:
        require(document["source_hashes"][filename] == digest(HERE / filename),
                "source hash mismatch: " + filename)
    require(document["library_sha256"] == library_policy["retained_sha256"],
            "retained library identity changed from original certificate")
    require(digest(library) == library_policy["evaluated_sha256"], "loaded library file changed during replay")
    if not library_policy["allow_rebuilt"]:
        require(document["library_sha256"] == library_policy["evaluated_sha256"],
                "native library hash mismatch; rebuilt binaries require explicit --allow-rebuilt-library")

    # Independent reconstruction of the requested exact printed-input cube.
    source = (HERE / SOURCES[0]).read_text(encoding="utf-8")
    match = re.search(r"DATA\[200\]\s*=\s*\{([^}]+)\}", source)
    require(match is not None, "source DATA array absent")
    data = [int(token.strip()) for token in match.group(1).split(",")]
    require(len(data) == 200 and all(0 < value < 10 for value in data), "invalid source DATA")
    centres = [Q(value**3, 1000) for value in data]
    expected_lo = [outward(c-radius, False) for c in centres]
    expected_hi = [outward(c+radius, True) for c in centres]
    lo, hi = document["input_intervals"]["lower"], document["input_intervals"]["upper"]
    require(isinstance(lo, list) and isinstance(hi, list) and len(lo) == len(hi) == 200,
            "200 complete input intervals required")
    require(lo == expected_lo and hi == expected_hi, "input cube differs from exact requested cube enclosure")
    for c, low, high in zip(centres, lo, hi):
        require(0 < exact_binary(low) <= c-radius <= c+radius <= exact_binary(high) < 1,
                "input interval lacks positivity or complete cube inclusion")
    slack = min(Q(1)-exact_binary(hi[i])-exact_binary(hi[i+1]) for i in range(0, 200, 2))
    require(slack > 0, "enclosing cube is not strictly inside the canonical simplices")

    leaves, nodes = document["leaf_cover"], document["all_nodes"]
    require(isinstance(leaves, list) and leaves and isinstance(nodes, list), "missing tree/cover")
    require(type(document["positive_leaves"]) is int and document["positive_leaves"] == len(leaves),
            "positive leaf count mismatch")
    require(type(document["nodes"]) is int and document["nodes"] == len(nodes) == 2*len(leaves)-1,
            "complete binary-tree node count mismatch")
    by_id = {}
    for row in nodes:
        check_row(row)
        require(row["id"] not in by_id, "duplicate tree identifier")
        by_id[row["id"]] = row
    require("r" in by_id, "root absent")
    leaf_ids, margins = set(), []
    for row in leaves:
        current = check_row(row)
        require(current is not None and row["positive"] is True, "nonpositive retained leaf")
        require(row["id"] not in leaf_ids, "duplicate retained leaf")
        require(by_id.get(row["id"]) == row, "leaf record differs from full tree")
        leaf_ids.add(row["id"])
        margins.extend(current)
    for identifier, row in by_id.items():
        if identifier != "r":
            require(identifier[:-1] in by_id, "orphan tree node")
        children = (identifier+"0", identifier+"1")
        if identifier in leaf_ids:
            require(not any(child in by_id for child in children), "retained leaf has descendants")
        else:
            require(row["positive"] is False and all(child in by_id for child in children),
                    "nonleaf lacks both exact children")
    lower_edges = [exact_text(row["t_lower_exact"]) for row in leaves]
    require(lower_edges == sorted(lower_edges), "leaf cover is not ordered")
    require(lower_edges[0] == 0 and exact_text(leaves[-1]["t_upper_exact"]) == Q(1, 3),
            "cover does not span closed interval [0,1/3]")
    require(all(exact_text(left["t_upper_exact"]) == exact_text(right["t_lower_exact"])
                for left, right in zip(leaves, leaves[1:])), "gap/interior overlap in exact cover")
    require(document["uniform_margin_lower_exact"] == str(min(margins)) and
            document["uniform_margin_floor_12dp"] == floor12(min(margins)),
            "stored uniform margin inconsistent")
    return leaves, lo, hi, {"minimum_input_lower_exact": str(min(map(exact_binary, lo))),
                            "minimum_simplex_slack_exact": str(slack)}


def replay(document, library, merge, radius, native, library_policy, progress=False):
    leaves, lo, hi, domain = structure(document, library, merge, radius, library_policy)
    minima = []
    for index, row in enumerate(leaves):
        scores = native.evaluate(lo, hi, exact_text(row["t_lower_exact"]),
                                 exact_text(row["t_upper_exact"]), merge)
        for stored, fresh in zip(row["score_intervals"], scores):
            require(stored[0] <= fresh[0] <= fresh[1] <= stored[1],
                    "fresh 256-bit enclosure escapes stored interval at " + row["id"])
        margins = scores_and_margins(scores)
        require(min(margins) > 0, "fresh leaf has nonpositive winner-margin lower bound")
        minima.extend(margins)
        if progress and (index+1) % 32 == 0:
            print(f"Replayed {index+1}/{len(leaves)} retained leaves", flush=True)
    return {"positive_leaves_replayed": len(leaves), "fresh_margin_lower_exact": str(min(minima)),
            "fresh_margin_floor_12dp": floor12(min(minima)), "input_domain_checks": domain}


def negative_controls(document, library, merge, radius, native, library_policy):
    mutations = [
        ("missing_leaf", lambda d: d["leaf_cover"].pop()),
        ("changed_t_endpoint", lambda d: d["leaf_cover"][0].__setitem__("t_lower_exact", "1/10000")),
        ("forged_margin", lambda d: d["leaf_cover"][0]["margin_lowers_exact"].__setitem__(0, "1")),
        ("changed_source_hash", lambda d: d["source_hashes"].__setitem__(SOURCES[2], "0"*64)),
        ("changed_library_hash", lambda d: d.__setitem__("library_sha256", "0"*64)),
        ("wrong_radius", lambda d: d.__setitem__("canonical_radius_exact", str(radius*2))),
        ("shrunk_input_cube", lambda d: d["input_intervals"]["lower"].__setitem__(0, d["input_intervals"]["upper"][0])),
        ("missing_tree_child", lambda d: d["all_nodes"].pop()),
        ("unresolved_node", lambda d: d["unresolved"].append(d["leaf_cover"][0])),
        ("wrong_model", lambda d: d.__setitem__("merge", 1-merge)),
        ("changed_domain", lambda d: d.__setitem__("t_domain_exact", ["0", "1/4"])),
        ("nonfinite_input", lambda d: d["input_intervals"]["upper"].__setitem__(0, math.inf)),
    ]
    records = []
    for name, mutate in mutations:
        altered = copy.deepcopy(document)
        mutate(altered)
        try:
            structure(altered, library, merge, radius, library_policy)
        except (InvalidCertificate, KeyError, TypeError, IndexError, ZeroDivisionError) as error:
            records.append({"test": name, "rejected": True, "reason": str(error)})
        else:
            raise InvalidCertificate("negative control unexpectedly accepted: " + name)

    # Alter both leaf/tree copies and all summaries consistently. Structural
    # checking alone must pass; fresh numerical replay must reject the forgery.
    altered = copy.deepcopy(document)
    target = altered["leaf_cover"][0]
    low, high = target["score_intervals"][0]
    target["score_intervals"][0][0] = math.nextafter((low+high)/2, high)
    target["margin_lowers_exact"] = list(map(str, scores_and_margins(target["score_intervals"])))
    for index, row in enumerate(altered["all_nodes"]):
        if row["id"] == target["id"]:
            altered["all_nodes"][index] = copy.deepcopy(target)
    uniform = min(Q(value) for row in altered["leaf_cover"] for value in row["margin_lowers_exact"])
    altered["uniform_margin_lower_exact"] = str(uniform)
    altered["uniform_margin_floor_12dp"] = floor12(uniform)
    structure(altered, library, merge, radius, library_policy)
    try:
        replay(altered, library, merge, radius, native, library_policy)
    except InvalidCertificate as error:
        records.append({"test": "self_consistent_forged_narrow_score", "rejected": True,
                        "reason": str(error), "structural_check_passed": True})
    else:
        raise InvalidCertificate("fresh replay accepted unjustifiably narrowed score interval")
    return records


def endpoint_checks(document, merge, native):
    lo, hi = document["input_intervals"]["lower"], document["input_intervals"]["upper"]
    compact = native.evaluate(lo, hi, Q(0), Q(0), merge)
    baseweight = native.evaluate(lo, hi, Q(0), Q(0), merge, limit=True)
    margin_compact, margin_limit = map(scores_and_margins, (compact, baseweight))
    require(min(margin_compact) > 0 and min(margin_limit) > 0,
            "t=0 / direct limit cube evaluation failed to certify winner")
    require(all(max(a[0], b[0]) <= min(a[1], b[1]) for a, b in zip(compact, baseweight)),
            "compact and direct-limit enclosures are disjoint")
    return {"t_zero_full_cube_scores": compact, "baseweight_limit_full_cube_scores": baseweight,
            "t_zero_margin_floor_12dp": floor12(min(margin_compact)),
            "baseweight_limit_margin_floor_12dp": floor12(min(margin_limit)),
            "enclosure_overlap_is_only_implementation_diagnostic": True,
            "equality_justification": "analytic positive-input t=0 argument in independent audit; not interval overlap"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--runtime-dir", type=Path, action="append", default=[])
    parser.add_argument("--expected-merge", choices=(0, 1), type=int, required=True)
    parser.add_argument("--expected-radius", type=Q, required=True)
    parser.add_argument("--negative-controls", action="store_true")
    parser.add_argument("--allow-rebuilt-library", action="store_true",
                        help="explicit fresh validation with a rebuilt binary; all source hashes remain mandatory")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    def reject_constant(value):
        raise InvalidCertificate("nonfinite JSON value: " + value)
    document = json.loads(args.certificate.read_text(encoding="utf-8"), parse_constant=reject_constant)
    retained_hash = document["library_sha256"]
    require(isinstance(retained_hash, str) and re.fullmatch(r"[0-9a-f]{64}", retained_hash),
            "original certificate has invalid retained-library hash")
    library_policy = {"retained_sha256": retained_hash, "evaluated_sha256": digest(args.library),
                      "allow_rebuilt": args.allow_rebuilt_library}
    native = Native(args.library, args.runtime_dir)
    result = replay(document, args.library, args.expected_merge, args.expected_radius,
                    native, library_policy, progress=True)
    result.update(status="PASS", bits=256, mpfr_version=native.version,
                  expected_model="H" if args.expected_merge == 0 else "A",
                  expected_radius_exact=str(args.expected_radius),
                  certificate_sha256=digest(args.certificate), checker_sha256=digest(Path(__file__)),
                  library_sha256=digest(args.library),
                  validation_mode="rebuilt-library-fresh-validation" if args.allow_rebuilt_library else "retained-binary-replay",
                  historical_library_sha256=retained_hash,
                  evaluated_library_sha256=library_policy["evaluated_sha256"],
                  binary_identity_matches_historical=retained_hash == library_policy["evaluated_sha256"],
                  rebuilt_library_explicitly_allowed=args.allow_rebuilt_library,
                  source_hashes={name: digest(HERE/name) for name in SOURCES},
                  coverage="complete exact closed rational-halving cover of t in [0,1/3]",
                  finite_q_consequence="all finite real q>=3, by bijection t=1/q",
                  L2_ball_included=f"canonical L2 ball of radius {args.expected_radius} is contained in the checked L-infinity cube")
    result["endpoint_checks"] = endpoint_checks(document, args.expected_merge, native)
    if args.negative_controls:
        result["negative_controls"] = negative_controls(document, args.library, args.expected_merge,
                                                         args.expected_radius, native, library_policy)
    result["runtime_library_hashes"] = {
        str(path.resolve()): digest(path)
        for directory in args.runtime_dir
        for filename in ("libmpfr-6.dll", "libgmp-10.dll", "libwinpthread-1.dll")
        if (path := directory / filename).exists()}
    result["total_native_evaluations"] = native.calls
    result["elapsed_seconds"] = time.perf_counter()-started
    result["trust_boundary"] = [
        "same native interval kernel and arithmetic libraries as generator, not independent arithmetic implementation",
        "source hashes establish identity, not proof that binary implements source or source transcribes publisher",
        "analytic interval inclusion and positive-input limit argument are reviewed, not proof-assistant verified",
        "point/enclosure cross-checks are implementation diagnostics, not proof of formula equality"]
    if args.allow_rebuilt_library:
        result["trust_boundary"].append(
            "explicit rebuilt-binary validation: historical certificate bytes/hash retained; not historical binary replay")
    encoded = json.dumps(result, indent=2, allow_nan=False) + "\n"
    args.output.write_text(encoded, encoding="utf-8")
    print(encoded)


if __name__ == "__main__":
    main()
