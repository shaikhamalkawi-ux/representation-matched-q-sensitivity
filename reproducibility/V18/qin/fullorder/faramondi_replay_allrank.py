"""Standalone full-order cover checker; imports no generation/replay helpers.

Requires Python's standard library and a compiled continuous_q_kernel library.
All model/source paths are explicit CLI arguments. Rebuilt-library validation
is labelled separately from replay of the originally recorded native binary.
"""
from __future__ import annotations

import argparse
import copy
import ctypes
from fractions import Fraction as F
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time

CPP_NAMES = ("continuous_q_kernel.cpp", "mpfr_interval_kernel.cpp")
HISTORICAL_PY_NAMES = ("certify_continuous_q.py", "replay_joint_certificate.py")


class InvalidCertificate(ValueError):
    pass


def require(value, message):
    if not value:
        raise InvalidCertificate(message)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_endpoint(value):
    require(type(value) in (int, float) and math.isfinite(value), "finite real endpoint required")
    return F.from_float(float(value))


def floor12(value):
    integer = value.numerator * 10**12 // value.denominator
    return ("-" if integer < 0 else "") + f"{abs(integer)//10**12}.{abs(integer)%10**12:012d}"


def outward(value, upper):
    result = float(value)
    if (upper and F.from_float(result) < value) or (not upper and F.from_float(result) > value):
        result = math.nextafter(result, math.inf if upper else -math.inf)
    return result


def interval_for_id(identifier, root):
    require(isinstance(identifier, str) and re.fullmatch(r"r[01]*", identifier), "invalid dyadic node id")
    low, high = root
    for bit in identifier[1:]:
        middle = (low+high)/2
        if bit == "0":
            high = middle
        else:
            low = middle
    return low, high


def check_scores(scores):
    require(isinstance(scores, list) and len(scores) == 5, "five score intervals required")
    require(all(isinstance(p, list) and len(p) == 2 for p in scores), "two endpoints required")
    exact = [[exact_endpoint(x) for x in p] for p in scores]
    require(all(l <= h for l, h in exact), "reversed score interval")
    return exact


def margins(scores, pairs):
    exact = check_scores(scores)
    return [exact[i][0]-exact[j][1] for i, j in pairs]


class NativeEvaluation:
    """Direct ctypes binding independent of both previously used Python wrappers."""
    def __init__(self, library, runtime_dirs):
        self.handles = []
        if os.name == "nt":
            for directory in (library.parent, *runtime_dirs):
                self.handles.append(os.add_dll_directory(str(directory.resolve())))
        self.lib = ctypes.CDLL(str(library.resolve()))
        ptr = ctypes.POINTER(ctypes.c_double)
        self.lib.kernel_precision.argtypes = [ctypes.c_uint]
        self.lib.kernel_precision.restype = None
        self.lib.kernel_precision(256)
        self.lib.kernel_scores_real.argtypes = [ptr, ptr, ctypes.c_double, ctypes.c_double,
                                                ctypes.c_double, ctypes.c_double, ctypes.c_int, ptr]
        self.lib.kernel_scores_real.restype = ctypes.c_int
        self.lib.kernel_error.restype = ctypes.c_char_p
        self.lib.kernel_version.restype = ctypes.c_char_p
        self.version = self.lib.kernel_version().decode("ascii")
        self.calls = 0

    def scores(self, low, high, qlow, qhigh, merge):
        low_values = (ctypes.c_double*200)(*low)
        high_values = (ctypes.c_double*200)(*high)
        output = (ctypes.c_double*10)()
        ql, qh = float(qlow), float(qhigh)
        require(F.from_float(ql) == qlow and F.from_float(qh) == qhigh,
                "q endpoints must be exactly representable by the native API")
        status = self.lib.kernel_scores_real(low_values, high_values, ql, qh, ql, qh, merge, output)
        self.calls += 1
        require(status == 0, "native error: " + self.lib.kernel_error().decode("utf-8"))
        result = [[output[2*i], output[2*i+1]] for i in range(5)]
        check_scores(result)
        return result


def check_structure(document, args):
    expected_order = args.expected_order.split(",")
    require(set(expected_order) == {"A1", "A2", "A3", "A4", "A5"} and len(expected_order) == 5,
            "expected order must list every alternative once")
    index_order = [int(label[1:])-1 for label in expected_order]
    pairs = list(zip(index_order, index_order[1:]))
    root = (args.expected_q_min, args.expected_q_max)
    require(root[0] >= 1 and root[1] > root[0], "invalid independently declared q domain")
    require(document["status"] == "PASS", "certificate is not PASS")
    require(document["model"] == args.expected_model, "wrong model")
    expected_merge = 0 if args.expected_model == "H" else 1
    require(type(document["merge"]) is int and document["merge"] == expected_merge, "wrong native merge")
    require(document["q_domain_exact"] == [str(q) for q in root], "wrong root q domain")
    require(F(document["canonical_cube_radius_exact"]) == args.expected_radius > 0, "wrong canonical radius")
    require(document["order"] == expected_order, "wrong full order")
    require(document["adjacent_pairs"] == [[f"A{i+1}", f"A{j+1}"] for i, j in pairs], "wrong adjacent comparisons")
    require(document["all_200_coordinates_perturbed"] is True and document["contains_equal_radius_L2_ball"] is True,
            "wrong input-uncertainty scope")
    require(document["generation_bits"] == 192 and document["replay_bits"] == 256, "wrong retained precision metadata")
    require(document["unresolved_count"] == 0 and document["unresolved"] == [], "unresolved nodes present")
    require(set(document["source_hashes"]) == set(CPP_NAMES+HISTORICAL_PY_NAMES), "wrong source manifest")
    for name in CPP_NAMES:
        require(document["source_hashes"][name] == sha256(args.kernel_sources/name), "kernel source hash mismatch: " + name)
    require(re.fullmatch(r"[0-9a-f]{64}", document["native_library_sha256"]) is not None, "invalid historical library hash")
    if args.validation_mode == "retained-binary":
        require(document["native_library_sha256"] == sha256(args.library), "native library hash mismatch")

    # Recover exact canonical centres directly from the hashed C++ source.
    match = re.search(r"DATA\[200\]\s*=\s*\{([^}]+)\}", (args.kernel_sources/CPP_NAMES[1]).read_text(encoding="utf-8"))
    require(match is not None, "source DATA array absent")
    data = [int(v.strip()) for v in match.group(1).split(",")]
    require(len(data) == 200 and all(0 <= v <= 10 for v in data), "invalid source DATA array")
    center = [F(v**3, 1000) for v in data]
    expected_low = [outward(c-args.expected_radius, False) for c in center]
    expected_high = [outward(c+args.expected_radius, True) for c in center]
    low, high = document["input_intervals"]["lower"], document["input_intervals"]["upper"]
    require(len(low) == len(high) == 200 and low == expected_low and high == expected_high,
            "input cube differs from independently reconstructed radius and exact centre")
    require(all(0 <= exact_endpoint(l) <= exact_endpoint(h) <= 1 for l, h in zip(low, high)), "cube outside [0,1]")
    require(all(exact_endpoint(high[i])+exact_endpoint(high[i+1]) <= 1 for i in range(0, 200, 2)),
            "cube outside qROF product simplex")

    nodes, leaves = document["all_nodes"], document["leaves"]
    require(type(document["node_count"]) is int and type(document["leaf_count"]) is int, "integer counts required")
    require(leaves and len(nodes) == document["node_count"] == 2*len(leaves)-1
            and len(leaves) == document["leaf_count"], "wrong full binary-tree counts")
    by_id = {}
    for row in nodes:
        require(row["id"] not in by_id, "duplicate tree node")
        ql, qh = interval_for_id(row["id"], root)
        require((F(row["q_lower"]), F(row["q_upper"])) == (ql, qh), "node interval differs from root and dyadic path")
        calculated = margins(row["score_intervals"], pairs)
        require(row["adjacent_lowers_exact"] == [str(m) for m in calculated], "false stored adjacent margins")
        by_id[row["id"]] = row
    require("r" in by_id, "root node absent")
    leaf_ids, gen_margins, old_replay_margins = set(), [], []
    for row in leaves:
        require(row["id"] not in leaf_ids and by_id.get(row["id"]) == row, "duplicate or inconsistent leaf")
        leaf_ids.add(row["id"])
        generated = margins(row["score_intervals"], pairs)
        old_replay = margins(row["replay_256_score_intervals"], pairs)
        require(min(generated) > 0 and min(old_replay) > 0, "nonpositive retained adjacent bound")
        require(row["replay_256_adjacent_lowers_exact"] == [str(m) for m in old_replay], "false stored replay margins")
        for (gl, gh), (rl, rh) in zip(row["score_intervals"], row["replay_256_score_intervals"]):
            require(gl <= rl <= rh <= gh, "retained replay not contained in generation interval")
        gen_margins.extend(generated)
        old_replay_margins.extend(old_replay)
    for identifier in by_id:
        require(identifier == "r" or identifier[:-1] in by_id, "orphan tree node")
        children = [identifier+"0", identifier+"1"]
        if identifier in leaf_ids:
            require(not any(child in by_id for child in children), "retained leaf has children")
        else:
            require(all(child in by_id for child in children), "internal node missing child")
    require(F(leaves[0]["q_lower"]) == root[0] and F(leaves[-1]["q_upper"]) == root[1], "cover does not span root")
    require(all(F(a["q_upper"]) == F(b["q_lower"]) for a, b in zip(leaves, leaves[1:])), "gap, overlap, or unsorted cover")
    require(document["generation_uniform_margin_floor_12dp"] == floor12(min(gen_margins)), "false generation uniform summary")
    require(document["replay_uniform_margin_floor_12dp"] == floor12(min(old_replay_margins)), "false replay uniform summary")
    return leaves, low, high, pairs, expected_merge


def replay(document, args, native):
    leaves, low, high, pairs, merge = check_structure(document, args)
    minima = []
    for row in leaves:
        fresh = native.scores(low, high, F(row["q_lower"]), F(row["q_upper"]), merge)
        for key in ("score_intervals", "replay_256_score_intervals"):
            for (sl, sh), (fl, fh) in zip(row[key], fresh):
                require(sl <= fl <= fh <= sh, "fresh interval not contained in stored " + key + ": " + row["id"])
        lower = margins(fresh, pairs)
        require(min(lower) > 0, "fresh full-order margin not strictly positive: " + row["id"])
        minima.extend(lower)
    return {"retained_leaves_replayed": len(leaves), "fresh_uniform_margin_exact": str(min(minima)),
            "fresh_uniform_margin_floor_12dp": floor12(min(minima)), "bits": 256}


def malformed_tests(document, args, native):
    changes = [
        ("missing_leaf", lambda d: d["leaves"].pop()),
        ("wrong_root_domain", lambda d: d.__setitem__("q_domain_exact", ["4", "10"])),
        ("wrong_root_node", lambda d: next(r for r in d["all_nodes"] if r["id"] == "r").__setitem__("q_lower", "4")),
        ("wrong_model", lambda d: d.__setitem__("model", "A" if args.expected_model == "H" else "H")),
        ("wrong_merge", lambda d: d.__setitem__("merge", 1-d["merge"])),
        ("wrong_radius", lambda d: d.__setitem__("canonical_cube_radius_exact", str(args.expected_radius*2))),
        ("wrong_order", lambda d: d["order"].reverse()),
        ("missing_tree_child", lambda d: d["all_nodes"].pop()),
        ("shrunk_input_cube", lambda d: d["input_intervals"]["lower"].__setitem__(0, d["input_intervals"]["upper"][0])),
        ("false_stored_margin", lambda d: d["leaves"][0]["adjacent_lowers_exact"].__setitem__(0, "1")),
        ("wrong_source_hash", lambda d: d["source_hashes"].__setitem__(CPP_NAMES[0], "0"*64)),
        ("unexpected_unresolved", lambda d: d["unresolved"].append(d["leaves"][0])),
    ]
    if args.validation_mode == "retained-binary":
        changes.append(("wrong_library_hash", lambda d: d.__setitem__("native_library_sha256", "0"*64)))
    results = []
    for name, change in changes:
        altered = copy.deepcopy(document)
        change(altered)
        try:
            check_structure(altered, args)
        except (InvalidCertificate, KeyError, TypeError, IndexError) as exc:
            results.append({"test": name, "rejected": True, "phase": "structure", "reason": str(exc)})
        else:
            raise InvalidCertificate("malformed certificate accepted: " + name)

    # Forge *both* stored precisions, matching tree copies, exact margins and
    # summaries, so rejection requires actual independent native evaluation.
    forged = copy.deepcopy(document)
    row = forged["leaves"][0]
    pairs = [(int(i[1:])-1, int(j[1:])-1) for i, j in forged["adjacent_pairs"]]
    for key, margin_key in (("score_intervals", "adjacent_lowers_exact"),
                            ("replay_256_score_intervals", "replay_256_adjacent_lowers_exact")):
        lo, hi = row[key][0]
        row[key][0][0] = math.nextafter((lo+hi)/2, hi)
        row[margin_key] = [str(m) for m in margins(row[key], pairs)]
    forged["all_nodes"] = [copy.deepcopy(row) if n["id"] == row["id"] else n for n in forged["all_nodes"]]
    forged["generation_uniform_margin_floor_12dp"] = floor12(min(F(m) for r in forged["leaves"] for m in r["adjacent_lowers_exact"]))
    forged["replay_uniform_margin_floor_12dp"] = floor12(min(F(m) for r in forged["leaves"] for m in r["replay_256_adjacent_lowers_exact"]))
    check_structure(forged, args)
    try:
        replay(forged, args, native)
    except InvalidCertificate as exc:
        results.append({"test": "self_consistent_false_narrow_bounds", "rejected": True,
                        "phase": "fresh_native_evaluation", "structure_passed": True, "reason": str(exc)})
    else:
        raise InvalidCertificate("fresh replay accepted false narrowed bounds")
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--library", type=Path, required=True)
    parser.add_argument("--kernel-sources", type=Path, required=True)
    parser.add_argument("--runtime-dir", type=Path, action="append", default=[])
    parser.add_argument("--expected-model", choices=("H", "A"), required=True)
    parser.add_argument("--expected-radius", type=F, required=True)
    parser.add_argument("--expected-order", required=True)
    parser.add_argument("--expected-q-min", type=F, default=F(3))
    parser.add_argument("--expected-q-max", type=F, default=F(10))
    parser.add_argument("--validation-mode", choices=("retained-binary", "rebuilt-binary"), default="retained-binary")
    parser.add_argument("--malformed-tests", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    document = json.loads(args.certificate.read_text(encoding="utf-8"))
    # Reject malformed structure and provenance before loading a native library.
    check_structure(document, args)
    native = NativeEvaluation(args.library, args.runtime_dir)
    result = replay(document, args, native)
    result.update({"status": "PASS", "certificate_sha256": sha256(args.certificate),
                   "checker_sha256": sha256(Path(__file__)), "validation_mode": args.validation_mode,
                   "evaluated_library_sha256": sha256(args.library),
                   "historical_library_sha256": document["native_library_sha256"],
                   "kernel_source_hashes": {name: sha256(args.kernel_sources/name) for name in CPP_NAMES},
                   "expected_model": args.expected_model, "expected_radius_exact": str(args.expected_radius),
                   "expected_order": args.expected_order.split(","),
                   "expected_q_domain_exact": [str(args.expected_q_min), str(args.expected_q_max)],
                   "mpfr_version": native.version,
                   "coverage": "complete closed dyadic cover; independently reconstructed input cube",
                   "independence": "standalone Python checker, no generator/helper imports; shared C++ interval/model kernel",
                   "historical_python_hashes": "recorded in original certificate; not executed or required by this independent checker",
                   "limits": "conditional on source-model transcription, compiled kernel and MPFR/GMP; not an exact radius or proof-assistant result"})
    if args.malformed_tests:
        result["malformed_certificate_tests"] = malformed_tests(document, args, native)
    result["total_native_evaluations"] = native.calls
    result["elapsed_seconds"] = time.perf_counter()-started
    encoded = json.dumps(result, indent=2)+"\n"
    if args.output:
        require(args.output.resolve() != args.certificate.resolve(), "output must not overwrite retained certificate")
        args.output.write_text(encoded, encoding="utf-8")
    print(encoded)


if __name__ == "__main__":
    main()
