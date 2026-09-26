"""Zhang--Wang (2024), DOI 10.1016/j.eswa.2023.122574, source audit.

Run with Python + mpmath. All decimal source entries are exact decimal inputs.
Reads no private data. Writes zhang_*.json/csv outputs, beside this script by
default or under --output-dir. A publisher PDF is not needed for replay.
Source: Tables 8--10 p.14, Eq.(7) p.7, Definition 5.1 p.8,
Eq.(13) p.10, full algorithm pp.10--11, worked outputs p.15.
The source's weights are POSITION weights after ordering, not expert labels.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from fractions import Fraction
from pathlib import Path

import mpmath as mp

mp.mp.dps = 80
HERE = Path(__file__).resolve().parent
RECORDED_SOURCE_PDF_SHA256 = "e0905b42d032de9beca882b62a50c784f5429045b8b0a711ca2a19894b88c1f3"

# Each line is one alternative (Z1 through Z4); each row has P1 through P5.
# Three blocks are experts e1, e2, e3 (Tables 8, 9, 10, printed p.14).
SOURCE = [
    [".6 .3 .7 .4 .5 .4 .8 .4 .6 .6", ".8 .2 .7 .1 .7 .3 .9 .4 .7 .2",
     ".5 .4 .6 .5 .8 .6 .7 .5 .6 .3", ".6 .5 .8 .4 .7 .4 .8 .3 .8 .4"],
    [".4 .6 .8 .5 .7 .6 .6 .5 .8 .4", ".8 .1 .9 .3 .8 .3 .7 .2 .7 .1",
     ".7 .2 .7 .4 .8 .5 .6 .4 .6 .2", ".6 .3 .7 .3 .7 .2 .8 .2 .7 .2"],
    [".5 .6 .7 .5 .8 .6 .6 .3 .7 .7", ".8 .3 .9 .4 .8 .2 .7 .1 .8 .3",
     ".6 .5 .7 .7 .6 .6 .7 .4 .7 .6", ".5 .6 .8 .5 .7 .3 .6 .1 .5 .1"],
]
PRINTED_H = [
    ".5294 .4547 .7469 .4729 .7113 .5206 .7057 .4025 .7240 .5032",
    ".8000 .1747 .8575 .2195 .7711 .2551 .8138 .2219 .7469 .1978",
    ".6269 .3205 .6792 .5145 .7512 .5578 .6702 .4229 .6459 .3577",
    ".5701 .4691 .7712 .3824 .7000 .2821 .7512 .1569 .7082 .2071",
]
PRINTED_FINAL = ".7059 .4609 .8063 .2159 .6809 .4270 .7206 .2591"
PRINTED_SCORES = ".3723 .2780 .3842 .3491"
EXPERT_WEIGHTS = list(map(mp.mpf, [".40", ".25", ".35"]))
CRITERION_WEIGHTS = list(map(mp.mpf, [".20", ".30", ".15", ".20", ".15"]))


def pairs(s):
    a = list(map(mp.mpf, s.split()))
    return list(zip(a[::2], a[1::2]))


X = [[pairs(row) for row in expert] for expert in SOURCE]
BASE_Q = mp.mpf(3)
Z0 = [[[(mu**BASE_Q, nu**BASE_Q) for mu, nu in row] for row in e] for e in X]


def source_score(a, q):
    """Direct Eq.(7), including explicit distance and hesitancy factors."""
    mu, nu = a
    pq = 1 - mu**q - nu**q
    assert pq >= -mp.mpf("1e-70")
    distance = mp.sqrt(((1-mu**q)**2 + nu**(2*q)) / 2)
    return distance / mp.sqrt(1+pq)


def canonical_score(z):
    u, v = z
    assert 0 <= u and 0 <= v and u+v <= 1+mp.mpf("1e-70")
    return mp.sqrt(((1-u)**2+v*v)/(2*(2-u-v)))


def source_key(a, q):
    return source_score(a, q), 1-a[0]**q-a[1]**q


def canonical_key(z):
    return canonical_score(z), 1-z[0]-z[1]


def source_owa(a, w, q):
    """Direct Eq.(13); source order has smaller inverse score first."""
    order = sorted(range(len(a)), key=lambda i: source_key(a[i], q))
    b = [a[i] for i in order]
    mu = (1-mp.fprod((1-p[0]**q)**weight for p, weight in zip(b, w) if weight))**(1/q)
    nu = mp.fprod(p[1]**weight for p, weight in zip(b, w) if weight)
    return (mu, nu), order


def canonical_owa(z, w):
    order = sorted(range(len(z)), key=lambda i: canonical_key(z[i]))
    b = [z[i] for i in order]
    return (1-mp.fprod((1-p[0])**weight for p, weight in zip(b, w) if weight),
            mp.fprod(p[1]**weight for p, weight in zip(b, w) if weight)), order


def source_pipeline(x, q):
    h, expert_orders = [], []
    for i in range(4):
        row, orders = [], []
        for j in range(5):
            a, order = source_owa([x[k][i][j] for k in range(3)], EXPERT_WEIGHTS, q)
            row.append(a)
            orders.append(order)
        h.append(row)
        expert_orders.append(orders)
    final_orders = [source_owa(row, CRITERION_WEIGHTS, q) for row in h]
    final = [a for a, order in final_orders]
    scores = [source_score(a, q) for a in final]
    ranking = sorted(range(4), key=lambda i: source_key(final[i], q))
    return dict(h=h, final=final, scores=scores, ranking=ranking,
                expert_orders=expert_orders, criterion_orders=[o for a, o in final_orders])


def canonical_pipeline(z):
    h, expert_orders = [], []
    for i in range(4):
        result = [canonical_owa([z[k][i][j] for k in range(3)], EXPERT_WEIGHTS) for j in range(5)]
        h.append([a for a, o in result])
        expert_orders.append([o for a, o in result])
    final_orders = [canonical_owa(row, CRITERION_WEIGHTS) for row in h]
    final = [a for a, o in final_orders]
    return dict(h=h, final=final, scores=[canonical_score(a) for a in final],
                ranking=sorted(range(4), key=lambda i: canonical_key(final[i])),
                expert_orders=expert_orders, criterion_orders=[o for a, o in final_orders])


def n(x):
    return mp.nstr(x, 35)


def encode(obj):
    if isinstance(obj, mp.mpf):
        return n(obj)
    if isinstance(obj, dict):
        return {k: encode(v) for k, v in obj.items()}
    if isinstance(obj, (tuple, list)):
        return [encode(v) for v in obj]
    return obj


def ranking_text(order, prefix="Z"):
    return ">".join(f"{prefix}{i+1}" for i in order)


def rounded_agrees(a, printed):
    # Half-unit acceptance around the source's four-decimal displays.
    return abs(a-printed) <= mp.mpf(".00005")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdf", type=Path,
                        help="Optional local source PDF; verify its SHA-256 against the recorded source.")
    parser.add_argument("--output-dir", type=Path, default=HERE,
                        help="Output directory; default is the script directory.")
    args = parser.parse_args()
    source_verified = False
    if args.source_pdf is not None:
        actual_hash = hashlib.sha256(args.source_pdf.read_bytes()).hexdigest()
        if actual_hash != RECORDED_SOURCE_PDF_SHA256:
            raise SystemExit("Source PDF SHA-256 does not match the recorded source; audit aborted.")
        source_verified = True
    args.output_dir.mkdir(parents=True, exist_ok=True)
    baseline = source_pipeline(X, BASE_Q)
    canonical = canonical_pipeline(Z0)
    assert max(abs(a-b) for a, b in zip(baseline["scores"], canonical["scores"])) < mp.mpf("1e-70")
    assert baseline["ranking"] == [1, 3, 0, 2]
    source_compare = []
    for i, row in enumerate(map(pairs, PRINTED_H)):
        for j, a in enumerate(row):
            for c in range(2):
                actual = baseline["h"][i][j][c]
                source_compare.append(dict(stage="H", alternative=i+1, criterion=j+1,
                                           component="mu" if c==0 else "nu", printed=n(a[c]),
                                           calculated=n(actual), rounded_agrees=rounded_agrees(actual, a[c])))
    for i, a in enumerate(pairs(PRINTED_FINAL)):
        for c in range(2):
            actual = baseline["final"][i][c]
            source_compare.append(dict(stage="final", alternative=i+1, criterion="",
                                       component="mu" if c==0 else "nu", printed=n(a[c]),
                                       calculated=n(actual), rounded_agrees=rounded_agrees(actual, a[c])))
    for i, p in enumerate(map(mp.mpf, PRINTED_SCORES.split())):
        actual = baseline["scores"][i]
        source_compare.append(dict(stage="score", alternative=i+1, criterion="", component="rho",
                                   printed=n(p), calculated=n(actual), rounded_agrees=rounded_agrees(actual, p)))
    sweep = []
    for q_int in range(1, 51):
        q = mp.mpf(q_int)
        matched_x = [[[(mu**(BASE_Q/q), nu**(BASE_Q/q)) for mu, nu in row] for row in e] for e in X]
        matched = source_pipeline(matched_x, q)
        error = max(abs(a-b) for a,b in zip(matched["scores"], baseline["scores"]))
        assert error < mp.mpf("1e-70")
        assert matched["ranking"] == baseline["ranking"]
        assert matched["criterion_orders"] == baseline["criterion_orders"]
        # Identical inputs can exchange numerical order; their OWA contribution is identical.
        raw_admissible = all(mu**q+nu**q <= 1 for e in X for row in e for mu,nu in row)
        raw = source_pipeline(X, q) if raw_admissible else None
        sweep.append(dict(q=q_int, matched_score_error=n(error), matched_ranking=ranking_text(matched["ranking"]),
                          raw_admissible=raw_admissible,
                          raw_ranking=ranking_text(raw["ranking"]) if raw else "inadmissible",
                          raw_scores=";".join(n(x) for x in raw["scores"]) if raw else "",
                          matched_scores=";".join(n(x) for x in matched["scores"])))
    # Every input in source Table 6, p.13; no example selected on result.
    table6_mu = list(map(mp.mpf, ".12 .25 .36 .45 .54 .60 .63 .66 .68 .70 .72 .74 .76 .78 .79".split()))
    table6_printed = list(map(mp.mpf, ".4996 .4961 .4888 .4791 .4669 .4592 .4564 .4551 .4555 .4572 .4607 .4666 .4755 .4884 .4966".split()))
    table6_scores = [source_score((x,x), BASE_Q) for x in table6_mu]
    table6_ranks = sorted(range(15), key=lambda i: source_key((table6_mu[i],table6_mu[i]), BASE_Q))
    printed_rank_positions = [15,13,12,10,8,5,3,1,2,4,6,7,9,11,14]
    assert [table6_ranks.index(i)+1 for i in range(15)] == printed_rank_positions
    assert all(rounded_agrees(a,b) for a,b in zip(table6_scores,table6_printed))
    table6 = []
    table6_fractions = [Fraction(s) for s in ".12 .25 .36 .45 .54 .60 .63 .66 .68 .70 .72 .74 .76 .78 .79".split()]
    for k in range(3,51):
        raw_scores = [source_score((x,x), mp.mpf(k)) for x in table6_mu]
        raw_rank = sorted(range(15), key=lambda i: source_key((table6_mu[i],table6_mu[i]), mp.mpf(k)))
        transported = [(x**(BASE_Q/mp.mpf(k)),)*2 for x in table6_mu]
        matched_scores = [source_score(a, mp.mpf(k)) for a in transported]
        matched_rank = sorted(range(15), key=lambda i: source_key(transported[i], mp.mpf(k)))
        error = max(abs(a-b) for a,b in zip(table6_scores,matched_scores))
        assert error < mp.mpf("1e-70") and matched_rank == table6_ranks
        # Ranking by squared score avoids all irrational operations. The exact
        # rational margin below certifies every raw winner against its runner-up.
        powers = [x**k for x in table6_fractions]
        square_scores = [((1-t)**2+t*t)/(4*(1-t)) for t in powers]
        exact_order = sorted(range(15), key=lambda i: (square_scores[i],1-2*powers[i]))
        assert exact_order == raw_rank
        margin = square_scores[exact_order[1]]-square_scores[exact_order[0]]
        assert margin > 0
        table6.append(dict(q=k,raw_ranking=ranking_text(raw_rank,"a"),
                           matched_ranking=ranking_text(matched_rank,"a"),
                           raw_winner=f"a{raw_rank[0]+1}", raw_runner_up=f"a{raw_rank[1]+1}",
                           raw_squared_score_margin_exact=str(margin),
                           raw_squared_score_margin_decimal=n(mp.mpf(margin.numerator)/margin.denominator),
                           matched_score_error=n(error), raw_scores=";".join(n(x) for x in raw_scores)))
    # Exact tie audit on Delta: (0,0) and (1/2,1/2) both have rho=1/2.
    assert canonical_score((mp.mpf(0),mp.mpf(0))) == mp.mpf(".5")
    assert canonical_score((mp.mpf(".5"),mp.mpf(".5"))) == mp.mpf(".5")
    tie_control = {"canonical_points": [[0,0],[".5",".5"]], "rho": ".5",
                   "powered_hesitancies": [1,0], "preferred_index": 1,
                   "note": "Definition 5.1 prefers smaller powered hesitancy; preserved at every q."}
    result = dict(source_doi="10.1016/j.eswa.2023.122574", source_pdf_sha256=RECORDED_SOURCE_PDF_SHA256,
                  source_pdf_freshly_verified=source_verified,
                  source_pdf_provenance_status="fresh local hash verification" if source_verified else "recorded hash only; source PDF not freshly verified",
                  arithmetic="mpmath, 80 decimal digits; numerical regression, not directed-rounding certificate",
                  input_records=60, input_coordinates=120, source_integer_minimum_q=2,
                  baseline_q=3, baseline=baseline, canonical=canonical,
                  source_comparison_count=len(source_compare),
                  source_comparison_agreements=sum(x["rounded_agrees"] for x in source_compare),
                  source_mismatches=[x for x in source_compare if not x["rounded_agrees"]],
                  all_matched_rungs=list(range(1,51)),
                  maximum_matched_score_error=max(mp.mpf(row["matched_score_error"]) for row in sweep),
                  raw_rankings={row["q"]:row["raw_ranking"] for row in sweep},
                  table6=dict(input_count=15, scores=table6_scores, printed_rank_positions=printed_rank_positions,
                              printed_rank_agreements=15, printed_score_agreements=sum(rounded_agrees(a,b) for a,b in zip(table6_scores,table6_printed)),
                              score_disagreements=[dict(index=i+1, actual=a, printed=b) for i,(a,b) in enumerate(zip(table6_scores,table6_printed)) if not rounded_agrees(a,b)],
                              ranking=ranking_text(table6_ranks,"a"), sweep=table6), tie_control=tie_control)
    (args.output_dir/"zhang_pipeline_results.json").write_text(json.dumps(encode(result), indent=2)+"\n", encoding="utf-8")
    for name, data in [("zhang_pipeline_sweep.csv",sweep),("zhang_source_comparison.csv",source_compare),
                       ("zhang_table6_sweep.csv",table6)]:
        with (args.output_dir/name).open("w", newline="", encoding="utf-8") as f:
            writer=csv.DictWriter(f, fieldnames=list(data[0]))
            writer.writeheader()
            writer.writerows(data)
    print(json.dumps(encode({"baseline_scores":baseline["scores"],"baseline_final":baseline["final"],
                            "source_agreements":result["source_comparison_agreements"],
                            "source_mismatches":result["source_mismatches"],
                            "max_matched_error":result["maximum_matched_score_error"],
                            "table6_printed_scores_reproduced":result["table6"]["printed_score_agreements"],
                            "table6_winner_transitions":[row for i,row in enumerate(table6) if i==0 or row["raw_winner"] != table6[i-1]["raw_winner"]],
                            "raw_rankings":sorted(set(x["raw_ranking"] for x in sweep))}),indent=2))


if __name__ == "__main__":
    main()
