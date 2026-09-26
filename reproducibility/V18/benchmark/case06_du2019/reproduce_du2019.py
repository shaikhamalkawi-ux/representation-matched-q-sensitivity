from qrung_iso_control import *
import csv, json
from pathlib import Path

OUT = Path("results")
OUT.mkdir(exist_ok=True)

weights = [0.2, 0.1, 0.3, 0.15, 0.15, 0.1]
X = {
    "A1": [(0.5,0.2),(0.8,0.3),(0.8,0.3),(0.7,0.3),(0.4,0.2),(0.4,0.8)],
    "A2": [(0.6,0.3),(0.5,0.8),(0.6,0.5),(0.6,0.5),(0.7,0.4),(0.5,0.6)],
    "A3": [(0.3,0.4),(0.8,0.5),(0.7,0.6),(0.6,0.4),(0.6,0.2),(0.4,0.7)],
    "A4": [(0.7,0.4),(0.5,0.6),(0.7,0.4),(0.5,0.5),(0.7,0.6),(0.6,0.5)],
    "A5": [(0.7,0.6),(0.6,0.4),(0.4,0.7),(0.4,0.3),(0.7,0.7),(0.5,0.4)],
}
published_q3 = {
    "A1": (0.7159, 0.2846, 0.3438),
    "A2": (0.6102, 0.4618, 0.1287),
    "A3": (0.6557, 0.4360, 0.1991),
    "A4": (0.6646, 0.4667, 0.1918),
    "A5": (0.6083, 0.5274, 0.0784),
}

all_values = [x for row in X.values() for x in row]
q_min = minimum_integer_rung(all_values)

# Baseline reproduction
baseline = {}
baseline_pass = True
for alt, row in X.items():
    agg = qrof_wpm(row, weights, q=3, power=2)
    sc = q_score(agg, 3)
    baseline[alt] = {"mu": agg[0], "nu": agg[1], "q_score": sc}
    exp = published_q3[alt]
    baseline_pass &= (round(agg[0],4), round(agg[1],4), round(sc,4)) == exp

rows = []
baseline_rank_score = None
baseline_rank_ci = None
max_ctrl_score_dev = 0.0

for q in range(q_min, 51):
    raw_score, ctrl_score, raw_ci, ctrl_ci = {}, {}, {}, {}
    for alt, row in X.items():
        raw_agg = qrof_wpm(row, weights, q=q, power=2)
        ctrl_row = transport_row(row, r=3, q=q)
        ctrl_agg = qrof_wpm(ctrl_row, weights, q=q, power=2)

        raw_score[alt] = q_score(raw_agg, q)
        ctrl_score[alt] = q_score(ctrl_agg, q)
        raw_ci[alt] = euclidean_closeness(raw_agg)
        ctrl_ci[alt] = euclidean_closeness(ctrl_agg)

    if baseline_rank_score is None and q == 3:
        baseline_rank_score = rank_desc(ctrl_score)
        baseline_rank_ci = rank_desc(ctrl_ci)
        baseline_ctrl_scores = dict(ctrl_score)

    rows.append({
        "q": q,
        "raw_qscore_rank": ">".join(rank_desc(raw_score)),
        "controlled_qscore_rank": ">".join(rank_desc(ctrl_score)),
        "raw_CI_rank": ">".join(rank_desc(raw_ci)),
        "controlled_CI_rank": ">".join(rank_desc(ctrl_ci)),
        **{f"raw_qscore_{a}": raw_score[a] for a in X},
        **{f"ctrl_qscore_{a}": ctrl_score[a] for a in X},
        **{f"raw_CI_{a}": raw_ci[a] for a in X},
        **{f"ctrl_CI_{a}": ctrl_ci[a] for a in X},
    })

# score invariance is evaluated against the q0=3 controlled scores
for row in rows:
    for alt in X:
        max_ctrl_score_dev = max(
            max_ctrl_score_dev,
            abs(row[f"ctrl_qscore_{alt}"] - baseline_ctrl_scores[alt])
        )

def transitions(key):
    out = []
    prev = None
    for row in rows:
        cur = row[key]
        if prev is not None and cur != prev:
            out.append({"q": row["q"], "from": prev, "to": cur})
        prev = cur
    return out

with (OUT / "du2019_sweep_q2_q50.csv").open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

verification = {
    "case": "Du 2019 q-ROFWPM, r=2",
    "source_doi": "10.1002/int.22167",
    "q_min": q_min,
    "q3_published_rounding_reproduction": baseline_pass,
    "q3_recomputed": baseline,
    "raw_qscore_transitions_q2_q50": transitions("raw_qscore_rank"),
    "controlled_qscore_transitions_q2_q50": transitions("controlled_qscore_rank"),
    "raw_CI_transitions_q2_q50": transitions("raw_CI_rank"),
    "controlled_CI_transitions_q2_q50": transitions("controlled_CI_rank"),
    "max_controlled_qscore_deviation_vs_q3": max_ctrl_score_dev,
    "interpretation": {
        "q_score_endpoint": "Class I / transport-invariant under the q-natural pipeline.",
        "euclidean_closeness_endpoint": "Raw-coordinate endpoint; controlled stability is empirical, not theorem-guaranteed."
    }
}
(OUT / "du2019_verification.json").write_text(
    json.dumps(verification, indent=2), encoding="utf-8"
)

print(json.dumps(verification, indent=2))
