from pathlib import Path
import json

ROOT=Path(__file__).resolve().parent
def J(rel):
    return json.loads((ROOT/rel).read_text(encoding="utf-8"))

c1=J("case01_pinar_boran_thesis/results/verification.json")
c2=J("case02_khan_cosine_topsis/results/verification.json")
c3=J("case03_liu_wang_qrofwa/results/verification.json")
c4=J("case04_seikh_mandal_archimedean/results/verification.json")
c5=J("case05_du2021_einstein/results/verification.json")
c6=J("case06_du2019/results/du2019_verification.json")
c7=J("case07_alkan_kahraman_2021/results/verification.json")
proof=J("theory/proof_check_results.json")

errors=[]
def require(cond,msg):
    if not cond: errors.append(msg)
def rungs(xs): return [x["q"] for x in xs]

# all case directories
for i in range(1,8):
    require(any(ROOT.glob(f"case{i:02d}_*")), f"case {i} directory missing")

# Case 1
require(all(c1["baseline_checks"].values()),"C1 baseline")
require(c1["figure_7_4a_all_90_printed_values_reproduced"],"C1 Figure 7.4a")
require(rungs(c1["raw_rank_transitions_q2_q50_high_precision"])==[6,10],"C1 raw rungs")
require(c1["controlled_rank_transitions_q2_q50_high_precision"]==[],"C1 controlled")

# Case 2
require(c2["table4_exact_to_printed_6dp"],"C2 Table4")
require(c2["table5_all_q1_q10_within_one_4dp_unit"],"C2 Table5")
require(rungs(c2["raw_rank_transitions_q1_q50"])==[3,4],"C2 raw rungs")
require(c2["controlled_rank_transitions_q1_q50"]==[],"C2 controlled")

# Case 3
require(c3["q3_table4_score_reproduction"],"C3 q3")
require(c3["table5_all_published_rows_reproduced_to_4dp"],"C3 Table5")
require(rungs(c3["raw_rank_transitions_q2_q50"])==[8,9],"C3 raw rungs")
require(c3["controlled_rank_transitions_q2_q50"]==[],"C3 controlled")
require(c3["qrofwg_audit"]["status"]=="HOLD_LOCALIZED_X1_MISMATCH","C3 WG hold")

# Case 4
require(c4["publication_matched_q4_using_printed_3dp_weights"]["pass"],"C4 baseline")
require(rungs(c4["raw_rank_transitions_q4_q50"])==[10],"C4 raw rungs")
require(c4["controlled_rank_transitions_q4_q50"]==[],"C4 controlled")

# Case 5
require(c5["q3_reproduction"]["pass"],"C5 q3")
require(c5["table2_all_published_q_values_reproduced_to_4dp"],"C5 Table2")
require(rungs(c5["raw_rank_transitions_q2_q50"])==[3,9,11],"C5 raw rungs")
require(c5["controlled_rank_transitions_q2_q50"]==[],"C5 controlled")

# Case 6
require(c6["q3_published_rounding_reproduction"],"C6 q3")
require(rungs(c6["raw_qscore_transitions_q2_q50"])==[7],"C6 raw q-score")
require(c6["controlled_qscore_transitions_q2_q50"]==[],"C6 controlled q-score")

# Case 7
require(all(c7["baseline_checks"].values()),"C7 baseline")
require(all(v for d in c7["published_table19_reproduction"].values() for v in d.values()),"C7 Table19")
require(rungs(c7["raw_transitions_q2_q50"]["method1"])==[3,8,15],"C7 M1 raw")
require(rungs(c7["raw_transitions_q2_q50"]["method2"])==[3,11],"C7 M2 raw")
require(c7["controlled_transitions_q2_q50"]["method1"]==[],"C7 M1 controlled")
require(c7["controlled_transitions_q2_q50"]["method2"]==[],"C7 M2 controlled")
require(c7["table14_audit"]["formula_weights_reproduce_table17_18"],"C7 Eq34 weights")
require(not c7["table14_audit"]["printed_table14_weights_reproduce_table17_18"],"C7 Table14 audit")

require(proof["status"]=="PASS","theory proof checks")

out={
    "release":"QRUNG_ISO_CONTROL_BENCHMARK_v0_4",
    "status":"PASS" if not errors else "FAIL",
    "errors":errors,
    "executable_cases":7,
    "raw_transition_cases":7,
    "controlled_transition_cases":0,
    "theorem_level_primary_cases":[2,3,5,6,7],
    "residual_cases":[1,4]
}
(ROOT/"RELEASE_VERIFICATION.json").write_text(json.dumps(out,indent=2),encoding="utf-8")
if errors:
    raise SystemExit("RELEASE VERIFICATION FAIL: "+"; ".join(errors))
print("RELEASE VERIFICATION: PASS")
