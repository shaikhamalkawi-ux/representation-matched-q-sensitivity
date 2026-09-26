from pathlib import Path
import importlib.util, math, json

ROOT=Path(__file__).resolve().parents[1]
CASE=ROOT/"case07_alkan_kahraman_2021"
spec=importlib.util.spec_from_file_location("ak",CASE/"reproduce_alkan2021.py")
ak=importlib.util.module_from_spec(spec)
spec.loader.exec_module(ak)

pairs=list(ak.ascale.values())+list(ak.wscale.values())
max_score_dev=0.0
max_comp_dev=0.0
max_wg_dev=0.0

# Exact-data verification across all declared scale pairs and q=2..50.
for x in pairs:
    base=ak.score(x,5)
    for q in range(2,51):
        tx=ak.transport(x,5,q)
        max_score_dev=max(max_score_dev,abs(ak.score(tx,q)-base))
        # complement commutes with transport
        c1=ak.transport((x[1],x[0]),5,q)
        c2=(tx[1],tx[0])
        max_comp_dev=max(max_comp_dev,abs(c1[0]-c2[0]),abs(c1[1]-c2[1]))

# q-ROFWG equivariance on every raw case-7 cell
for c in ak.criteria:
    terms=ak.raw[c]
    for i,a in enumerate(ak.alts):
        vals5=[ak.ascale[terms[d*7+i]] for d in range(3)]
        agg5=ak.qrofwg(vals5,ak.dmw,5)
        for q in range(2,51):
            valsq=[ak.transport(x,5,q) for x in vals5]
            lhs=ak.qrofwg(valsq,ak.dmw,q)
            rhs=ak.transport(agg5,5,q)
            max_wg_dev=max(max_wg_dev,abs(lhs[0]-rhs[0]),abs(lhs[1]-rhs[1]))

result={
    "max_q_score_transport_error":max_score_dev,
    "max_complement_commutation_error":max_comp_dev,
    "max_qrofwg_equivariance_error_case7_cells":max_wg_dev,
    "status":"PASS" if max(max_score_dev,max_comp_dev,max_wg_dev)<1e-12 else "FAIL"
}
(ROOT/"theory"/"proof_check_results.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
print(json.dumps(result,indent=2))
