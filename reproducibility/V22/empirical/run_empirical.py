"""Replay the unchanged eleven-table comparison from any working directory."""
import argparse
from fractions import Fraction
import json
from pathlib import Path
import sys
import time
import traceback
sys.dont_write_bytecode=True
import numpy as np
import model_core as core
import portable_io as io


def historical_wa(key, cohort, array):
    # Distinct original S12 six/full-grid pair diagnostic, not R15 seven-node maximum.
    subset=(1,2,3,4,8,16) if key.startswith(("cms_","epa_")) else None
    indices=[int((q-1)*64) for q in subset] if subset else range(961)
    counts=[core.pair_categories(array[0],array[k])["strict_discordance"] for k in indices]
    result={"pair_scope":"six_nodes" if subset else "all961_nodes","max_strict_discordance":max(counts)}
    if key.startswith("pisa_"):
        winner=cohort["names"].index("Singapore"); u,v=cohort["fractions_u"],cohort["fractions_v"]
        assert all(a>0 and b>0 and a+b<1 for uu,vv in zip(u,v) for a,b in zip(uu,vv))
        for j in range(len(u)):
            if j==winner:continue
            assert all(a>=b for a,b in zip(u[winner],u[j])) and all(a<=b for a,b in zip(v[winner],v[j]))
            assert u[winner]!=u[j] or v[winner]!=v[j]
        result["standard_all_finite_q_wa_dominance"]="Singapore; exact rational strict-domain coordinate test"
    if key.startswith("epa_"):
        sat=[name for name,uu,vv in zip(cohort["names"],cohort["fractions_u"],cohort["fractions_v"])
             if any(a==1 for a in uu) and any(b==0 for b in vv)]
        result["standard_all_finite_q_wa_saturation_names"]=sat
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args(); integrity=io.verify_module();out=io.new_output(args.out); start=time.perf_counter()
    result={"status":"RUNNING","integrity":integrity,"no_network_or_installation":True,"cohorts":[]}
    try:
        core.synthetic_tests(out/"synthetic")
        expected_historical=io.load(io.HERE/"expected/historical_wa.json")
        keys=io.load(io.HERE/"data/INDEX.json")["cohorts"]
        for key in keys:
            c=io.cohort(key);dest=out/key;core.evaluate_cohort(c,dest)
            report=io.comparison_report(io.load(dest/"RESULTS_R1.json"))
            assert report==io.load(io.HERE/"expected"/(key+"_summary.json")),"Summary differs: "+key
            arrays=np.load(dest/"ALL_METHOD_SCORES_R1.npz",allow_pickle=False)
            ref=np.load(io.HERE/"expected"/(key+"_seven_nodes.npz"),allow_pickle=False)
            errors={}; checks=0
            for name in ref.files:
                actual=arrays[name][list(core.KEY_INDEX)] if name=="weights_dynamic" or name.endswith("__dynamic") else arrays[name]
                if np.issubdtype(ref[name].dtype,np.number):
                    error=float(np.max(np.abs(actual-ref[name])));assert error<=3e-14,(key,name,error);errors[name]=error
                else:assert np.array_equal(actual,ref[name]),(key,name)
                checks+=actual.size
            old=historical_wa(key,c,arrays["wa__dynamic"])
            assert old==expected_historical[key],"Historical WA diagnostic differs: "+key
            # The source-defined fixed/equal controls are algebraic broadcasts, not repeated trials.
            assert all(report["methods"][m]["dynamic"]["consecutive_top_set_changes"]==0 for m in core.METHODS)
            io.save(dest/"PUBLIC_REPLAY_CHECK.json",{"status":"PASS","seven_node_scalar_checks":int(checks),"max_errors":errors,"historical_wa":old})
            result["cohorts"].append({"key":key,"shape":c["inventory"]["shape"],"all_grid_summary":"PASS","seven_node_reference":"PASS","historical_wa":old})
        assert len(result["cohorts"])==11
        result.update(status="PASS",grid_nodes=961,cohort_analyses=11,dynamic_method_paths=33,
                      top_changes=0,full_grid_model_evaluation=True,independent_model_implementation=False,
                      limitations=["Purpose-selected related encodings, not prevalence or held-out validation.","Grid checks are not continuum certificates.","No accuracy, statistical, causal or policy inference.","Original-source acquisition/extraction is not rerun."])
        assert io.verify_module()==integrity
    except Exception:
        result.update(status="FAIL",traceback=traceback.format_exc());raise
    finally:
        result["elapsed_seconds"]=time.perf_counter()-start;io.save(out/"EMPIRICAL_REPLAY.json",result)
        print(json.dumps({"status":result["status"],"cohorts":len(result["cohorts"]),"elapsed_seconds":result["elapsed_seconds"]}),flush=True)


if __name__=="__main__":main()
