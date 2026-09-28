"""Independent original Decimal90 arithmetic at the seven specified nodes.

Only I/O was ported. Does not import the NumPy primary arithmetic evaluator.
Ordinary high precision is not outward rounding or a continuum certificate.
"""
import argparse
from decimal import Decimal as D, localcontext
import json
from pathlib import Path
import sys
import time
import traceback
sys.dont_write_bytecode=True
import numpy as np
import decimal_core as dc
import portable_io as io

NODES=(1,2,3,4,8,12,16)
METHODS=("wa","wg","cosine_topsis")


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args();integrity=io.verify_module();out=io.new_output(args.out);started=time.perf_counter()
    result={"status":"RUNNING","integrity":integrity,"precision":90,"outward":False,"cohorts":[]}
    try:
        with localcontext() as context:
            context.prec=90;io.save(out/"SYNTHETIC_GATES.json",dc.synthetic_gates());pi=dc.decimal_pi()
            for key in io.load(io.HERE/"data/INDEX.json")["cohorts"]:
                c=io.cohort(key);raw=io.load(io.HERE/"data"/(key+".json"))["rows"]
                rows=[{k:r[k] for k in ("id","mu","nu")} for r in raw]
                lu,lv,logs,cells=dc.caches_for(rows,pi)
                fixed=dc.entropy_weights(lu,lv,D(1));equal=[D(1)/len(fixed)]*len(fixed)
                fs=dc.scores(logs,cells,fixed);es=dc.scores(logs,cells,equal)
                reference=np.load(io.HERE/"expected"/(key+"_seven_nodes.npz"),allow_pickle=False)
                original_decimal=np.load(io.HERE/"expected"/(key+"_decimal_reference.npz"),allow_pickle=False)
                original_decimal_scores=original_decimal["scores"]
                maximum=D(0);old_maximum=0.;count=0
                # Explicitly retain all 60 fixed/equal scalar-weight checks across eleven cohorts.
                for label,values in (("weights_frozen_q1",fixed),("weights_equal",equal)):
                    refweights=reference[label]
                    for j,value in enumerate(values):
                        error=abs(value-D(str(refweights[j])));assert error<=D("3e-14")
                        maximum=max(maximum,error);count+=1
                for node,q in enumerate(NODES):
                    w=dc.entropy_weights(lu,lv,D(q));dyn=dc.scores(logs,cells,w)
                    for j,value in enumerate(w):
                        error=abs(value-D(str(reference["weights_dynamic"][node,j])));assert error<=D("3e-14");maximum=max(maximum,error);count+=1
                    for arm,scores in enumerate((dyn,fs,es)):
                        arm_name=("dynamic","frozen_q1","equal")[arm]
                        for method,vals in enumerate(scores):
                            target=reference[METHODS[method]+"__"+arm_name]
                            if arm==0:target=target[node]
                            for i,value in enumerate(vals):
                                error=abs(value-D(str(target[i])));assert error<=D("3e-14"),(key,node,method,i)
                                maximum=max(maximum,error);count+=1
                                old_error=abs(float(value)-float(original_decimal_scores[method,arm,node,i]));assert old_error<=2e-15
                                old_maximum=max(old_maximum,old_error)
                            floatvals=np.array([float(v) for v in vals]);tol=1e-12
                            oldr=1+len(target)-np.searchsorted(np.sort(target),target+tol,side="right")
                            newr=1+len(floatvals)-np.searchsorted(np.sort(floatvals),floatvals+tol,side="right")
                            assert np.array_equal(oldr,newr),(key,node,method,arm,"rank")
                            assert np.array_equal(floatvals>=floatvals.max()-tol,target>=target.max()-tol)
                item={"cohort":key,"status":"PASS","scalar_comparisons":count,"maximum_error_against_primary":str(maximum),
                      "maximum_error_against_original_Decimal_reference":old_maximum,"rank_top_agreement":True}
                result["cohorts"].append(item);print(json.dumps(item),flush=True)
        result.update(status="PASS",independent_arithmetic=True,nodes=list(NODES),source_extraction_rerun=False,
                      full_grid_model_evaluation=False,continuum_certificate=False,
                      scalar_comparisons=sum(v["scalar_comparisons"] for v in result["cohorts"]))
        assert io.verify_module()==integrity
    except Exception:
        result.update(status="FAIL",traceback=traceback.format_exc());raise
    finally:
        result["elapsed_seconds"]=time.perf_counter()-started;io.save(out/"DECIMAL_REPLAY.json",result)


if __name__=="__main__":main()
