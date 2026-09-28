"""Portable I/O and hash gates; no network, source acquisition or installation."""
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import numpy as np

if not __debug__:
    raise RuntimeError("Optimized Python (-O/-OO) is unsupported: integrity assertions must execute")

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def save(path, obj):
    with Path(path).open("x", encoding="utf-8") as f:
        json.dump(obj,f,indent=2,ensure_ascii=False,allow_nan=False); f.write("\n")


def digest(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,ensure_ascii=False,separators=(",",":"),allow_nan=False).encode()).hexdigest()


def verify_module():
    manifest=load(HERE/"MANIFEST.json")
    for name,item in manifest["files"].items():
        p=HERE/name
        assert p.resolve().is_relative_to(HERE.resolve()) and p.is_file() and not p.is_symlink(),name
        assert p.stat().st_size==item["bytes"] and sha(p)==item["sha256"],name
    actual={p.relative_to(HERE).as_posix() for p in HERE.rglob("*") if p.is_file() and p.name!="MANIFEST.json" and "__pycache__" not in p.parts}
    assert actual==set(manifest["files"]),"Unexpected module payload"
    return {"status":"PASS","payloads":len(actual),"manifest_sha256":sha(HERE/"MANIFEST.json")}


def cohort(key):
    obj=load(HERE/"data"/(key+".json")); rows=obj["rows"]
    canonical=[{"id":r["id"],"mu":[str(Fraction(x)) for x in r["mu"]],"nu":[str(Fraction(x)) for x in r["nu"]]} for r in rows]
    ids=[r["id"] for r in rows]; fu=[[Fraction(x) for x in r["mu"]] for r in canonical]; fv=[[Fraction(x) for x in r["nu"]] for r in canonical]
    assert digest(ids)==obj["ids_sha256"] and digest(canonical)==obj["canonical_fraction_rows_sha256"]
    assert len(ids)==len(set(ids)) and [len(rows),len(obj["criteria"])]==obj["shape"]
    assert all(0<=a<=1 and 0<=b<=1 and a+b<=1 for u,v in zip(fu,fv) for a,b in zip(u,v))
    return {"key":key,"inventory":{k:obj[k] for k in ("shape","criteria","source_id","canonical_fraction_rows_sha256","ids_sha256")},
            "ids":ids,"names":[r["name"] for r in rows],"fractions_u":fu,"fractions_v":fv,
            "u":np.array(fu,dtype=float),"v":np.array(fv,dtype=float)}


def comparison_report(report):
    omit={"started_utc","completed_utc","inventory","score_array_pin","exact_dominance_array_pin"}
    return {k:v for k,v in report.items() if k not in omit}


def new_output(path):
    path=Path(path).resolve()
    assert not path.is_relative_to(HERE.resolve()),"Outputs must be outside immutable release module"
    path.mkdir(parents=True,exist_ok=False)
    return path
