"""Portable frozen-cover replay using unchanged original QOL interval functions.

Replaces only workspace-bound input/output and dependency location bookkeeping.
It does not rerun source extraction or modify any mathematical evaluator.
"""
import argparse
import copy
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
if os.name=='nt' and not str(HERE).startswith('\\\\?\\'):
    HERE=Path('\\\\?\\'+str(HERE))
ROOT=HERE/'original/outputs/FSS_REAL_DATA_20260928'
DATA=ROOT/'qol2023/root_results_R1/EXACT_INPUTS_R1.json'
INDEPENDENT=ROOT/'protocol_review/INDEPENDENT_INPUTS_R1.json'
COVER=ROOT/'qol2023/continuum_decimal_R1.json'
PINS={DATA:'4d739810763b3093e62756e35832018c095c1a9bd43c1f1453040290c983c5ee',
      INDEPENDENT:'fed011e72692db5822087162301ea91f0876b9c66ce7a03e666384a958925d18',
      COVER:'06c0bdbf87f8771f42c2f1ef05f0b95361368ab56f58d1dc9bfec36294075b56'}

def require(ok,message):
    if not ok:raise ValueError(message)

def intervals_of(leaves):
    values=[(F(v['p_lo']),F(v['p_hi'])) for v in leaves]
    require(0<len(values)<=4096,'cover size')
    require(values==sorted(values) and values[0][0]==F(1,16) and values[-1][1]==1,'endpoints/order')
    require(all(a<b for a,b in values) and all(a[1]==b[0] for a,b in zip(values,values[1:])),'complete contiguous positive-width cover')
    return values

def main():
    if not __debug__:raise RuntimeError('Original QOL kernels use assertions; do not use -O/-OO')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--arithmetic',choices=['decimal','mpfr'],required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();require(not args.output.exists(),'fresh output only')
    start=time.perf_counter()
    manifest=json.loads((HERE/'THEORY_MANIFEST.json').read_text())['payloads']
    code=ROOT/('qol2023/certify_window_decimal_r1.py' if args.arithmetic=='decimal' else 'protocol_review/independent_qol_interval_mpfr_r1.py')
    require(hashlib.sha256(code.read_bytes()).hexdigest()==manifest[code.relative_to(HERE).as_posix()]['sha256'],'original code pin')
    for path,digest in PINS.items():require(hashlib.sha256(path.read_bytes()).hexdigest()==digest,'exact input pin')
    spec=importlib.util.spec_from_file_location('original_qol_'+args.arithmetic,code)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    cover=json.loads(COVER.read_text(encoding='utf-8'))
    source=json.loads((DATA if args.arithmetic=='decimal' else INDEPENDENT).read_text(encoding='utf-8'))
    groups={}
    for label in ('primary','secondary'):
        old=cover['groups'][label];keys=old['criteria']
        require(old['status']=='PASS' and old['unresolved_leaves']==0 and old['unresolved']==[],'frozen cover complete')
        require(keys==(['Q1a_1','Q1b_1'] if label=='primary' else ['Q1a_1','Q1a_2','Q1a_3','Q1a_4','Q1a_5','Q1a_6','Q1a_7','Q1b_1','Q1b_2','Q1b_3']),'fixed criteria')
        intervals=intervals_of(old['leaves']);require(len(intervals)==(32 if label=='primary' else 61),'frozen leaf count')
        checked=[]
        if args.arithmetic=='decimal':
            prepared=module.prepare(source['rows'],keys)
            for prior,(left,right) in zip(old['leaves'],intervals):
                value=module.evaluate(left,right,prepared)
                require(value is not None and F(value['gap_lower'])>0,'positive all-rival Decimal gap')
                require(value=={k:v for k,v in prior.items() if k!='depth'},'frozen Decimal parity')
                checked.append(value)
            minimum=min(F(v['gap_lower']) for v in checked)
        else:
            rows=[r for r in source['rows'] if all(c in r['criteria'] for c in keys)];require(len(rows)==83,'all cities retained')
            model=module.Model(rows,keys)
            for left,right in intervals:
                value=model.eval(left,right);require(F(value['minimum_gap_lower'])>0,'positive all-rival MPFR gap')
                checked.append({'p_lo':str(left),'p_hi':str(right),**value})
            minimum=min(F(v['minimum_gap_lower']) for v in checked)
        groups[label]={'criteria':keys,'winner':'3102 Zürich','rows':83,'rivals':82,'accepted_leaves':len(checked),
                       'complete_contiguous_cover':True,'minimum_certified_gap_lower':str(minimum),'leaves':checked}
        print(label,len(checked),'positive leaves',flush=True)
    rejections=[]
    for kind in ('omit','duplicate','reverse'):
        leaves=copy.deepcopy(cover['groups']['primary']['leaves'])
        if kind=='omit':leaves.pop(0)
        elif kind=='duplicate':leaves.insert(0,dict(leaves[0]))
        else:leaves[0]['p_lo'],leaves[0]['p_hi']=leaves[0]['p_hi'],leaves[0]['p_lo']
        try:intervals_of(leaves)
        except ValueError:rejections.append(kind)
        else:raise ValueError('malformed cover accepted')
    for path,digest in PINS.items():require(hashlib.sha256(path.read_bytes()).hexdigest()==digest,'unchanged exact inputs')
    result={'status':'PASS_QOL_'+args.arithmetic.upper()+'_FROZEN_COVER_REPLAY','arithmetic':'Outward Decimal80' if args.arithmetic=='decimal' else 'Directed MPFR256',
            'original_code_sha256':hashlib.sha256(code.read_bytes()).hexdigest(),
            'exact_inputs':{p.relative_to(HERE).as_posix():s for p,s in PINS.items()},
            'source_spreadsheet_sha256_provenance_only':'49986134cde7cddd1858ad1536558e66698ae86d28252f4d9b911ee63d79c03c',
            'groups':groups,'q_closed_interval':['1','16'],'malformed_cover_rejections':rejections,
            'adapter':'Original evaluator functions unchanged; exact cover validation; portable input/output/dependency bookkeeping',
            'source_extraction_rerun':False,'dependency_binaries_bundled_or_identity_verified':False,
            'scope':'Post-screen fixed analyst-encoded table winner over q in[1,16]. No source uncertainty, all-q theorem, official city ranking or prevalence.',
            'elapsed_seconds':time.perf_counter()-start}
    with args.output.open('x',encoding='utf-8') as f:json.dump(result,f,ensure_ascii=False,indent=2);f.write('\n')
    print(result['status'],flush=True)

if __name__=='__main__':main()
