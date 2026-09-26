"""Standard-library exact replay of the three V17 point witnesses.

No optimizer, NumPy or SciPy is imported. Acceptance requires all stored exact
fields to match a fresh integer-interval evaluation, plus rational feasibility.
The optional --cross-rung table reuses each same rational raw-baseline vector
unchanged; it is supplemental pointwise evidence, not a radius comparison proof.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
DEFAULT_MODEL = HERE.parent / 'vendor'
EXPECTED_HASHES = {
    'model_inputs.json':'71351416bb73b22f3781bad962a023b6ea7362e8020429aab97b47d3f36d1797',
    'exact_models.py':'8d5941c32cd7cf7140c3273ac6805a370f73bcea627d321cd2ba3c2bfc62f427',
    'dyadic_interval.py':'38cae911d1feb49562dcfd3693bb56a8e1bb8219c8afa51f89d44aff5998852e'
}


def require(condition, text):
    if not condition:
        raise ValueError(text)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--model-dir',type=Path,default=DEFAULT_MODEL)
    parser.add_argument('--record','--input',dest='record',type=Path,
                        default=HERE/'FIXED_RUNG_WITNESSES.json')
    parser.add_argument('--output',type=Path,default=HERE/'EXACT_WITNESS_REPLAY.json',
                        help='Fresh output path; existing files are never overwritten.')
    parser.add_argument('--cross-rung',action='store_true',
                        help='Include supplemental exact same-vector comparisons across q=4,8,16.')
    args=parser.parse_args()
    args.output=args.output.resolve()
    if args.output.exists():
        raise FileExistsError(f'Refusing to overwrite existing results: {args.output}')
    if not args.output.parent.is_dir():
        raise ValueError(f'Output parent does not exist: {args.output.parent}')
    model=args.model_dir.resolve()
    hashes={name:hashlib.sha256((model/name).read_bytes()).hexdigest() for name in EXPECTED_HASHES}
    require(hashes==EXPECTED_HASHES,'immutable model hash mismatch')
    sys.path.insert(0,str(model))
    from dyadic_interval import IV,PREC
    from exact_models import INPUTS,seikh
    center=[F(n,INPUTS['seikh']['denominator']) for n in INPUTS['seikh']['numerators']]
    record=json.loads(args.record.read_text(encoding='utf-8'))
    require(record['immutable_input_sha256']==hashes,'record model hash mismatch')
    require(record['plan']['q_values']==[4,8,16],'unexpected rung plan')
    require(len(record['witnesses'])==3,'expected exactly three point witnesses')
    require([w['q_exact'] for w in record['witnesses']]==['4','8','16'],'witness rung identity')

    def check(w):
        raw=[F(x) for x in w['raw_baseline_normalized_inputs']]
        require(len(raw)==40,'wrong source dimension')
        require(all(0<x<1 for x in raw),'not a positive interior point')
        can=[x**4 for x in raw]
        pair_max=max(can[i]+can[i+1] for i in range(0,40,2))
        require(pair_max<=1,'inadmissible pair')
        require(str(pair_max)==w['maximum_canonical_pair_sum'],'stored pair maximum mismatch')
        radius=max(abs(x-c) for x,c in zip(raw,center))
        require(str(radius)==w['actual_linf_change'],'stored actual radius mismatch')
        require(radius<=F(w['radius_limit']),'radius exceeds declared limit')
        require(w['admissible'] is True,'admissibility flag mismatch')
        q=F(w['q_exact'])
        require(q in (4,8,16),'wrong rung')
        require(w['winner']=='Y1','wrong witness target')
        sc,we=seikh(IV.frac(1/q),[IV.frac(x) for x in can],tight=True)
        require(w['scores_exact']==[x.endpoints() for x in sc],'stored scores mismatch')
        require(w['weights_exact']==[x.endpoints() for x in we],'stored weights mismatch')
        gap=IV(sc[1].lo-sc[0].hi,sc[1].hi-sc[0].lo)
        margins=[IV(sc[0].lo-sc[j].hi,sc[0].hi-sc[j].lo) for j in (1,2,3)]
        require(w['margin_Y2_minus_rival']==gap.endpoints(),'stored reversal margin mismatch')
        require(w['rival_margins_over_other_alternatives']==[x.endpoints() for x in margins],
                'stored winner margins mismatch')
        require(gap.hi<0 and min(x.lo for x in margins)>0,'not a strict unique Y1 winner')
        require(w['rival_is_unique_winner'] is True and w['status']=='PASS_WITNESS',
                'stored acceptance flags mismatch')
        require(w['bits']==PREC,'wrong arithmetic precision')
        return {'q':str(q),'radius':str(radius),'unique_winner':'Y1',
                'margin_Y2_minus_Y1':gap.endpoints(),'status':'PASS'}

    checked=[check(w) for w in record['witnesses']]
    negative=[]
    first=record['witnesses'][0]
    mutations=[
        ('actual radius',lambda w:w.update(actual_linf_change='1/100')),
        ('rung',lambda w:w.update(q_exact='8')),
        ('score',lambda w:w['scores_exact'][0].__setitem__(0,'0')),
        ('weight',lambda w:w['weights_exact'][0].__setitem__(0,'0')),
        ('margin',lambda w:w['margin_Y2_minus_rival'].__setitem__(1,'0')),
        ('winner',lambda w:w.update(winner='Y2')),
        ('dimension',lambda w:w['raw_baseline_normalized_inputs'].pop()),
        ('source coordinate',lambda w:w['raw_baseline_normalized_inputs'].__setitem__(0,'4/5')),
        ('radius limit',lambda w:w.update(radius_limit='1/100')),
        ('pair maximum',lambda w:w.update(maximum_canonical_pair_sum='0')),
        ('precision',lambda w:w.update(bits=53)),
    ]
    for name,mutate in mutations:
        altered=deepcopy(first);mutate(altered)
        try:
            check(altered)
        except (ValueError,KeyError,IndexError) as e:
            negative.append({'mutation':name,'rejected':True,'reason':str(e)})
        else:
            raise ValueError('negative test accepted: '+name)
    cross=[]
    for w in record['witnesses'] if args.cross_rung else []:
        can=[F(x)**4 for x in w['raw_baseline_normalized_inputs']]
        row={'source_witness_rung':w['q_exact'],'source_radius':w['actual_linf_change'],'evaluations':[]}
        for q in (4,8,16):
            sc,_=seikh(IV.frac(F(1,q)),[IV.frac(x) for x in can],tight=True)
            winners=[i for i in range(4) if all(sc[i].lo>sc[j].hi for j in range(4) if j!=i)]
            require(len(winners)==1,'cross-rung unique winner unresolved')
            gap=IV(sc[1].lo-sc[0].hi,sc[1].hi-sc[0].lo)
            row['evaluations'].append({'q':q,'unique_winner':f'Y{winners[0]+1}',
                                       'margin_Y2_minus_Y1':gap.endpoints()})
        cross.append(row)
    output={'status':'PASS','arithmetic':'Python integer224 interval arithmetic; exact rational source feasibility and L-infinity distance',
            'immutable_input_sha256':hashes,'checked_witnesses':checked,
            'negative_tests':negative,'negative_test_count':len(negative),
            'cross_rung_enabled':args.cross_rung,'cross_rung_same_raw_input_checks':cross,
            'scope':'Point witnesses prove upper bounds only. Cross-rung checks of individual points do not establish ordered optimal radii.'}
    with args.output.open('x',encoding='utf-8') as stream:
        stream.write(json.dumps(output,indent=2)+'\n')
    print(json.dumps({'status':'PASS','witness_count':len(checked),'negative_test_count':len(negative),
                      'cross_rung_winners':[[v['unique_winner'] for v in x['evaluations']] for x in cross]}))


if __name__=='__main__':
    main()
