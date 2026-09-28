"""Independent rounding fragility audit of the frozen seven-row compact source.

No import of either MPFR compact model evaluator. Uses the separate Decimal100
outward arithmetic primitive, source-form score identity, exact raw boxes, exact
half-even Fraction rounding, and the original seven unchanged rational q nodes.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
SOURCE=ROOT/'outputs/FSS_EIGHT_HOUR_RESEARCH_20260926/external_case/compact_multiwinner/COMPACT_SOURCE_AND_NODES.json'
SOURCE_SHA='2f68987986dc54e3981732e59b8ad08d636531ef696912344bb04d9b571cfbe0'
KERNEL=HERE.parent/'pinar_independent/replay_decimal.py'
KERNEL_SHA='82fc5a524072a0da7dac34aafa53e3b96db39df727d42c3bc8158b3adfd36235'
if hashlib.sha256(KERNEL.read_bytes()).hexdigest()!=KERNEL_SHA:raise ValueError('Decimal primitive pin')
spec=importlib.util.spec_from_file_location('compact_decimal_primitive',KERNEL)
kernel=importlib.util.module_from_spec(spec);spec.loader.exec_module(kernel)
Box=kernel.Box
require=kernel.require
LABELS=['C1','C2','C3','D1','D2','D3','D4']


def point(d):return Box.direct(d,d)


def entropy_atom(y):
    require(y.hi<0,'entropy argument must be negative')
    endpoints=[point(d)*point(d).exp() for d in [y.lo,y.hi]]
    lo=min(x.lo for x in endpoints)
    if y.lo<=-1<=y.hi:lo=min(lo,(-Box(-1).exp()).lo)
    return Box.direct(lo,max(x.hi for x in endpoints))


def make_box(source,h):
    logs=[];complements=[];nonmembers=[]
    checked=0
    for row in source:
        rowlogs=[];rowcomp=[];rownon=[]
        for mu,nu in row:
            ml,mh=mu-h,mu+h;nl,nh=nu-h,nu+h
            require(F(1,16)<ml<=mh<F(7,8) and F(1,16)<nl<=nh<F(7,8),'source raw compact domain')
            require(mh**4+nh**4<1,'entire raw source pair box inadmissible')
            logmu=Box(ml,mh).log();lognu=Box(nl,nh).log()
            rowlogs.append((logmu,lognu))
            # Monotone fourth power at exact rational raw endpoints.
            rowcomp.append(Box(1-mh**4,1-ml**4).log())
            rownon.append(4*lognu)
            checked+=2
        logs.append(rowlogs);complements.append(rowcomp);nonmembers.append(rownon)
    require(checked==28,'coordinate count')
    return logs,complements,nonmembers


def evaluate(precomputed,q):
    logs,complements,nonmembers=precomputed
    p=F(4)/q
    require(q>=4,'q below r4')
    factors=[1+sum(entropy_atom(x*p) for row in logs for x in row[j])/7 for j in range(2)]
    require(all(e.lo>0 for e in factors),'entropy factor positivity')
    # The exact range of e1/(e1+e2) over the positive factor rectangle.
    e1,e2=factors
    low=point(e1.lo)/(point(e1.lo)+point(e2.hi))
    high=point(e1.hi)/(point(e1.hi)+point(e2.lo))
    w1=Box.direct(low.lo,high.hi)
    weights=[w1,1-w1]
    require(all(w.lo>0 and w.hi<1 for w in weights),'weight domain')
    scores=[1-sum(w*x for w,x in zip(weights,co)).exp()-sum(w*x for w,x in zip(weights,no)).exp()
            for co,no in zip(complements,nonmembers)]
    return scores,factors,weights


def check_nodes(source,nodes,h):
    precomputed=make_box(source,h)
    results=[]
    for node in nodes:
        q=F(node['q']);expected=LABELS.index(node['winner'])
        scores,entropy,weights=evaluate(precomputed,q)
        required_gaps=[(i,scores[expected]-s) for i,s in enumerate(scores) if i!=expected]
        good=all(g.lo>0 for i,g in required_gaps)
        refutations=[{'rival':LABELS[i],'required_minus_rival':g.pack()} for i,g in required_gaps if g.hi<0]
        actual=[]
        for i,score in enumerate(scores):
            if all((score-other).lo>0 for j,other in enumerate(scores) if i!=j):actual.append(LABELS[i])
        require(len(actual)<=1,'two certified strict winners')
        results.append({'q':str(q),'required_winner':node['winner'],
                        'required_status':'PASS' if good else ('STRICTLY_REFUTED' if refutations else 'ENCLOSURE_INCONCLUSIVE'),
                        'actual_strict_winner':actual[0] if actual else None,
                        'minimum_required_margin_lower':str(min(g.lo for i,g in required_gaps)),
                        'strict_refutations':refutations,
                        'all_scores':[s.pack() for s in scores],
                        'entropy_factors':[x.pack() for x in entropy],
                        'weights':[x.pack() for x in weights]})
    return results


def round_source(source,digits):
    den=10**digits
    return [[[F(round(x*den),den) for x in pair] for pair in row] for row in source]


def main(output):
    raw=SOURCE.read_bytes();require(hashlib.sha256(raw).hexdigest()==SOURCE_SHA,'frozen source file changed')
    data=json.loads(raw)
    canonical=json.dumps(data['source'],sort_keys=True,separators=(',',':')).encode()
    require(hashlib.sha256(canonical).hexdigest()==data['source_sha256'],'embedded exact source hash')
    require(data['baseline_r']=='4' and data['labels']==LABELS,'model baseline and label order')
    nodes=data['q_nodes'];h=F(data['uniform_raw_halfwidth'])
    require(h==F(1,10**11) and len(nodes)==7,'frozen box and node count')
    original=[[[F(x) for x in pair] for pair in row] for row in data['source']]
    kernel.selftests()
    trials=[]
    for digits in range(6,17):
        current=round_source(original,digits)
        if digits==16:require(current==original,'sixteen-place baseline control altered')
        center=check_nodes(current,nodes,F(0))
        box=check_nodes(current,nodes,h)
        center_pass=all(x['required_status']=='PASS' for x in center)
        box_pass=all(x['required_status']=='PASS' for x in box)
        delta=max(abs(x-y) for oldrow,row in zip(original,current) for oldpair,pair in zip(oldrow,row) for x,y in zip(oldpair,pair))
        trial={'digits':digits,'baseline_control_only':digits==16,
               'rounded_source':[[list(map(str,pair)) for pair in row] for row in current],
               'maximum_raw_coordinate_change':str(delta),
               'center_all_seven_required_winners':center_pass,
               'full_raw_box_all_seven_required_winners':box_pass,
               'center_records':center,'box_records':box}
        trials.append(trial)
        print(json.dumps({'digits':digits,'center_pass':center_pass,'box_pass':box_pass,
                          'center_winners':[x['actual_strict_winner'] for x in center],
                          'center_statuses':[x['required_status'] for x in center],
                          'minimum_center_margin_lower':min(F(x['minimum_required_margin_lower']) for x in center).__str__(),
                          'minimum_box_margin_lower':min(F(x['minimum_required_margin_lower']) for x in box).__str__()}),flush=True)
    tested=trials[:-1]
    first_center=next((x['digits'] for x in tested if x['center_all_seven_required_winners']),None)
    first_box=next((x['digits'] for x in tested if x['full_raw_box_all_seven_required_winners']),None)
    require(trials[-1]['center_all_seven_required_winners'] and trials[-1]['full_raw_box_all_seven_required_winners'],'baseline control replay failed')
    require(hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA,'source changed during read-only audit')
    out={'status':'COMPLETE_DECIMAL100_SPECIFIC_ROUNDING_PROTOCOL_AUDIT',
         'source_file_sha256_before_and_after':SOURCE_SHA,'source_internal_sha256':data['source_sha256'],
         'kernel_sha256':KERNEL_SHA,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'rounding':'Exact rational nearest decimal with ties to even; every raw mu/nu coordinate, q nodes unchanged.',
         'tested_precision_range':[6,15],'baseline_control_precision':16,
         'box_halfwidth':str(h),'box_geometry':'All28 baseline raw coordinates independently +/-h, then transport and recompute entropy weights.',
         'first_passing_center_precision_in_tested_protocol':first_center,
         'first_passing_box_precision_in_tested_protocol':first_box,
         'trials':trials,
         'scope':'First tested success under this particular simultaneous rounding protocol and frozen q nodes. Not a universal precision lower bound, exact radius, or optimal decimal encoding.'}
    with output.open('x',encoding='utf8') as stream:json.dump(out,stream,indent=2);stream.write('\n')
    print('FIRST_CENTER',first_center,'FIRST_BOX',first_box,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
