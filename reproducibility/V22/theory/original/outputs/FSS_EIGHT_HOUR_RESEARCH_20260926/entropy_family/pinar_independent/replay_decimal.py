"""Independent outward Decimal100 replay via the analytic p=1 TOPSIS reduction.

No imports from the Pinar witness, frozen evaluator, MPFR interval kernel, or
discovery code. Only the pinned source data and output-to-audit are read.
Decimal ln/exp are correctly rounded HALF_EVEN; one neighboring representable
number is added on each side. Other operations use explicit floor/ceiling.
https://docs.python.org/3/library/decimal.html#decimal.Decimal.exp
https://docs.python.org/3/library/decimal.html#decimal.Decimal.ln
"""
from decimal import Decimal as D, Context, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import sys

PREC = 100
DN = Context(prec=PREC,rounding=ROUND_FLOOR)
UP = Context(prec=PREC,rounding=ROUND_CEILING)
RN = Context(prec=PREC,rounding=ROUND_HALF_EVEN)
ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT/'outputs/FSS_V18_RELEASE_20260926/public_artifact/benchmark/case01_pinar_boran_thesis/case_data.json'
SOURCE_SHA = 'e6f8729f9d503e44b5d456cfbeaa6be378721592cc73f2edac40fc82b50f82c9'
MPFR_RECORD = ROOT/'outputs/FSS_EIGHT_HOUR_RESEARCH_20260926/external_case/PINAR_PATH_MPFR256_R2.json'
MPFR_SHA = 'e624df800e884efbbdbc3192e24704531796c5758b5efc378148088cba146058'


def require(condition,message):
    if not condition:
        raise ValueError(message)


class Box:
    def __init__(self,a,b=None):
        a=F(a); b=a if b is None else F(b)
        require(a<=b,'reversed interval')
        self.lo=DN.divide(D(a.numerator),D(a.denominator))
        self.hi=UP.divide(D(b.numerator),D(b.denominator))

    @classmethod
    def direct(cls,a,b):
        require(a.is_finite() and b.is_finite() and a<=b,'invalid endpoint')
        z=object.__new__(cls);z.lo=a;z.hi=b
        return z

    def __add__(self,other):
        other=box(other)
        return Box.direct(DN.add(self.lo,other.lo),UP.add(self.hi,other.hi))
    __radd__=__add__

    def __neg__(self):
        return Box.direct(self.hi.copy_negate(),self.lo.copy_negate())

    def __sub__(self,other):return self+-box(other)
    def __rsub__(self,other):return box(other)+-self

    def __mul__(self,other):
        other=box(other)
        return Box.direct(min(DN.multiply(x,y) for x in [self.lo,self.hi] for y in [other.lo,other.hi]),
                          max(UP.multiply(x,y) for x in [self.lo,self.hi] for y in [other.lo,other.hi]))
    __rmul__=__mul__

    def __truediv__(self,other):
        other=box(other)
        require(other.lo>0,'this checker only divides by positive intervals')
        reciprocal=Box.direct(DN.divide(D(1),other.hi),UP.divide(D(1),other.lo))
        return self*reciprocal
    def __rtruediv__(self,other):return box(other)/self

    def log(self):
        require(self.lo>0,'log nonpositive interval')
        low=RN.ln(self.lo); high=RN.ln(self.hi)
        return Box.direct(RN.next_minus(low),RN.next_plus(high))

    def exp(self):
        low=RN.exp(self.lo); high=RN.exp(self.hi)
        return Box.direct(RN.next_minus(low),RN.next_plus(high))

    def power(self,p):
        return (self.log()*p).exp()

    def contains(self,value):
        return F(self.lo)<=F(value)<=F(self.hi)

    def pack(self):
        return {'lower':str(self.lo),'upper':str(self.hi),
                'exact_endpoints':[str(F(self.lo)),str(F(self.hi))]}


def box(x):return x if isinstance(x,Box) else Box(x)
def bmin(values):return Box.direct(min(x.lo for x in values),min(x.hi for x in values))
def bmax(values):return Box.direct(max(x.lo for x in values),max(x.hi for x in values))
def rational_pair(pair):return tuple(F(str(x)) for x in pair)


def source_and_weights():
    raw=SOURCE.read_bytes()
    require(hashlib.sha256(raw).hexdigest()==SOURCE_SHA,'source hash mismatch')
    data=json.loads(raw)
    require(data['baseline_q']==3 and data['p']==1,'source q/distance exponent')
    require(data['cost_criteria']==['X3'],'cost convention')
    expertise=list(map(rational_pair,data['dm_expertise']))
    importance={key:list(map(rational_pair,val)) for key,val in data['criterion_importance_pairs'].items()}
    scale={key:rational_pair(val) for key,val in data['linguistic_scale'].items()}
    metadata=expertise+[v for vals in importance.values() for v in vals]
    require(len(metadata)==80,'metadata pair count')
    require(all(0<u<1 and 0<v<1 and u**3+v**3<1 for u,v in metadata),'metadata feasibility')
    dwraw=[1+u**3-v**3 for u,v in expertise]
    dw=[v/sum(dwraw) for v in dwraw]
    cwraw=[sum(w*(1+u**3-v**3) for w,(u,v) in zip(dw,importance[f'X{j+1}'])) for j in range(15)]
    cw=[v/sum(cwraw) for v in cwraw]
    require(sum(dw)==sum(cw)==1 and all(x>0 for x in dw+cw),'weight simplex')
    return data,scale,dw,cw


def canonical_matrix(data,scale,dw,cw,lam):
    require(type(lam) is F and 0<=lam<=1,'lambda must be exact rational in[0,1]')
    source=[];changes=[];unchanged=0
    for ai in range(10):
        row=[]
        for j in range(15):
            vals=[]
            for di in range(5):
                key=f'KV{di+1}'
                u,v=scale[data['supplier_linguistic_rows'][key][f'A{ai+1}'][j]]
                original=(u**3,v**3)
                u,v=original
                if ai==7:
                    target=scale[data['supplier_linguistic_rows'][key]['A6'][j]]
                    u=(1-lam)*u+lam*target[0]**3
                    v=(1-lam)*v+lam*target[1]**3
                require(0<u<1 and 0<v<1 and u+v<1,'source orthopair infeasible')
                for ci,(old,new) in enumerate(zip(original,(u,v))):
                    if old!=new:
                        require(ai==7,'changed non-A8 source')
                        changes.append({'alternative':'A8','criterion':f'X{j+1}',
                                        'expert':key,'component':ci,'change':str(new-old)})
                    else:unchanged+=1
                vals.append((u,v))
            row.append(vals)
        source.append(row)
    # Merge expert aggregation and scalar criterion weighting algebraically.
    matrix=[]
    for ai in range(10):
        row=[]
        for j in range(15):
            log_complement=sum(Box(1-u).log()*(dw[di]*cw[j]) for di,(u,v) in enumerate(source[ai][j]))
            log_nonmember=sum(Box(v).log()*(dw[di]*cw[j]) for di,(u,v) in enumerate(source[ai][j]))
            U=1-log_complement.exp(); V=log_nonmember.exp()
            require(U.lo>0 and V.lo>0 and (U+V).hi<1,'weighted pair feasibility')
            row.append((U,V))
        matrix.append(row)
    return matrix,{'checked_pairs':750,'changed_coordinates':len(changes),
                   'unchanged_coordinates':unchanged,'maximum_change':str(max([abs(F(x['change'])) for x in changes],default=F(0))),
                   'coordinate_changes':changes}


def reduced_topsis(matrix,q):
    require(type(q) is int and q>=1,'integer endpoint required')
    t=F(1,q)
    k=F(1,2)-F(833,1000)/(q*q+3*q+1)
    require(0<k<1,'coefficient convexity')
    def utility(U,V):
        return (1-k)*(U.power(t)-V.power(t))+k*((1-V).power(t)-(1-U).power(t))
    low=Box(0);high=Box(0)
    for j in range(15):
        umax=bmax([row[j][0] for row in matrix]);umin=bmin([row[j][0] for row in matrix])
        vmax=bmax([row[j][1] for row in matrix]);vmin=bmin([row[j][1] for row in matrix])
        require((umax+vmin).hi<1 and (umin+vmax).hi<1,'ideal feasibility')
        if j==2:
            low-=utility(umax,vmin);high-=utility(umin,vmax)
        else:
            low+=utility(umin,vmax);high+=utility(umax,vmin)
    totals=[sum((-1 if j==2 else 1)*utility(U,V) for j,(U,V) in enumerate(row)) for row in matrix]
    span=high-low
    require(span.lo>0,'zero ideal span')
    cc=[(total-low)/span for total in totals]
    dp=[(high-total)/30 for total in totals]
    dn=[(total-low)/30 for total in totals]
    return cc,dp,dn,span


def check_rounding(interval,printed,digits):
    center=F(str(printed));radius=F(1,2*10**digits)
    require(F(interval.lo)>center-radius and F(interval.hi)<center+radius,'printed rounding mismatch')


def selftests():
    checks=0
    for a in [F(-7,9),F(-1,3),F(0),F(1,7),F(13,11)]:
        for b in [F(-4,7),F(1,13),F(9,5)]:
            require((Box(a)+Box(b)).contains(a+b),'addition containment')
            require((Box(a)-Box(b)).contains(a-b),'subtraction containment')
            require((Box(a)*Box(b)).contains(a*b),'multiplication containment')
            checks+=3
            if b>0:
                require((Box(a)/Box(b)).contains(a/b),'division containment');checks+=1
    for value in [F(1,1000),F(1,3),F(1),F(8,3),F(15)]:
        require(Box(value).log().exp().contains(value),'exp/log enclosure');checks+=1
    require(Box(0).exp().contains(1),'exp zero')
    return checks+1


def main(output):
    primitive_checks=selftests()
    data,scale,dw,cw=source_and_weights()
    raw=MPFR_RECORD.read_bytes()
    require(hashlib.sha256(raw).hexdigest()==MPFR_SHA,'audited MPFR record hash')
    audited=json.loads(raw)
    matrix,geometry=canonical_matrix(data,scale,dw,cw,F(3887,100000))
    require(geometry['maximum_change']=='2685917/100000000','unexpected source displacement')
    comparisons=[]
    for expected in audited['records']:
        q=expected['q'];winner=int(expected['winner'][1:])-1
        cc,dp,dn,span=reduced_topsis(matrix,q)
        gaps=[cc[winner]-value for i,value in enumerate(cc) if i!=winner]
        require(all(x.lo>0 for x in gaps),'wrong endpoint winner')
        for current,recorded in zip(cc,expected['all_closeness_intervals']):
            lo,hi=map(F,recorded['exact_dyadic_endpoints'])
            require(lo<=F(current.lo)<=F(current.hi)<=hi,'independent interval not inside MPFR enclosure')
        comparisons.append({'q':q,'winner':expected['winner'],
                            'minimum_winner_gap':str(min(x.lo for x in gaps)),
                            'A8_minus_A9':(cc[7]-cc[8]).pack(),
                            'closeness':[x.pack() for x in cc],
                            'all_10_intervals_contained_in_recorded_MPFR':True})
    baseline,zero_geometry=canonical_matrix(data,scale,dw,cw,F(0))
    cc,dp,dn,span=reduced_topsis(baseline,3)
    require(zero_geometry['changed_coordinates']==0,'baseline changed')
    for intervals,key in [(cc,'cc'),(dp,'d_plus'),(dn,'d_minus')]:
        for current,printed in zip(intervals,data['published_q3'][key]):check_rounding(current,printed,3)
    for current,printed in zip(dw,[.2374,.1766,.1766,.2047,.2047]):check_rounding(Box(current),printed,4)
    for current,printed in zip(cw,data['published_criterion_weights_q3']):check_rounding(Box(current),printed,5)
    control,_,_,_=reduced_topsis(baseline,8)
    require(all((control[7]-other).lo>0 for i,other in enumerate(control) if i!=7),'original q8 not A8')
    result={'status':'PASS_INDEPENDENT_DECIMAL100_OUTWARD_REPLAY','source_sha256':SOURCE_SHA,
            'audited_record_sha256':MPFR_SHA,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'precision_decimal_digits':PREC,'python':sys.version,
            'primitive_selftest_checks':primitive_checks,'geometry':geometry,
            'metadata_pairs_admissible':80,'expert_weights_exact':list(map(str,dw)),
            'criterion_weights_exact':list(map(str,cw)),
            'baseline_30_distance_and_closeness_rounding_cells':'PASS',
            'baseline_20_weight_rounding_cells':'PASS','baseline_q8_winner':'A8',
            'records':comparisons,
            'scope':'One exact canonical preaggregation perturbation of A8; endpoint witness only. Not original published centre, calibrated uncertainty, continuum proof, or ambient/source radius.'}
    with output.open('x',encoding='utf8') as stream:
        json.dump(result,stream,indent=2);stream.write('\n')
    print(json.dumps({'status':result['status'],'changed_coordinates':geometry['changed_coordinates'],
                      'maximum_change':geometry['maximum_change'],'records':[
                          {k:row[k] for k in ['q','winner','minimum_winner_gap']} for row in comparisons]},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    main(parser.parse_args().output)
