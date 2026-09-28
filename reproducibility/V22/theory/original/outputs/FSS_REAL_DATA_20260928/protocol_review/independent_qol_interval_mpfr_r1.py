"""Independent directed MPFR256 verification of a frozen QOL p-leaf cover.

Only root cover endpoints are reused; all arithmetic and all source scores are
recomputed from the independently extracted exact rational source closures.
"""
import argparse
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[2]
DEPS=ROOT/'outputs/FSS_STRENGTHENING_20260926/radius_refinement/mpfr_audit/deps'
sys.dont_write_bytecode=True
sys.path.insert(0,str(DEPS))
import gmpy2 as g

BITS=256


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def down(fn):
    with g.context(precision=BITS,round=g.RoundDown,trap_invalid=True,trap_divzero=True,
                   trap_overflow=True,trap_underflow=True):
        x=fn()
        assert g.is_finite(x)
        return x


def up(fn):
    with g.context(precision=BITS,round=g.RoundUp,trap_invalid=True,trap_divzero=True,
                   trap_overflow=True,trap_underflow=True):
        x=fn()
        assert g.is_finite(x)
        return x


def rat(x):
    x=F(x); exact=g.mpq(x.numerator,x.denominator)
    return down(lambda:g.mpfr(exact)),up(lambda:g.mpfr(exact))


def add(a,b):
    return down(lambda:a[0]+b[0]),up(lambda:a[1]+b[1])


def sub(a,b):
    return down(lambda:a[0]-b[1]),up(lambda:a[1]-b[0])


def mul(a,b):
    return (min(down(lambda x=x,y=y:x*y) for x in a for y in b),
            max(up(lambda x=x,y=y:x*y) for x in a for y in b))


def mul_positive(a,b):
    assert a[0]>=0 and b[0]>=0
    return down(lambda:a[0]*b[0]),up(lambda:a[1]*b[1])


def divide_positive(a,b):
    assert a[0]>=0 and b[0]>0
    return down(lambda:a[0]/b[1]),up(lambda:a[1]/b[0])


def logi(a):
    assert a[0]>0
    return down(lambda:g.log(a[0])),up(lambda:g.log(a[1]))


def expi(a):
    return down(lambda:g.exp(a[0])),up(lambda:g.exp(a[1]))


def exact_endpoint(a):
    n,d=a.as_integer_ratio()
    return F(int(n),int(d))


def pack(a):
    assert a[0]<=a[1]
    return [str(exact_endpoint(a[0])),str(exact_endpoint(a[1]))]


class Model:
    def __init__(self,rows,criteria):
        self.rows=rows; self.criteria=criteria
        self.m=len(rows);self.d=len(criteria)
        self.zero=rat(0);self.one=rat(1)
        self.log_source=[];self.score_logs=[]
        for row in rows:
            raw=[];loss=[]
            for c in criteria:
                mu=F(row['criteria'][c]['mu']);nu=F(row['criteria'][c]['nu'])
                assert 0<mu<1 and 0<nu<1 and mu+nu<=1
                raw.append((sub(self.zero,logi(rat(mu))),sub(self.zero,logi(rat(nu)))))
                loss.append((logi(rat(1-mu)),logi(rat(nu))))
            self.log_source.append(raw);self.score_logs.append(loss)
        self.winner=next(i for i,row in enumerate(rows) if row['name']=='3102 Zürich')

    def eval(self,left,right):
        p=(rat(left)[0],rat(right)[1])
        factors=[]
        for j in range(self.d):
            total=self.zero
            for row in self.log_source:
                for rate in row[j]:
                    y=mul_positive(rate,p)
                    atom=mul_positive(y,expi(sub(self.zero,y)))
                    total=add(total,atom)
            H=divide_positive(total,rat(self.m))
            e=sub(self.one,H)
            assert e[0]>0, 'Root accepted leaf fails positive factor with independent enclosure'
            factors.append(e)
        denom=self.zero
        for e in factors:denom=add(denom,e)
        weights=[divide_positive(e,denom) for e in factors]
        scores=[]
        for row in self.score_logs:
            a,b=self.zero,self.zero
            for w,(la,lb) in zip(weights,row):
                a=add(a,mul(w,la));b=add(b,mul(w,lb))
            scores.append(sub(self.one,add(expi(a),expi(b))))
        lower=scores[self.winner][0]
        gaps=[None if i==self.winner else down(lambda i=i:lower-scores[i][1]) for i in range(self.m)]
        assert all(x is None or x>0 for x in gaps), 'Independent MPFR separation failed'
        return {'weights':[pack(w) for w in weights],
                'scores':[pack(s) for s in scores],
                'all_rival_gap_lower':[None if x is None else str(exact_endpoint(x)) for x in gaps],
                'minimum_gap_lower':str(min(exact_endpoint(x) for x in gaps if x is not None))}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--cover',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    assert not args.output.exists(),'Refuse overwrite'
    source_input=BASE/'INDEPENDENT_INPUTS_R1.json'
    freeze=BASE.parent/'qol2023/CONTINUUM_VERIFICATION_FREEZE_R1.md'
    protocol=BASE.parent/'qol2023/PROTOCOL_FREEZE_R1.json'
    source_xlsx=ROOT/'outputs/FSS_AUTONOMOUS_RESEARCH_20260927/data_candidate_QOL2023_R1/quality_of_life_european_cities_2023_aggregated_data.xlsx'
    assert sha(source_input)=='fed011e72692db5822087162301ea91f0876b9c66ce7a03e666384a958925d18'
    assert sha(freeze)=='3e5319463273210e505294e2af33a2791d66316edb6c290ec4546a2648e012be'
    assert sha(protocol)=='3f763601184361e76894a1a926b0536512faf656f76d62e5142f6e74e362dd14'
    assert sha(source_xlsx)=='49986134cde7cddd1858ad1536558e66698ae86d28252f4d9b911ee63d79c03c'
    pins={str(p):sha(p) for p in (source_input,freeze,protocol,source_xlsx,args.cover,Path(__file__))}
    depfiles=[p for folder in (DEPS/'gmpy2',DEPS/'gmpy2.libs') for p in folder.rglob('*') if p.is_file() and p.suffix in ('.py','.pyd','.dll')]
    deppins={str(p.relative_to(DEPS)):sha(p) for p in depfiles}
    assert Path(g.__file__).resolve().is_relative_to(DEPS.resolve())
    source=json.loads(source_input.read_text(encoding='utf-8'))
    plan=json.loads(protocol.read_text(encoding='utf-8'))
    cover=json.loads(args.cover.read_text(encoding='utf-8'))
    groups={}
    for label in ('primary','secondary'):
        criteria=plan[label+'_criteria']
        rows=[r for r in source['rows'] if all(c in r['criteria'] for c in criteria)]
        assert len(rows)==83
        model=Model(rows,criteria)
        leaves=cover['groups'][label]['leaves']
        intervals=[(F(leaf['p_lo']),F(leaf['p_hi'])) for leaf in leaves]
        assert len(intervals)<=4096
        assert intervals==sorted(intervals) and intervals[0][0]==F(1,16) and intervals[-1][1]==1
        assert all(a<b for a,b in intervals)
        assert all(a[1]==b[0] for a,b in zip(intervals,intervals[1:]))
        reviewed=[]
        for n,(left,right) in enumerate(intervals):
            result=model.eval(left,right)
            reviewed.append({'p_lo':str(left),'p_hi':str(right),**result})
            if n%20==0:print(label,n+1,'/',len(intervals),flush=True)
        minimum=min(F(leaf['minimum_gap_lower']) for leaf in reviewed)
        groups[label]={'criteria':criteria,'cities':[r['name'] for r in rows],
                       'q_closed_interval':['1','16'],'p_closed_interval':['1/16','1'],
                       'winner':'3102 Zürich','leaves_count':len(reviewed),
                       'all_82_rivals_positive_on_every_leaf':True,'exact_contiguous_cover':True,
                       'minimum_certified_gap_lower':str(minimum),
                       'minimum_gap_decimal_display_not_certificate':float(minimum),'leaves':reviewed}
    assert all(sha(Path(p))==digest for p,digest in pins.items())
    assert all(sha(DEPS/p)==digest for p,digest in deppins.items())
    out={'record':'Independent fixed-table QOL directed MPFR256 leaf verification R1',
         'created_utc':datetime.now(timezone.utc).isoformat(),'status':'PASS',
         'arithmetic':{'precision_bits':256,'gmpy2':g.version(),'mpfr':g.mpfr_version(),
                       'rounding':'Explicit MPFR RoundDown/RoundUp for exact rational conversion and every arithmetic/log/exp operation; exact dyadic endpoint serialization.'},
         'pins':pins,'dependency_pins':deppins,'protected_inputs_and_dependencies_unchanged':True,
         'groups':groups,
         'scope':'Post-screen verification of one declared model on exact analyst-closed survey estimates over q in [1,16], retaining all83cities and both products. Not all finiteq, new data validation, sampling/input uncertainty, a new theorem, or officialcity ranking.'}
    args.output.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'output':str(args.output),'sha256':sha(args.output),'status':'PASS',
                      'groups':{k:{s:v[s] for s in ('leaves_count','minimum_gap_decimal_display_not_certificate')} for k,v in groups.items()}}),flush=True)


if __name__=='__main__':main()
