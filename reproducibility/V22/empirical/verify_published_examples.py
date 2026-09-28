"""Replay already-exposed printed numerical examples; no publisher PDF needed."""
import argparse
from pathlib import Path
from datetime import datetime, timezone
from decimal import Decimal as D, getcontext
from functools import reduce
from operator import mul
import hashlib, json, sys
sys.dont_write_bytecode=True
import portable_io
portable_io.verify_module()
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out",type=Path,required=True)
args=parser.parse_args()
OUT=portable_io.new_output(args.out)
HERE=Path(__file__).resolve().parent
def pin(p):
    b=p.read_bytes()
    return {"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()}
pins={"script":pin(Path(__file__)),"scope":"Already-exposed printed numerical examples; no source reauthentication."}
portable_io.save(OUT/"INPUT_PIN.json",pins)
getcontext().prec=90
PI=D('3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066')
def prod(xs):return reduce(mul,xs,D(1))
def cos(x):
    term=total=D(1)
    for k in range(1,200):
        term *= -x*x/D((2*k-1)*(2*k))
        total += term
        if abs(term)<D('1e-100'):return total
    raise ArithmeticError('Cosine failed to converge')
def row(s):return [(D(a),D(b)) for a,b in (v.split(',') for v in s.split())]
def wa(x,w,q):return (1-prod((1-a**q)**z for (a,b),z in zip(x,w)))**(1/q),prod(b**z for (a,b),z in zip(x,w))
def wg(x,w,q):return prod(a**z for (a,b),z in zip(x,w)),(1-prod((1-b**q)**z for (a,b),z in zip(x,w)))**(1/q)
def score(x,q):return x[0]**q-x[1]**q
def compare(cs,ps,tol):
    return [{'computed':str(c),'printed':p,'absolute_error':str(abs(c-D(p))),
             'compatible_half_last_decimal':abs(c-D(p))<=D(tol)} for c,p in zip(cs,ps)]
q=D(3)
basic=row('.7,.5 .3,.5 .6,.7'); bw=list(map(D,['.3','.3','.4']))
basic_result={'Example1_WA':compare(wa(basic,bw,q),['.591','.572'],'.0005'),
              'Example2_WG':compare(wg(basic,bw,q),['.510','.603'],'.0005'),
              'source_pages_prepublication':[9,13],'q':'3','weights':[str(x) for x in bw],
              'input':[[str(a),str(b)] for a,b in basic]}

lw=[row(x) for x in ['.5,.2 .8,.3 .8,.3 .7,.3 .4,.2 .4,.8',
 '.6,.3 .5,.8 .6,.5 .6,.5 .7,.4 .5,.6',
 '.3,.4 .8,.5 .7,.6 .6,.4 .6,.2 .4,.7',
 '.7,.4 .5,.6 .7,.4 .5,.5 .7,.6 .6,.5',
 '.7,.6 .6,.4 .4,.7 .4,.3 .7,.7 .5,.4']]
w=list(map(D,['.20','.10','.30','.15','.15','.10']))
av=[score(wa(x,w,q),q) for x in lw]; gv=[score(wg(x,w,q),q) for x in lw]
ar=compare(av,['.3015','.1184','.1610','.1791','.0404'],'.00005')
gr=compare(gv,['.1701','.0499','.0164','.1349','-.0860'],'.00005')
lw_result={'q':'3','source':'Author full-text pp17-19, TableIII, TableIV and q3 TableVI; no original PDF visual authentication.',
 'input':[[[str(a),str(b)] for a,b in x] for x in lw],'weights':[str(x) for x in w],
 'WA':ar,'WG':gr,'WA_computed_order':[f'X{i+1}' for i in sorted(range(5),key=lambda i:-av[i])],
 'WG_computed_order':[f'X{i+1}' for i in sorted(range(5),key=lambda i:-gv[i])],
 'WG_printed_order':['X1','X4','X2','X3','X5']}

khan=[row(x) for x in ['.8,.1 .3,.4 .5,.3 .2,.6','.7,.2 .6,.4 .5,.2 .3,.3',
 '.5,.5 .6,.3 .7,.1 .8,.2','.3,.5 .5,.4 .4,.4 .6,.2','.6,.2 .7,.1 .3,.5 .5,.4']]
kw=list(map(D,['.22','.20','.33','.25']))
def similarity(x,ideal):
    a0,b0=ideal; h0=1-a0**q-b0**q
    return sum((z*cos(PI/4*(abs(a**q-a0**q)+abs(b**q-b0**q))*(1-abs(1-a**q-b**q-h0)/2)) for (a,b),z in zip(x,kw)),D(0))
dp=[similarity(x,(D(1),D(0))) for x in khan]; dm=[similarity(x,(D(0),D(1))) for x in khan]
cc=[a/(a+b) for a,b in zip(dp,dm)]
kr={'Dplus':compare(dp,['.903296','.922302','.927977','.904749','.908783'],'.0000005'),
    'Dminus':compare(dm,['.842578','.853451','.767745','.885070','.858689'],'.0000005'),
    'closeness':compare(cc,['.517389','.519386','.547246','.505497','.514171'],'.0000005')}
ko=[f'y{i+1}' for i in sorted(range(5),key=lambda i:-cc[i])]
records=basic_result['Example1_WA']+basic_result['Example2_WG']+ar+gr+[r for rr in kr.values() for r in rr]
checks={'all_four_basic_coordinates_compatible':all(r['compatible_half_last_decimal'] for r in records[:4]),
        'all_five_WA_scores_compatible':all(r['compatible_half_last_decimal'] for r in ar),
        'WG_X1_failure_X2toX5_compatible':[r['compatible_half_last_decimal'] for r in gr]==[False,True,True,True,True],
        'all_fifteen_Khan_values_compatible':all(r['compatible_half_last_decimal'] for rr in kr.values() for r in rr),
        'Khan_order':ko==['y3','y2','y1','y5','y4'],
        'source_code_unchanged':pin(Path(__file__))==pins['script']}
result={'created_utc':datetime.now(timezone.utc).isoformat(),'scope':'Published source examples only; no real-data table or sweep read/executed.',
 'implementation':'Standard-library Decimal90; Taylor cosine, 112-decimal supplied pi constant; ordinary high precision, not outward interval certification.',
 'rounding':'Closed half-last-decimal analyst compatibility convention; not a source rounding guarantee.',
 'liu_wang_basic_examples':basic_result,'liu_wang_example4':lw_result,
 'khan_table4':{'q':'3','input':[[[str(a),str(b)] for a,b in x] for x in khan],'weights':[str(x) for x in kw],'records':kr,'order':ko},
 'printed_value_comparisons':len(records),'compatible_values':sum(r['compatible_half_last_decimal'] for r in records),
 'incompatible_values':sum(not r['compatible_half_last_decimal'] for r in records),
 'checks':checks,'gate_checks_all_pass':all(checks.values()),
 'interpretation':'Printed WG Example4 X1 remains incompatible and reverses its top relative to printed q3 order; formula-only benchmark qualification is essential.'}

expected=portable_io.load(HERE/"expected/published_examples.json")
assert {k:v for k,v in result.items() if k!="created_utc"} == {k:v for k,v in expected.items() if k!="created_utc"}, "Archived printed-example numerical result changed"
assert result["gate_checks_all_pass"]
result["status"]="PASS"
result["original_numeric_report_agrees_except_timestamp"]=True
portable_io.save(OUT/"PUBLISHED_EXAMPLE_RESULTS_R2.json",result)
print(json.dumps({k:result[k] for k in ("status","printed_value_comparisons","compatible_values","incompatible_values","gate_checks_all_pass")}))
