from pathlib import Path
import csv, json
from mpmath import mp

mp.dps=80
HERE=Path(__file__).resolve().parent
OUT=HERE/'results'; OUT.mkdir(exist_ok=True)
weights=[mp.mpf('0.20'),mp.mpf('0.10'),mp.mpf('0.30'),mp.mpf('0.15'),mp.mpf('0.15'),mp.mpf('0.10')]
X={
'X1':[(.5,.2),(.8,.3),(.8,.3),(.7,.3),(.4,.2),(.4,.8)],
'X2':[(.6,.3),(.5,.8),(.6,.5),(.6,.5),(.7,.4),(.5,.6)],
'X3':[(.3,.4),(.8,.5),(.7,.6),(.6,.4),(.6,.2),(.4,.7)],
'X4':[(.7,.4),(.5,.6),(.7,.4),(.5,.5),(.7,.6),(.6,.5)],
'X5':[(.7,.6),(.6,.4),(.4,.7),(.4,.3),(.7,.7),(.5,.4)],
}
X={a:[(mp.mpf(str(mu)),mp.mpf(str(nu))) for mu,nu in row] for a,row in X.items()}
alts=list(X)

published_table4_q3={'X1':.3015,'X2':.1184,'X3':.1610,'X4':.1791,'X5':.0404}
published_table5={
2:[.3746,.1443,.1839,.2064,.0356],
3:[.3015,.1184,.1610,.1791,.0404],
5:[.1755,.0609,.0972,.1044,.0336],
8:[.0799,.0183,.0391,.0381,.0162],
10:[.0487,.0079,.0213,.0188,.0088],
12:[.0302,.0034,.0118,.0092,.0046],
15:[.0149,.0010,.0051,.0031,.0016],
20:[.0047,.0001,.0014,.0005,.0003],
}
# Table VI q-ROFWG values. They are used only as an audit check.
published_wg={
2:[.2480,.0758,.0352,.1644,-.0883],
3:[.1701,.0499,.0164,.1349,-.0860],
5:[.0627,.0106,-.0021,.0694,-.0591],
8:[.0081,-.0066,-.0043,.0205,-.0246],
10:[-.0007,-.0071,-.0027,.0086,-.0126],
12:[-.0028,-.0056,-.0015,.0036,-.0063],
15:[-.0025,-.0033,-.0005,.0009,-.0022],
}

def transport(x,r,q):
    e=mp.mpf(r)/q
    return x[0]**e,x[1]**e

def wa(row,q):
    q=mp.mpf(q); pm=mp.mpf(1); pn=mp.mpf(1)
    for (mu,nu),w in zip(row,weights):
        pm *= (1-mu**q)**w; pn *= nu**w
    return (1-pm)**(1/q),pn

def wg(row,q):
    q=mp.mpf(q); pm=mp.mpf(1); pn=mp.mpf(1)
    for (mu,nu),w in zip(row,weights):
        pm *= mu**w; pn *= (1-nu**q)**w
    return pm,(1-pn)**(1/q)

def score(x,q): return x[0]**q-x[1]**q

def run(q,controlled_from=None):
    sc={}
    for a,row in X.items():
        rr=[transport(x,controlled_from,q) for x in row] if controlled_from else row
        sc[a]=score(wa(rr,q),q)
    rank=tuple(sorted(sc,key=lambda a:(-sc[a],a)))
    return sc,rank

q_min=next(q for q in range(1,101) if all(mu**q+nu**q<=1 for row in X.values() for mu,nu in row))

q3,rank3=run(3)
q3_pass=all(abs(q3[a]-mp.mpf(str(published_table4_q3[a])))<=mp.mpf('0.00005') for a in alts) and rank3==('X1','X4','X3','X2','X5')

table5_checks={}; max_table5=mp.mpf(0)
for q,pvals in published_table5.items():
    sc,rank=run(q)
    err=max(abs(sc[a]-mp.mpf(str(pvals[i]))) for i,a in enumerate(alts))
    max_table5=max(max_table5,err)
    table5_checks[str(q)]={'within_4dp_rounding_interval':bool(err<=mp.mpf('0.00005')),'max_abs_error':mp.nstr(err,30),'rank':'>'.join(rank)}

# Standard printed q-ROFWG formula/data audit. X2-X5 reproduce; X1 does not.
wg_audit={}
for q,pvals in published_wg.items():
    calc={a:score(wg(row,q),q) for a,row in X.items()}
    wg_audit[str(q)]={a:{'computed':mp.nstr(calc[a],30),'published':pvals[i],'abs_error':mp.nstr(abs(calc[a]-mp.mpf(str(pvals[i]))),30)} for i,a in enumerate(alts)}

rows=[]; raw_trans=[]; ctrl_trans=[]; pr=pc=None
base_ctrl=run(3,controlled_from=3)[0]; max_ctrl=mp.mpf(0)
for q in range(q_min,51):
    raw,rr=run(q); ctrl,cr=run(q,controlled_from=3)
    if pr is not None and rr!=pr: raw_trans.append({'q':q,'from':'>'.join(pr),'to':'>'.join(rr)})
    if pc is not None and cr!=pc: ctrl_trans.append({'q':q,'from':'>'.join(pc),'to':'>'.join(cr)})
    pr,pc=rr,cr
    max_ctrl=max(max_ctrl,max(abs(ctrl[a]-base_ctrl[a]) for a in alts))
    rec={'q':q,'raw_rank':'>'.join(rr),'controlled_rank':'>'.join(cr)}
    for a in alts:
        rec[f'raw_score_{a}']=mp.nstr(raw[a],30); rec[f'ctrl_score_{a}']=mp.nstr(ctrl[a],30)
    rows.append(rec)
with (OUT/'q2_q50_raw_controlled_qrofwa.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

verification={
'case':'Liu & Wang 2018 q-ROFWA Example 4','doi':'10.1002/int.21927','correction_doi':'10.1155/int/9864340',
'baseline_q':3,'q_min':q_min,'q3_table4_score_reproduction':q3_pass,
'table5_all_published_rows_reproduced_to_4dp':all(v['within_4dp_rounding_interval'] for v in table5_checks.values()),
'table5_max_abs_error':mp.nstr(max_table5,30),'table5_checks':table5_checks,
'raw_rank_transitions_q2_q50':raw_trans,'controlled_rank_transitions_q2_q50':ctrl_trans,
'max_controlled_score_deviation_vs_q3':mp.nstr(max_ctrl,30),
'qrofwg_audit':{'status':'HOLD_LOCALIZED_X1_MISMATCH','details':wg_audit,
'note':'The standard printed q-ROFWG operator applied to Table III reproduces X2-X5 at the published q values but not X1. q-ROFWG is therefore excluded from the clean benchmark endpoint.'},
'correction_scope_note':'The 2026 correction changes Definition 1 indeterminacy degree only; it does not replace the q-ROFWA operator, Example 4 data, or Table V.',
'decision':'EXECUTABLE_ADMIT_QROFWA_ONLY'}
(OUT/'verification.json').write_text(json.dumps(verification,indent=2),encoding='utf-8')
(OUT/'qrofwg_mismatch_audit.json').write_text(json.dumps(verification['qrofwg_audit'],indent=2),encoding='utf-8')
print('CASE 3 OVERALL:', 'PASS' if q3_pass and all(v['within_4dp_rounding_interval'] for v in table5_checks.values()) and not ctrl_trans else 'FAIL')
print('q_min =',q_min)
print('Table V max error =',mp.nstr(max_table5,12))
print('raw transitions =',raw_trans)
print('controlled transitions =',ctrl_trans)
print('max controlled score deviation =',mp.nstr(max_ctrl,12))
print('q-ROFWG audit = HOLD_LOCALIZED_X1_MISMATCH')
