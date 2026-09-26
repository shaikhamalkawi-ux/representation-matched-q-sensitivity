from pathlib import Path
import csv, json
from mpmath import mp

mp.dps=80
HERE=Path(__file__).resolve().parent
OUT=HERE/'results'; OUT.mkdir(exist_ok=True)
weights=[mp.mpf('0.145'),mp.mpf('0.2'),mp.mpf('0.355'),mp.mpf('0.3')]
X={
'A1':[(.8,.3),(.6,.6),(.4,.8),(.5,.4)],
'A2':[(.3,.8),(.8,.1),(.7,.2),(.4,.4)],
'A3':[(.7,.4),(.6,.5),(.7,.4),(.6,.4)],
'A4':[(.6,.3),(.6,.2),(.7,.4),(.5,.4)],
}
X={a:[(mp.mpf(str(mu)),mp.mpf(str(nu))) for mu,nu in row] for a,row in X.items()}
alts=list(X)
published_q3={
'A1':((.5690,.5423),.0248),'A2':((.6360,.2660),.2384),
'A3':((.6548,.4185),.2074),'A4':((.6178,.3346),.1983)}
published_table2={
4:{'A1':(.5842,.5396),'A2':(.6521,.2644),'A3':(.6562,.4184),'A4':(.6219,.3342)},
6:{'A1':(.6155,.5360),'A2':(.6775,.2630),'A3':(.6595,.4183),'A4':(.6303,.3340)},
8:{'A1':(.6433,.5341),'A2':(.6951,.2625),'A3':(.6628,.4183),'A4':(.6380,.3340)},
10:{'A1':(.6660,.5332),'A2':(.7079,.2623),'A3':(.6659,.4183),'A4':(.6447,.3340)},
15:{'A1':(.7043,.5324),'A2':(.7290,.2621),'A3':(.6726,.4183),'A4':(.6575,.3340)},
20:{'A1':(.7265,.5322),'A2':(.7424,.2621),'A3':(.6777,.4183),'A4':(.6662,.3340)},
}

def transport(x,r,q):
    e=mp.mpf(r)/q
    return x[0]**e,x[1]**e

def ewa(row,q):
    q=mp.mpf(q); pplus=mp.mpf(1); pminus=mp.mpf(1); pn=mp.mpf(1); pden=mp.mpf(1)
    for (mu,nu),w in zip(row,weights):
        pplus *= (1+mu**q)**w
        pminus *= (1-mu**q)**w
        pn *= nu**(w*q)
        pden *= (2-nu**q)**w
    mu=((pplus-pminus)/(pplus+pminus))**(1/q)
    nu=(2*pn/(pden+pn))**(1/q)
    return mu,nu

def score(x,q): return x[0]**q-x[1]**q

def run(q,controlled_from=None):
    agg={}; sc={}
    for a,row in X.items():
        rr=[transport(x,controlled_from,q) for x in row] if controlled_from else row
        agg[a]=ewa(rr,q); sc[a]=score(agg[a],q)
    rank=tuple(sorted(sc,key=lambda a:(-sc[a],a)))
    return agg,sc,rank

q_min=next(q for q in range(1,101) if all(mu**q+nu**q<=1 for row in X.values() for mu,nu in row))

agg3,sc3,rank3=run(3)
q3_pair_err=max(abs(agg3[a][k]-mp.mpf(str(published_q3[a][0][k]))) for a in alts for k in (0,1))
q3_score_err=max(abs(sc3[a]-mp.mpf(str(published_q3[a][1]))) for a in alts)
q3_pass=q3_pair_err<=mp.mpf('0.00005') and q3_score_err<=mp.mpf('0.00005') and rank3==('A2','A3','A4','A1')

table2_checks={}; max_t2=mp.mpf(0)
for q,pairs in published_table2.items():
    agg,sc,rank=run(q)
    err=max(abs(agg[a][k]-mp.mpf(str(pairs[a][k]))) for a in alts for k in (0,1))
    max_t2=max(max_t2,err)
    table2_checks[str(q)]={'within_4dp_rounding_interval':bool(err<=mp.mpf('0.00005')),'max_pair_abs_error':mp.nstr(err,30),'rank':'>'.join(rank)}

rows=[]; raw_trans=[]; ctrl_trans=[]; pr=pc=None
base_ctrl=run(3,controlled_from=3)[1]; max_ctrl=mp.mpf(0)
for q in range(q_min,51):
    ra,rs,rr=run(q); ca,cs,cr=run(q,controlled_from=3)
    if pr is not None and rr!=pr: raw_trans.append({'q':q,'from':'>'.join(pr),'to':'>'.join(rr)})
    if pc is not None and cr!=pc: ctrl_trans.append({'q':q,'from':'>'.join(pc),'to':'>'.join(cr)})
    pr,pc=rr,cr
    max_ctrl=max(max_ctrl,max(abs(cs[a]-base_ctrl[a]) for a in alts))
    rec={'q':q,'raw_rank':'>'.join(rr),'controlled_rank':'>'.join(cr)}
    for a in alts:
        rec[f'raw_score_{a}']=mp.nstr(rs[a],30); rec[f'ctrl_score_{a}']=mp.nstr(cs[a],30)
    rows.append(rec)
with (OUT/'q2_q50_raw_controlled.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

verification={
'case':'Du 2021 q-ROFEWA blockchain traceability design case','doi':'10.3233/JIFS-210548','baseline_q':3,'q_min':q_min,
'q3_reproduction':{'pass':q3_pass,'max_pair_abs_error':mp.nstr(q3_pair_err,30),'max_score_abs_error':mp.nstr(q3_score_err,30)},
'table2_all_published_q_values_reproduced_to_4dp':all(v['within_4dp_rounding_interval'] for v in table2_checks.values()),
'table2_max_pair_abs_error':mp.nstr(max_t2,30),'table2_checks':table2_checks,
'raw_rank_transitions_q2_q50':raw_trans,'controlled_rank_transitions_q2_q50':ctrl_trans,
'max_controlled_score_deviation_vs_q3':mp.nstr(max_ctrl,30),
'class':'CLASS_I_Q_NATURAL_EINSTEIN_AVERAGING_PLUS_Q_SCORE',
'interpretation':'The Einstein aggregation formulas depend on the orthopair through powered coordinates before the final inverse q-root; q-score is powered-coordinate invariant. Controlled invariance is theorem-compatible.',
'decision':'EXECUTABLE_ADMIT'}
(OUT/'verification.json').write_text(json.dumps(verification,indent=2),encoding='utf-8')
print('CASE 5 OVERALL:', 'PASS' if q3_pass and all(v['within_4dp_rounding_interval'] for v in table2_checks.values()) and not ctrl_trans else 'FAIL')
print('q_min =',q_min)
print('q3 max pair error =',mp.nstr(q3_pair_err,12),'score error =',mp.nstr(q3_score_err,12))
print('Table 2 max error =',mp.nstr(max_t2,12))
print('raw transitions =',raw_trans)
print('controlled transitions =',ctrl_trans)
print('max controlled score deviation =',mp.nstr(max_ctrl,12))
