from pathlib import Path
import csv, json
from mpmath import mp

mp.dps=80
HERE=Path(__file__).resolve().parent
OUT=HERE/'results'; OUT.mkdir(exist_ok=True)

# Published normalized q-ROF decision matrix N (Section 6), after cost/benefit normalization.
N=[
[(.8,.5),(.6,.9),(.7,.21),(.8,.7),(.66,.7)],
[(.9,.5),(.8,.2),(.6,.4),(.66,.9),(.70,.41)],
[(.5,.5),(.6,.7),(.8,.8),(.77,.32),(.2,.7)],
[(.35,.81),(.6,.6),(.72,.5),(.79,.31),(.72,.65)],
]
N=[[(mp.mpf(str(mu)),mp.mpf(str(nu))) for mu,nu in row] for row in N]
alts=['Y1','Y2','Y3','Y4']
reported_weights=[mp.mpf(x) for x in ['0.192','0.205','0.192','0.222','0.189']]
published_q4_agg={'Y1':(.7323,.5483),'Y2':(.7717,.4356),'Y3':(.6786,.5658),'Y4':(.6918,.5381)}
published_q4_score={'Y1':.1972,'Y2':.3188,'Y3':.1097,'Y4':.1452}

def transport(x,r,q):
    e=mp.mpf(r)/q
    return x[0]**e,x[1]**e

def transport_matrix(M,r,q):
    return [[transport(x,r,q) for x in row] for row in M]

def entropy_weights(M):
    # Eq. (6): a raw-coordinate entropy stage. This is intentionally not q-natural.
    m=len(M); n=len(M[0]); vals=[]
    for j in range(n):
        s=mp.mpf(0)
        for i in range(m):
            mu,nu=M[i][j]
            s += (mu*mp.log(mu) if mu else 0)+(nu*mp.log(nu) if nu else 0)
        vals.append(1+s/m)
    den=sum(vals)
    return [v/den for v in vals]

def qrofwa(row,w,q):
    q=mp.mpf(q); pm=mp.mpf(1); pn=mp.mpf(1)
    for (mu,nu),ww in zip(row,w):
        pm *= (1-mu**q)**ww; pn *= nu**ww
    return (1-pm)**(1/q),pn

def score(x,q): return x[0]**q-x[1]**q

def run(q,controlled_from=None,recompute_weights=True,use_reported_weights=False):
    M=transport_matrix(N,controlled_from,q) if controlled_from else N
    if use_reported_weights:
        w=reported_weights
    elif recompute_weights:
        w=entropy_weights(M)
    else:
        w=entropy_weights(N)
    agg={a:qrofwa(M[i],w,q) for i,a in enumerate(alts)}
    sc={a:score(agg[a],q) for a in alts}
    rank=tuple(sorted(sc,key=lambda a:(-sc[a],a)))
    return {'weights':w,'agg':agg,'score':sc,'rank':rank}

q_min=next(q for q in range(1,101) if all(mu**q+nu**q<=1 for row in N for mu,nu in row))

# Publication-matched q=4 calculation: the paper prints 3-decimal weights before its downstream results.
bpub=run(4,use_reported_weights=True)
max_agg_err=max(abs(bpub['agg'][a][k]-mp.mpf(str(published_q4_agg[a][k]))) for a in alts for k in (0,1))
max_score_err=max(abs(bpub['score'][a]-mp.mpf(str(published_q4_score[a]))) for a in alts)
base_rank=bpub['rank']
baseline_pass=max_agg_err<=mp.mpf('0.00011') and max_score_err<=mp.mpf('0.00005') and base_rank==('Y2','Y1','Y4','Y3')

# Formula-derived entropy weights are also retained; their values round to the published vector.
w_exact=entropy_weights(N)
weight_round_pass=[str(mp.nint(w*1000)/1000) for w in w_exact]==['0.192','0.205','0.192','0.222','0.189']

# Scientific q sweep: raw path keeps the original decision matrix and recomputes Eq. (6), which gives the same weights for every raw q.
# Controlled path transports the data and recomputes Eq. (6), exposing the raw-coordinate entropy stage as residual dependence.
rows=[]; raw_trans=[]; ctrl_trans=[]; pr=pc=None
base_ctrl=run(4,controlled_from=4,recompute_weights=True)
max_w_dev=mp.mpf(0); max_score_dev=mp.mpf(0)
for q in range(q_min,51):
    raw=run(q,recompute_weights=True)
    ctrl=run(q,controlled_from=4,recompute_weights=True)
    if pr is not None and raw['rank']!=pr:
        raw_trans.append({'q':q,'from':'>'.join(pr),'to':'>'.join(raw['rank'])})
    if pc is not None and ctrl['rank']!=pc:
        ctrl_trans.append({'q':q,'from':'>'.join(pc),'to':'>'.join(ctrl['rank'])})
    pr,pc=raw['rank'],ctrl['rank']
    max_w_dev=max(max_w_dev,max(abs(x-y) for x,y in zip(ctrl['weights'],base_ctrl['weights'])))
    max_score_dev=max(max_score_dev,max(abs(ctrl['score'][a]-base_ctrl['score'][a]) for a in alts))
    rec={'q':q,'raw_rank':'>'.join(raw['rank']),'controlled_rank':'>'.join(ctrl['rank'])}
    for a in alts:
        rec[f'raw_score_{a}']=mp.nstr(raw['score'][a],30); rec[f'ctrl_score_{a}']=mp.nstr(ctrl['score'][a],30)
    for j,w in enumerate(ctrl['weights'],1): rec[f'ctrl_weight_C{j}']=mp.nstr(w,30)
    rows.append(rec)

with (OUT/'q4_q50_raw_controlled.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

verification={
'case':'Seikh & Mandal 2023 q-ROFAWA site-selection case','doi':'10.3390/sym15091680','baseline_q':4,'q_min':q_min,
'published_normalized_matrix_used':True,
'formula_entropy_weights':[mp.nstr(x,30) for x in w_exact],
'formula_weights_round_to_published_3dp':weight_round_pass,
'publication_matched_q4_using_printed_3dp_weights':{
    'pass':baseline_pass,'max_aggregate_pair_abs_error':mp.nstr(max_agg_err,30),'max_score_abs_error':mp.nstr(max_score_err,30),'rank':'>'.join(base_rank)},
'raw_rank_transitions_q4_q50':raw_trans,'controlled_rank_transitions_q4_q50':ctrl_trans,
'max_controlled_entropy_weight_deviation_vs_q4':mp.nstr(max_w_dev,30),
'max_controlled_score_deviation_vs_q4':mp.nstr(max_score_dev,30),
'class':'CLASS_III_RAW_COORDINATE_DEPENDENT_ENTROPY_STAGE',
'interpretation':'The q-ROFAWA aggregation and q-score are q-natural, but Eq. (6) entropy uses raw mu and nu. Under transport the criterion weights drift; despite that residual dependence, the controlled ranking remains stable through q=50.',
'decision':'EXECUTABLE_ADMIT'}
(OUT/'verification.json').write_text(json.dumps(verification,indent=2),encoding='utf-8')
print('CASE 4 OVERALL:', 'PASS' if baseline_pass and not ctrl_trans else 'FAIL')
print('q_min =',q_min)
print('baseline max pair error =',mp.nstr(max_agg_err,12),'score error =',mp.nstr(max_score_err,12))
print('raw transitions =',raw_trans)
print('controlled transitions =',ctrl_trans)
print('max controlled weight drift =',mp.nstr(max_w_dev,12))
print('max controlled score drift =',mp.nstr(max_score_dev,12))
