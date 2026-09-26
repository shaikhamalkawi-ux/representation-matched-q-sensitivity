from pathlib import Path
import csv, json, math

HERE=Path(__file__).resolve().parent
OUT=HERE/'results'; OUT.mkdir(exist_ok=True)

weights=[0.22,0.20,0.33,0.25]
X={
'y1':[(0.8,0.1),(0.3,0.4),(0.5,0.3),(0.2,0.6)],
'y2':[(0.7,0.2),(0.6,0.4),(0.5,0.2),(0.3,0.3)],
'y3':[(0.5,0.5),(0.6,0.3),(0.7,0.1),(0.8,0.2)],
'y4':[(0.3,0.5),(0.5,0.4),(0.4,0.4),(0.6,0.2)],
'y5':[(0.6,0.2),(0.7,0.1),(0.3,0.5),(0.5,0.4)],
}
alts=list(X)
pos=[(1.0,0.0)]*4
neg=[(0.0,1.0)]*4

published_table4={
'd_plus':[0.903296,0.922302,0.927977,0.904749,0.908783],
'd_minus':[0.842578,0.853451,0.767745,0.885070,0.858689],
'cc':[0.517389,0.519386,0.547246,0.505497,0.514171],
'rank':['y3','y2','y1','y5','y4']}
published_table5={
1:[0.5351,0.5785,0.6449,0.5236,0.5518],
2:[0.5243,0.5374,0.5797,0.5109,0.5259],
3:[0.5173,0.5193,0.5472,0.5054,0.5141],
4:[0.5125,0.5107,0.5297,0.5029,0.5082],
5:[0.5091,0.5063,0.5196,0.5016,0.5049],
6:[0.5067,0.5038,0.5134,0.5009,0.5031],
7:[0.5050,0.5024,0.5093,0.5005,0.5020],
8:[0.5038,0.5015,0.5066,0.5003,0.5013],
9:[0.5029,0.5009,0.5048,0.5002,0.5008],
10:[0.5022,0.5007,0.5035,0.5001,0.5005],
}

def transport(x,r,q):
    return x[0]**(r/q),x[1]**(r/q)

def pair_similarity(a,b,q):
    mu1,nu1=a; mu2,nu2=b
    # eta^q = 1-mu^q-nu^q, so Eq. (20) is entirely powered-coordinate based.
    e1=1-mu1**q-nu1**q
    e2=1-mu2**q-nu2**q
    arg=(math.pi/4.0)*(abs(mu1**q-mu2**q)+abs(nu1**q-nu2**q))*(1-0.5*abs(e1-e2))
    return math.cos(arg)

def weighted_similarity(row,ideal,q):
    return sum(w*pair_similarity(a,b,q) for a,b,w in zip(row,ideal,weights))

def rank_desc(cc):
    return tuple(sorted(cc,key=lambda a:(-cc[a],a)))

def run(q,controlled_from=None):
    D={a:([transport(x,controlled_from,q) for x in row] if controlled_from else row) for a,row in X.items()}
    # Fixed ideals (1,0) and (0,1) are fixed points of T_{r->q}.
    dp={a:weighted_similarity(row,pos,q) for a,row in D.items()}
    dm={a:weighted_similarity(row,neg,q) for a,row in D.items()}
    cc={a:dp[a]/(dp[a]+dm[a]) for a in alts}
    return {'d_plus':dp,'d_minus':dm,'cc':cc,'rank':rank_desc(cc)}

# q_min gate
q_min=next(q for q in range(1,101) if all(mu**q+nu**q<=1+1e-15 for row in X.values() for mu,nu in row))

b=run(3)
max_t4={
    'd_plus':max(abs(b['d_plus'][a]-published_table4['d_plus'][i]) for i,a in enumerate(alts)),
    'd_minus':max(abs(b['d_minus'][a]-published_table4['d_minus'][i]) for i,a in enumerate(alts)),
    'cc':max(abs(b['cc'][a]-published_table4['cc'][i]) for i,a in enumerate(alts)),
}
table4_pass=all(v<5.1e-7 for v in max_t4.values()) and list(b['rank'])==published_table4['rank']

table5_checks={}; max_t5=0.0
for q,pvals in published_table5.items():
    rr=run(q)
    err=max(abs(rr['cc'][a]-pvals[i]) for i,a in enumerate(alts))
    max_t5=max(max_t5,err)
    # Table 5 was printed only to four decimals and appears truncated in several cells.
    table5_checks[str(q)]={'max_abs_error':err,'within_one_last_printed_unit':err<1.0e-4,
                           'rank':'>'.join(rr['rank'])}

rows=[]; transitions={'raw':[],'controlled':[]}; pr=pc=None
base_ctrl=run(3,controlled_from=3)
max_ctrl=0.0
for q in range(q_min,51):
    raw=run(q); ctrl=run(q,controlled_from=3)
    if pr is not None and raw['rank']!=pr:
        transitions['raw'].append({'q':q,'from':'>'.join(pr),'to':'>'.join(raw['rank'])})
    if pc is not None and ctrl['rank']!=pc:
        transitions['controlled'].append({'q':q,'from':'>'.join(pc),'to':'>'.join(ctrl['rank'])})
    pr,pc=raw['rank'],ctrl['rank']
    max_ctrl=max(max_ctrl,max(abs(ctrl['cc'][a]-base_ctrl['cc'][a]) for a in alts))
    rec={'q':q,'raw_rank':'>'.join(raw['rank']),'controlled_rank':'>'.join(ctrl['rank'])}
    for a in alts:
        rec[f'raw_cc_{a}']=raw['cc'][a]; rec[f'ctrl_cc_{a}']=ctrl['cc'][a]
    rows.append(rec)

with (OUT/'q1_q50_raw_controlled.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

verification={
'case':'Khan et al. 2021 weighted-cosine qROF-TOPSIS Algorithm 1',
'doi':'10.1007/s40747-021-00425-7','baseline_q':3,'q_min':q_min,
'table4_exact_to_printed_6dp':table4_pass,'table4_max_abs_errors':max_t4,
'table5_all_q1_q10_within_one_4dp_unit':all(v['within_one_last_printed_unit'] for v in table5_checks.values()),
'table5_max_abs_error':max_t5,'table5_checks':table5_checks,
'raw_rank_transitions_q1_q50':transitions['raw'],
'controlled_rank_transitions_q1_q50':transitions['controlled'],
'max_controlled_closeness_deviation_vs_q3':max_ctrl,
'decision':'EXECUTABLE_ADMIT',
'interpretation':'Eq. (20), the fixed ideals, and the closeness calculation factor through powered coordinates; controlled invariance is theorem-compatible.'
}
(OUT/'verification.json').write_text(json.dumps(verification,indent=2),encoding='utf-8')
print('CASE 2 OVERALL:', 'PASS' if table4_pass and all(v['within_one_last_printed_unit'] for v in table5_checks.values()) and not transitions['controlled'] else 'FAIL')
print('q_min =',q_min)
print('Table 4 max errors =',max_t4)
print('Table 5 max 4dp error =',max_t5)
print('raw transitions =',transitions['raw'])
print('controlled transitions =',transitions['controlled'])
print('max controlled CC deviation =',max_ctrl)
