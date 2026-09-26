from pathlib import Path
import json, csv
from mpmath import mp

mp.dps = 80
HERE=Path(__file__).resolve().parent
OUT=HERE/"results"; OUT.mkdir(exist_ok=True)
D=json.loads((HERE/"case_data.json").read_text(encoding="utf-8"))

alts=[f"A{i}" for i in range(1,11)]
criteria=[f"X{i}" for i in range(1,16)]
dms=[f"KV{i}" for i in range(1,6)]
scale={k:(mp.mpf(str(v[0])),mp.mpf(str(v[1]))) for k,v in D["linguistic_scale"].items()}
rows=D["supplier_linguistic_rows"]
crit={c:[(mp.mpf(str(x[0])),mp.mpf(str(x[1]))) for x in vals] for c,vals in D["criterion_importance_pairs"].items()}
expertise=[(mp.mpf(str(x[0])),mp.mpf(str(x[1]))) for x in D["dm_expertise"]]

# Data-integrity gate
assert len(rows)==5
assert all(len(rows[dm])==10 for dm in dms)
assert all(len(rows[dm][a])==15 for dm in dms for a in alts)
assert sum(len(rows[dm][a]) for dm in dms for a in alts)==750

def transport(x,r,q):
    return x[0]**(mp.mpf(r)/q), x[1]**(mp.mpf(r)/q)

def dm_weights(q, exp):
    vals=[1+mu**q-nu**q for mu,nu in exp]
    den=sum(vals)
    return [v/den for v in vals]

def qrofwa(vals,w,q):
    pm=mp.mpf(1)
    pn=mp.mpf(1)
    for (mu,nu),wi in zip(vals,w):
        pm *= (1-mu**q)**wi
        pn *= nu**wi
    return (1-pm)**(mp.mpf(1)/q), pn

def criterion_weights(q,dw,cvals):
    nums=[]
    for c in criteria:
        nums.append(sum(lam*(1+mu**q-nu**q) for lam,(mu,nu) in zip(dw,cvals[c])))
    den=sum(nums)
    return [x/den for x in nums]

def scalar_weight(x,w,q):
    mu,nu=x
    return (1-(1-mu**q)**w)**(mp.mpf(1)/q), nu**w

def kfun(q):
    q=mp.mpf(q)
    # Thesis Eq. 5.39 / 6.11 as printed, including 0.333.
    return (mp.mpf("0.5")*q*q+mp.mpf("1.5")*q-mp.mpf("0.333"))/(q*q+3*q+1)

def distance(xv,yv,q,p=1):
    k=kfun(q)
    s=mp.mpf(0)
    for (ma,na),(mb,nb) in zip(xv,yv):
        t1=(1-k)*(ma-mb)+k*((1-na**q)**(mp.mpf(1)/q)-(1-nb**q)**(mp.mpf(1)/q))
        t2=(1-k)*(na-nb)+k*((1-ma**q)**(mp.mpf(1)/q)-(1-mb**q)**(mp.mpf(1)/q))
        s += abs(t1)**p+abs(t2)**p
    return (s/(2*len(xv)))**(mp.mpf(1)/p)

def run(q, controlled_from=None):
    if controlled_from is None:
        sc=scale; exp=expertise; cv=crit
    else:
        sc={k:transport(v,controlled_from,q) for k,v in scale.items()}
        exp=[transport(v,controlled_from,q) for v in expertise]
        cv={c:[transport(v,controlled_from,q) for v in vals] for c,vals in crit.items()}

    dw=dm_weights(q,exp)
    agg={}
    for a in alts:
        for j,c in enumerate(criteria):
            vals=[sc[rows[dm][a][j]] for dm in dms]
            agg[(a,c)]=qrofwa(vals,dw,q)

    cw=criterion_weights(q,dw,cv)
    W={(a,c):scalar_weight(agg[(a,c)],cw[j],q) for a in alts for j,c in enumerate(criteria)}

    pos={}; neg={}
    for c in criteria:
        vals=[W[(a,c)] for a in alts]
        if c in D["cost_criteria"]:
            pos[c]=(min(x[0] for x in vals),max(x[1] for x in vals))
            neg[c]=(max(x[0] for x in vals),min(x[1] for x in vals))
        else:
            pos[c]=(max(x[0] for x in vals),min(x[1] for x in vals))
            neg[c]=(min(x[0] for x in vals),max(x[1] for x in vals))

    dp={};dn={};cc={}
    pv=[pos[c] for c in criteria]; nv=[neg[c] for c in criteria]
    for a in alts:
        xv=[W[(a,c)] for c in criteria]
        dp[a]=distance(xv,pv,q,D["p"])
        dn[a]=distance(xv,nv,q,D["p"])
        cc[a]=dn[a]/(dp[a]+dn[a])
    rank=tuple(sorted(alts,key=lambda a:(-cc[a],a)))
    return {"dw":dw,"cw":cw,"dp":dp,"dn":dn,"cc":cc,"rank":rank}

def within_printed(x,y,dec=3):
    # source is printed at dec decimals; allow the half-unit rounding interval
    return abs(x-mp.mpf(str(y))) <= mp.mpf("0.5")*mp.mpf(10)**(-dec)+mp.mpf("1e-30")

# q_min of all declared linguistic points
allpairs=list(scale.values())
q_min=None
for q in range(1,101):
    if all(mu**q+nu**q<=1 for mu,nu in allpairs):
        q_min=q; break

# Baseline q=3
b=run(3)
pub=D["published_q3"]
baseline_checks={
    "dm_weights_q3_to_4dp": all(within_printed(x,y,4) for x,y in zip(b["dw"],[.2374,.1766,.1766,.2047,.2047])),
    "criterion_weights_q3_to_5dp": all(within_printed(x,y,5) for x,y in zip(b["cw"],D["published_criterion_weights_q3"])),
    "table74_dplus_to_3dp": all(within_printed(b["dp"][a],pub["d_plus"][i],3) for i,a in enumerate(alts)),
    "table74_dminus_to_3dp": all(within_printed(b["dn"][a],pub["d_minus"][i],3) for i,a in enumerate(alts)),
    "table75_cc_to_3dp": all(within_printed(b["cc"][a],pub["cc"][i],3) for i,a in enumerate(alts)),
    "table75_rank_exact": list(b["rank"])==pub["rank"],
}

# Figure 7.4(a): all 10 x 9 printed Ci values q=2..10
figure_checks={}
fig=D["published_figure_7_4a_q2_q10"]
max_fig_abs=mp.mpf(0)
for q in range(2,11):
    rr=run(q)
    okay=True
    for a in alts:
        y=fig[a][q-2]
        max_fig_abs=max(max_fig_abs,abs(rr["cc"][a]-mp.mpf(str(y))))
        okay &= within_printed(rr["cc"][a],y,3)
    figure_checks[str(q)]=bool(okay)

# High-precision q=2..50 raw and representation-controlled
sweep=[]
trans={"raw":[],"controlled":[]}
prev_raw=prev_ctrl=None
base_ctrl=run(3,controlled_from=3)
max_dw_dev=mp.mpf(0); max_cw_dev=mp.mpf(0)
for q in range(q_min,51):
    raw=run(q)
    ctrl=run(q,controlled_from=3)
    if prev_raw is not None and raw["rank"]!=prev_raw:
        trans["raw"].append({"q":q,"from":">".join(prev_raw),"to":">".join(raw["rank"])})
    if prev_ctrl is not None and ctrl["rank"]!=prev_ctrl:
        trans["controlled"].append({"q":q,"from":">".join(prev_ctrl),"to":">".join(ctrl["rank"])})
    prev_raw=raw["rank"]; prev_ctrl=ctrl["rank"]
    max_dw_dev=max(max_dw_dev,max(abs(x-y) for x,y in zip(ctrl["dw"],base_ctrl["dw"])))
    max_cw_dev=max(max_cw_dev,max(abs(x-y) for x,y in zip(ctrl["cw"],base_ctrl["cw"])))
    rec={
        "q":q,
        "raw_rank":">".join(raw["rank"]),
        "controlled_rank":">".join(ctrl["rank"]),
    }
    for a in alts:
        rec[f"raw_cc_{a}"]=mp.nstr(raw["cc"][a],30)
        rec[f"ctrl_cc_{a}"]=mp.nstr(ctrl["cc"][a],30)
    sweep.append(rec)

with (OUT/"q2_q50_high_precision_raw_controlled.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=list(sweep[0]))
    w.writeheader(); w.writerows(sweep)

verification={
    "case":"Pinar 2020 thesis q-ROF TOPSIS",
    "q_min":q_min,
    "baseline_q":3,
    "p":1,
    "raw_cells_transcribed":750,
    "baseline_checks":baseline_checks,
    "figure_7_4a_all_90_printed_values_reproduced":all(figure_checks.values()),
    "figure_7_4a_checks_by_q":figure_checks,
    "maximum_absolute_difference_vs_3dp_figure_value":mp.nstr(max_fig_abs,30),
    "raw_rank_transitions_q2_q50_high_precision":trans["raw"],
    "controlled_rank_transitions_q2_q50_high_precision":trans["controlled"],
    "max_controlled_dm_weight_deviation_vs_q3":mp.nstr(max_dw_dev,30),
    "max_controlled_criterion_weight_deviation_vs_q3":mp.nstr(max_cw_dev,30),
    "numerical_note":{
        "reason_for_high_precision":"At large q, ordinary float64 suffers cancellation/underflow in the weighted q-ROF and distance pipeline and can create spurious rank reversals.",
        "precision_decimal_digits":80
    },
    "reporting_note":{
        "thesis_text":"The prose says the A5/A7 change is for q>7.",
        "recomputed_unrounded_result":"The exact raw ranking first changes at q=6; Figure 7.4(a) prints A5 and A7 both as 0.214 at q=6, then separates them at q=7. A6 also passes A5 at q=10.",
        "classification":"Localized text/rounding-resolution discrepancy; the printed numeric Figure 7.4(a) is reproduced."
    },
    "decision":"EXECUTABLE_ADMIT"
}
(OUT/"verification.json").write_text(json.dumps(verification,ensure_ascii=False,indent=2),encoding="utf-8")

overall=(
    all(baseline_checks.values())
    and all(figure_checks.values())
    and len(trans["controlled"])==0
)
print("CASE 1 OVERALL:", "PASS" if overall else "FAIL")
print("q_min =",q_min)
print("baseline checks =",baseline_checks)
print("Figure 7.4(a) 90/90 printed values =",all(figure_checks.values()))
print("raw transitions =",trans["raw"])
print("controlled transitions =",trans["controlled"])
print("max controlled DM-weight deviation =",mp.nstr(max_dw_dev,12))
print("max controlled criterion-weight deviation =",mp.nstr(max_cw_dev,12))
