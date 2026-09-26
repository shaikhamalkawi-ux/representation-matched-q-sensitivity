from pathlib import Path
import json, math, csv

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
DATA = json.loads((HERE / "case_data.json").read_text(encoding="utf-8"))

alts = [f"A{i}" for i in range(1,8)]
criteria = [f"C{i}" for i in range(1,9)]
dmw = DATA["decision_maker_weights"]
ascale = {k:tuple(v) for k,v in DATA["alternative_scale"].items()}
wscale = {k:tuple(v) for k,v in DATA["criterion_scale"].items()}
raw = DATA["table6_rows"]
crit_terms = DATA["table8_criterion_terms"]
cost = {c for c,t in DATA["criterion_types"].items() if t=="cost"}

def qrofwg(vals, weights, q):
    mu = math.prod(v[0]**w for v,w in zip(vals,weights))
    nu = (1.0 - math.prod((1.0-v[1]**q)**w for v,w in zip(vals,weights)))**(1.0/q)
    return mu,nu

def transport(x, r, q):
    return x[0]**(r/q), x[1]**(r/q)

def aggregate_decision(q, transport_from=None):
    out = {}
    for c in criteria:
        terms = raw[c]
        for ai_idx,a in enumerate(alts):
            vals = []
            for d in range(3):
                x = ascale[terms[d*7+ai_idx]]
                if transport_from is not None:
                    x = transport(x, transport_from, q)
                vals.append(x)
            out[(c,a)] = qrofwg(vals, dmw, q)
    return out

def aggregate_criterion_weights(q, transport_from=None):
    out={}
    for c,terms in crit_terms.items():
        vals=[]
        for term in terms:
            x=wscale[term]
            if transport_from is not None:
                x=transport(x,transport_from,q)
            vals.append(x)
        out[c]=qrofwg(vals,dmw,q)
    return out

def qmul(x,y,q):
    mu=x[0]*y[0]
    nu=(x[1]**q+y[1]**q-x[1]**q*y[1]**q)**(1.0/q)
    return mu,nu

def score(x,q): return x[0]**q-x[1]**q
def accuracy(x,q): return x[0]**q+x[1]**q

def ideal_solutions(N,q):
    pos={}; neg={}
    for c in criteria:
        items=[(a,N[(c,a)]) for a in alts]
        key=lambda item:(score(item[1],q),accuracy(item[1],q))
        pos[c]=max(items,key=key)[1]
        neg[c]=min(items,key=key)[1]
    return pos,neg

def entropy(x,q):
    mu,nu=x
    ke=(1.0/math.sqrt(2.0))*math.sqrt(
        (mu**q)**2+(nu**q)**2+(mu**q+nu**q)**2
    )
    return 1.0-ke

def entropy_weights(D,q):
    ent={c:[entropy(D[(c,a)],q) for a in alts] for c in criteria}
    total=sum(sum(vs) for vs in ent.values())
    xi={c:sum(ent[c])/total for c in criteria}
    denom=sum(1.0-xi[c] for c in criteria)
    return {c:(1.0-xi[c])/denom for c in criteria}

def rank_desc(cc):
    return tuple(sorted(alts,key=lambda a:(-cc[a],a)))

def method1(q, transport_from=None):
    D=aggregate_decision(q,transport_from)
    W=aggregate_criterion_weights(q,transport_from)
    R={(c,a):qmul(W[c],D[(c,a)],q) for c in criteria for a in alts}
    N={(c,a):((R[(c,a)][1],R[(c,a)][0]) if c in cost else R[(c,a)])
       for c in criteria for a in alts}
    pos,neg=ideal_solutions(N,q)
    dp={};dn={};cc={}
    n=len(criteria)
    for a in alts:
        sp=sm=0.0
        for c in criteria:
            x=N[(c,a)]; p=pos[c]; m=neg[c]
            sp+=(x[0]**q-p[0]**q)**2+(x[1]**q-p[1]**q)**2
            sm+=(x[0]**q-m[0]**q)**2+(x[1]**q-m[1]**q)**2
        dp[a]=math.sqrt(sp/(2.0*n))
        dn[a]=math.sqrt(sm/(2.0*n))
        cc[a]=dn[a]/(dn[a]+dp[a])
    return {"D":D,"W":W,"N":N,"pos":pos,"neg":neg,
            "d_plus":dp,"d_minus":dn,"cc":cc,"rank":rank_desc(cc)}

def method2(q, transport_from=None, override_weights=None):
    D=aggregate_decision(q,transport_from)
    W=entropy_weights(D,q) if override_weights is None else override_weights
    N={(c,a):((D[(c,a)][1],D[(c,a)][0]) if c in cost else D[(c,a)])
       for c in criteria for a in alts}
    pos,neg=ideal_solutions(N,q)
    dp={};dn={};cc={}
    for a in alts:
        sp=sm=0.0
        for c in criteria:
            x=N[(c,a)]; p=pos[c]; m=neg[c]
            sp+=W[c]*((x[0]**q-p[0]**q)**2+(x[1]**q-p[1]**q)**2)
            sm+=W[c]*((x[0]**q-m[0]**q)**2+(x[1]**q-m[1]**q)**2)
        dp[a]=math.sqrt(sp/2.0)
        dn[a]=math.sqrt(sm/2.0)
        cc[a]=dn[a]/(dn[a]+dp[a])
    return {"D":D,"W":W,"N":N,"pos":pos,"neg":neg,
            "d_plus":dp,"d_minus":dn,"cc":cc,"rank":rank_desc(cc)}

def max_round_error(actual, published, decimals):
    return max(abs(round(actual[a],decimals)-published[i])
               for i,a in enumerate(alts))

# q_min from the declared scales
pairs=list(ascale.values())+list(wscale.values())
q_min=None
for q in range(1,101):
    if all(mu**q+nu**q<=1.0+1e-15 for mu,nu in pairs):
        q_min=q; break

# Baseline reproduction q=5
m1=method1(5)
m2=method2(5)
p1=DATA["published_q5_method1"]
p2=DATA["published_q5_method2"]

baseline_checks={
    "method1_d_plus_printed_precision":
        max_round_error(m1["d_plus"],p1["d_plus"],5)==0,
    "method1_d_minus_printed_precision":
        max_round_error(m1["d_minus"],p1["d_minus"],5)==0,
    "method1_cc_printed_precision":
        max_round_error(m1["cc"],p1["cc"],4)==0,
    "method1_rank_exact":
        list(m1["rank"])==p1["rank"],
    "method2_d_plus_printed_precision":
        max_round_error(m2["d_plus"],p2["d_plus"],5)==0,
    "method2_d_minus_printed_precision":
        max_round_error(m2["d_minus"],p2["d_minus"],5)==0,
    "method2_cc_printed_precision":
        max_round_error(m2["cc"],p2["cc"],5)==0,
    "method2_rank_exact":
        list(m2["rank"])==p2["rank"],
}

# Table 19 full published q sweep
table19_checks={}
for qstr,expect in DATA["published_table19"].items():
    q=int(qstr)
    table19_checks[qstr]={
        "method1":list(method1(q)["rank"])==expect["method1"],
        "method2":list(method2(q)["rank"])==expect["method2"],
    }

# Localized Table 14 inconsistency
formula_w=m2["W"]
printed_w=dict(zip(criteria,p2["table14_weights"]))
max_table14_weight_difference=max(abs(formula_w[c]-printed_w[c]) for c in criteria)
m2_printed=method2(5,override_weights=printed_w)
table14_audit={
    "formula_derived_weights":formula_w,
    "printed_table14_weights":printed_w,
    "max_abs_difference":max_table14_weight_difference,
    "formula_weights_reproduce_table17_18":
        baseline_checks["method2_d_plus_printed_precision"]
        and baseline_checks["method2_d_minus_printed_precision"]
        and baseline_checks["method2_cc_printed_precision"],
    "printed_table14_weights_reproduce_table17_18":
        (max_round_error(m2_printed["d_plus"],p2["d_plus"],5)==0
         and max_round_error(m2_printed["d_minus"],p2["d_minus"],5)==0
         and max_round_error(m2_printed["cc"],p2["cc"],5)==0),
    "interpretation":
        "Localized reporting inconsistency: Eq. (34) applied to the published raw inputs "
        "reproduces Tables 17-18, while the printed Table 14 weight vector does not."
}

# Raw and controlled q=2..50
rows=[]
base1=m1["cc"]; base2=m2["cc"]; basew=m2["W"]
max_ctrl_cc1=max_ctrl_cc2=max_ctrl_w2=0.0
prev_raw1=prev_raw2=prev_ctrl1=prev_ctrl2=None
transitions={"raw_method1":[],"raw_method2":[],"controlled_method1":[],"controlled_method2":[]}

for q in range(q_min,51):
    r1=method1(q)
    r2=method2(q)
    c1=method1(q,transport_from=5)
    c2=method2(q,transport_from=5)
    max_ctrl_cc1=max(max_ctrl_cc1,max(abs(c1["cc"][a]-base1[a]) for a in alts))
    max_ctrl_cc2=max(max_ctrl_cc2,max(abs(c2["cc"][a]-base2[a]) for a in alts))
    max_ctrl_w2=max(max_ctrl_w2,max(abs(c2["W"][c]-basew[c]) for c in criteria))

    ranks = {
        "raw_method1":r1["rank"],"raw_method2":r2["rank"],
        "controlled_method1":c1["rank"],"controlled_method2":c2["rank"]
    }
    for key,cur in ranks.items():
        prev=locals().get("prev_"+key.replace("method",""))
    # explicit transition tracking
    current=[r1["rank"],r2["rank"],c1["rank"],c2["rank"]]
    previous=[prev_raw1,prev_raw2,prev_ctrl1,prev_ctrl2]
    keys=["raw_method1","raw_method2","controlled_method1","controlled_method2"]
    for key,pr,cu in zip(keys,previous,current):
        if pr is not None and pr!=cu:
            transitions[key].append({"q":q,"from":">".join(pr),"to":">".join(cu)})
    prev_raw1,prev_raw2,prev_ctrl1,prev_ctrl2=current

    row={
        "q":q,
        "raw_method1_rank":">".join(r1["rank"]),
        "controlled_method1_rank":">".join(c1["rank"]),
        "raw_method2_rank":">".join(r2["rank"]),
        "controlled_method2_rank":">".join(c2["rank"]),
    }
    for a in alts:
        row[f"raw_m1_cc_{a}"]=r1["cc"][a]
        row[f"ctrl_m1_cc_{a}"]=c1["cc"][a]
        row[f"raw_m2_cc_{a}"]=r2["cc"][a]
        row[f"ctrl_m2_cc_{a}"]=c2["cc"][a]
    rows.append(row)

with (OUT/"q2_q50_raw_controlled_sweep.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)

with (OUT/"q5_formula_vs_table14_weights.csv").open("w",newline="",encoding="utf-8") as f:
    w=csv.writer(f)
    w.writerow(["criterion","Eq34_recomputed","Table14_printed","difference"])
    for c in criteria:
        w.writerow([c,formula_w[c],printed_w[c],formula_w[c]-printed_w[c]])

verification={
    "case":"Alkan & Kahraman 2021 q-ROF TOPSIS",
    "doi":DATA["source"]["doi"],
    "baseline_q":5,
    "q_min":q_min,
    "q1_max_scale_power_sum":max(mu+nu for mu,nu in pairs),
    "q2_max_scale_power_sum":max(mu**2+nu**2 for mu,nu in pairs),
    "baseline_checks":baseline_checks,
    "published_table19_reproduction":table19_checks,
    "table14_audit":table14_audit,
    "raw_transitions_q2_q50": {
        "method1":transitions["raw_method1"],
        "method2":transitions["raw_method2"],
    },
    "controlled_transitions_q2_q50": {
        "method1":transitions["controlled_method1"],
        "method2":transitions["controlled_method2"],
    },
    "max_controlled_cc_deviation_vs_q5":{
        "method1":max_ctrl_cc1,
        "method2":max_ctrl_cc2
    },
    "max_controlled_entropy_weight_deviation_vs_q5":max_ctrl_w2,
    "decision":"ADMIT_WITH_LOCALIZED_TABLE14_REPORTING_INCONSISTENCY",
    "claim_boundary":
        "Both controlled pipelines are theorem-compatible q-natural pipelines. "
        "The Table 14 discrepancy is retained as an audit note and not silently corrected."
}
(OUT/"verification.json").write_text(json.dumps(verification,indent=2),encoding="utf-8")
(OUT/"table14_audit.json").write_text(json.dumps(table14_audit,indent=2),encoding="utf-8")

overall = (
    all(baseline_checks.values())
    and all(v for qd in table19_checks.values() for v in qd.values())
    and not transitions["controlled_method1"]
    and not transitions["controlled_method2"]
)
print("CASE 7 OVERALL:", "PASS" if overall else "FAIL")
print("q_min =", q_min)
print("raw Method I transitions =", transitions["raw_method1"])
print("raw Method II transitions =", transitions["raw_method2"])
print("controlled Method I transitions =", transitions["controlled_method1"])
print("controlled Method II transitions =", transitions["controlled_method2"])
print("max controlled CC deviation M1 =", max_ctrl_cc1)
print("max controlled CC deviation M2 =", max_ctrl_cc2)
print("max controlled entropy-weight deviation =", max_ctrl_w2)
print("Table 14 formula-derived max abs discrepancy =", max_table14_weight_difference)
print("Eq34 weights reproduce Tables 17-18 =", table14_audit["formula_weights_reproduce_table17_18"])
print("Printed Table14 weights reproduce Tables 17-18 =", table14_audit["printed_table14_weights_reproduce_table17_18"])
