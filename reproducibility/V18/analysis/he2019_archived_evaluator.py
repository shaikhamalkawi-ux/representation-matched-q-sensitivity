#!/usr/bin/env python3
"""He et al. (2019) q-rung picture-fuzzy Dombi-Hamy evaluator reconstruction.

Scope:
- Reconstructs the table-consistent Step-2 q-RPFDWA aggregation from the
  three published expert matrices (p=3 balancing factor times expert weights).
- Reconstructs Step-3 q-RPFDWHM (k=2, Dombi lambda=2) and the score.
- Reproduces the locked q=3 baseline, q=1..10 raw transition, continuous
  A2/A5 crossing, and matched-control invariance.
- Computes local conditioning and a numerically found feasible decision-boundary
  point in canonical coordinates. The numerical boundary distance is an upper
  bound on the true global reversal radius, NOT a certified global optimum.

The source has an internal formula/example ambiguity for q-RPFDWHM(k=1):
literal Eq. (21) reproduces Example 2 membership but not the printed eta/v.
For the project application, the table-consistent Step-2 implementation below
reproduces 59/60 Table-4 components at 3 decimals; the sole mismatch is
0.089464623... vs printed 0.090, consistent with the retained audit note.
"""
import itertools, math, csv, json
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, brentq

OUT=Path(__file__).resolve().parent

E1=np.array([
[[.5,.4,.1],[.8,.1,.1],[.4,.3,.2],[.1,.8,.1]],
[[.7,.1,.1],[.1,.7,.2],[.1,.7,.2],[.7,.1,.1]],
[[.8,.1,.1],[.1,.8,.1],[.1,.8,.1],[.6,.2,.1]],
[[.7,.1,.1],[.7,.2,.1],[.1,.7,.1],[.1,.8,.1]],
[[.7,.2,.1],[.6,.2,.1],[.8,.1,.1],[.1,.7,.1]]],float)
E2=np.array([
[[.5,.3,.1],[.8,.1,.1],[.5,.1,.3],[.1,.8,.1]],
[[.6,.1,.2],[.2,.5,.2],[.2,.6,.1],[.6,.2,.1]],
[[.8,.1,.1],[.1,.7,.1],[.2,.6,.1],[.6,.1,.2]],
[[.8,.1,.1],[.7,.1,.2],[.2,.6,.1],[.5,.3,.1]],
[[.6,.1,.1],[.8,.1,.1],[.6,.3,.1],[.1,.8,.1]]],float)
E3=np.array([
[[.7,.2,.1],[.3,.5,.1],[.8,.1,.1],[.5,.3,.1]],
[[.1,.7,.1],[.5,.2,.2],[.8,.1,.1],[.7,.1,.2]],
[[.5,.2,.2],[.3,.5,.1],[.1,.7,.1],[.6,.2,.1]],
[[.7,.2,.1],[.2,.7,.1],[.5,.2,.3],[.3,.5,.1]],
[[.5,.2,.1],[.6,.2,.1],[.7,.1,.1],[.2,.7,.1]]],float)
EXPERTS=np.stack([E1,E2,E3])
EW=np.array([.35,.20,.45])
CW=np.array([.4,.2,.3,.1])
TABLE4=np.array([
[[.729,.189,.083],[.831,.092,.083],[.820,.090,.095],[.522,.286,.083]],
[[.715,.092,.086],[.522,.190,.167],[.819,.095,.089],[.769,.086,.092]],
[[.832,.092,.092],[.315,.475,.083],[.186,.596,.083],[.686,.108,.086]],
[[.811,.092,.083],[.738,.108,.086],[.522,.190,.092],[.469,.322,.083]],
[[.721,.108,.083],[.782,.108,.083],[.822,.086,.083],[.211,.627,.083]]],float)
LOCKED_FINAL=np.array([
[.833894,.096167,.067670],[.822747,.082145,.078385],
[.747109,.102921,.067386],[.759373,.107586,.067838],
[.809737,.092902,.065290]],float)
LOCKED_SCORES=np.array([1.57956258,1.55644599,1.41670944,1.43757840,1.53064573])
LOCKED_RAW={q:('A1>A2>A5>A4>A3' if q<=6 else 'A1>A5>A2>A4>A3') for q in range(1,11)}
BASE_RANK='A1>A2>A5>A4>A3'


def rpfdwa_table_consistent(vals,w,q=3.0,lam=2.0):
    """Dombi weighted average used in the project application.

    The source Table 4 is reproduced by the p*w balancing coefficient (p=len(vals)).
    In powered coordinates this is the Dombi weighted-average addition law.
    """
    x=np.asarray(vals)**q
    p=len(vals)
    out=np.empty(3,dtype=np.result_type(vals,float))
    odds=x[:,0]/(1-x[:,0])
    s=np.sum(p*w*odds**lam)
    U=1-1/(1+s**(1/lam)); out[0]=U**(1/q)
    for d in (1,2):
        inv=(1-x[:,d])/x[:,d]
        s=np.sum(p*w*inv**lam)
        X=1/(1+s**(1/lam)); out[d]=X**(1/q)
    return out


def rpfdwhm(vals,w,q=3.0,lam=2.0,k=2):
    vals=np.asarray(vals)
    n=len(vals); combs=list(itertools.combinations(range(n),k)); out=[]
    for d in range(3):
        x=vals[:,d]**q
        ratio=(1-x)/x if d==0 else x/(1-x)
        terms=[]
        for comb in combs:
            inner=sum(w[j]*ratio[j]**lam for j in comb)/k
            terms.append(1/inner)
        A=np.mean(terms)
        X=1-1/(1+A**(1/lam)) if d==0 else 1/(1+A**(1/lam))
        out.append(X**(1/q))
    return np.array(out)


def full_run(q,matched=False,r=3.0,lam=2.0,k=2):
    X=EXPERTS**(r/q) if matched else EXPERTS
    collective=np.empty((5,4,3))
    for a in range(5):
        for c in range(4):
            collective[a,c]=rpfdwa_table_consistent(X[:,a,c,:],EW,q,lam)
    final=np.array([rpfdwhm(collective[a],CW,q,lam,k) for a in range(5)])
    score=final[:,0]**q + 1-final[:,2]**q
    order=np.argsort(-score)
    rank='>'.join(f'A{i+1}' for i in order)
    return collective,final,score,rank

# Canonical evaluator, q-independent for the matched q-natural implementation.
COMBS=list(itertools.combinations(range(4),2))
def canon_alt_score(Zalt,lam=2.0):
    coll=np.empty((4,3),dtype=np.result_type(Zalt))
    p=3.0
    for c in range(4):
        x=Zalt[:,c,:]
        odds=x[:,0]/(1-x[:,0]); s=np.sum(p*EW*odds**lam)
        coll[c,0]=1-1/(1+s**(1/lam))
        for d in (1,2):
            inv=(1-x[:,d])/x[:,d]; s=np.sum(p*EW*inv**lam)
            coll[c,d]=1/(1+s**(1/lam))
    out=np.empty(3,dtype=np.result_type(Zalt))
    for d in range(3):
        x=coll[:,d]; ratio=(1-x)/x if d==0 else x/(1-x)
        terms=[]
        for comb in COMBS:
            inner=sum(CW[j]*ratio[j]**lam for j in comb)/2.0
            terms.append(1/inner)
        A=sum(terms)/len(terms)
        out[d]=1-1/(1+A**(1/lam)) if d==0 else 1/(1+A**(1/lam))
    return out[0]+1-out[2]

Z0=EXPERTS**3

def _uv(za): return za[:,:,[0,2]].ravel()
def _inject(xuv,za0):
    za=za0.copy().astype(np.result_type(xuv)); uv=xuv.reshape(3,4,2)
    za[:,:,0]=uv[:,:,0]; za[:,:,2]=uv[:,:,1]; return za

def _margin(x,j,z1,zj):
    a=_inject(x[:24],z1); b=_inject(x[24:],zj)
    return canon_alt_score(a)-canon_alt_score(b)

def _csjac(x,j,z1,zj,h=1e-30):
    g=np.empty_like(x); xc=x.astype(complex)
    for t in range(len(x)):
        xx=xc.copy(); xx[t]+=1j*h; g[t]=np.imag(_margin(xx,j,z1,zj))/h
    return g

def _simp(x,z1,zj):
    a=x[:24].reshape(3,4,2); b=x[24:].reshape(3,4,2)
    et=np.r_[z1[:,:,1].ravel(),zj[:,:,1].ravel()]
    sm=np.r_[a.sum(2).ravel(),b.sum(2).ravel()]
    return 1-et-sm

def _simp_jac(x,z1,zj):
    J=np.zeros((24,48))
    for t in range(24):
        base=2*t if t<12 else 24+2*(t-12); J[t,base:base+2]=-1
    return J

def solve_boundary(j,nstarts=4):
    z1=Z0[:,0,:,:]; zj=Z0[:,j,:,:]; x0=np.r_[_uv(z1),_uv(zj)]
    g0=_margin(x0,j,z1,zj); grad=_csjac(x0,j,z1,zj); lin=g0/np.linalg.norm(grad)
    xlin=x0-g0*grad/np.dot(grad,grad)
    cons=[{'type':'eq','fun':lambda x:_margin(x,j,z1,zj),'jac':lambda x:_csjac(x,j,z1,zj)},
          {'type':'ineq','fun':lambda x:_simp(x,z1,zj),'jac':lambda x:_simp_jac(x,z1,zj)}]
    def fun(x): d=x-x0; return .5*np.dot(d,d)
    def jac(x): return x-x0
    rng=np.random.default_rng(100+j); inits=[xlin]
    for _ in range(nstarts-1):
        y=xlin+rng.normal(0,.0015,48)
        for __ in range(12):
            if _simp(y,z1,zj).min()>=0 and y.min()>1e-10 and y.max()<1-1e-10: break
            y=(y+x0)/2
        inits.append(np.clip(y,1e-10,1-1e-10))
    sols=[]
    for y in inits:
        rr=minimize(fun,y,jac=jac,method='SLSQP',bounds=[(1e-10,1-1e-10)]*48,
                    constraints=cons,options={'ftol':1e-13,'maxiter':1000,'disp':False})
        sols.append(rr)
    best=min(sols,key=lambda r:r.fun if abs(_margin(r.x,j,z1,zj))<1e-7 else 1e99)
    return {'j':j,'margin':g0,'grad_norm':np.linalg.norm(grad),'rho_lin':lin,
            'rho_found':math.sqrt(2*best.fun),'tie_residual':_margin(best.x,j,z1,zj),
            'simplex_slack':_simp(best.x,z1,zj).min(),'x0':x0,'x':best.x,'success':best.success}


def main():
    # Baseline / source-table verification
    coll,fin,sc,rank=full_run(3)
    rounded=np.round(coll,3)
    mismatch=np.argwhere(rounded!=TABLE4)
    assert len(mismatch)==1 and tuple(mismatch[0])==(0,2,1)
    assert abs(coll[0,2,1]-0.08946462325927318)<1e-12
    assert rank==BASE_RANK
    # q sweep and continuous crossing
    sweep=[]
    matched_scores=[]
    base_match=full_run(3,True)[2]
    for q in range(1,11):
        _,_,sraw,rraw=full_run(q,False)
        _,_,smat,rmat=full_run(q,True)
        assert rraw==LOCKED_RAW[q]
        assert rmat==BASE_RANK
        sweep.append([q,rraw,rmat,float(sraw[0]),float(sraw[1]),float(sraw[4]),float(np.max(np.abs(smat-base_match)))])
        matched_scores.append(smat)
    root=brentq(lambda qq: full_run(qq,False)[2][1]-full_run(qq,False)[2][4],6,7,xtol=1e-14)
    # Boundary candidates
    bnds=[solve_boundary(j) for j in range(1,5)]
    nearest=min(bnds,key=lambda d:d['rho_found'])
    # Transport nearest A1/A2 candidate across rungs
    j=nearest['j']; z1=Z0[:,0,:,:]; zj=Z0[:,j,:,:]
    z1p=_inject(nearest['x'][:24],z1); zjp=_inject(nearest['x'][24:],zj)
    transport=[]
    for q in (1,2,3,4,7,10):
        raw01=z1**(1/q); rawp1=z1p**(1/q); raw02=zj**(1/q); rawp2=zjp**(1/q)
        raw_euc=np.linalg.norm(np.r_[(rawp1-raw01).ravel(),(rawp2-raw02).ravel()])
        dq=np.linalg.norm(np.r_[(rawp1**q-raw01**q).ravel(),(rawp2**q-raw02**q).ravel()])
        # score via canonical map is exactly same; use canonical score to avoid duplicate formula noise
        margin=canon_alt_score(z1p)-canon_alt_score(zjp)
        transport.append([q,raw_euc,dq,margin])
    # CSVs
    with open(OUT/'HE2019_REBUILT_Q_SWEEP.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['q','raw_rank','matched_rank','raw_A1','raw_A2','raw_A5','max_matched_score_deviation']); w.writerows(sweep)
    with open(OUT/'HE2019_RADIUS_APPLICATION.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['challenger','baseline_margin','grad_norm','rho_lin','rho_found_upper_bound','tie_residual','min_simplex_slack'])
        for d in bnds: w.writerow([f'A{d["j"]+1}',d['margin'],d['grad_norm'],d['rho_lin'],d['rho_found'],d['tie_residual'],d['simplex_slack']])
    with open(OUT/'HE2019_RADIUS_TRANSPORT_CHECK.csv','w',newline='') as f:
        w=csv.writer(f); w.writerow(['q','raw_coordinate_L2_to_same_boundary','canonical_pulledback_dq','tie_margin']); w.writerows(transport)
    summary={
      'table4_rounded_matches':int(np.sum(rounded==TABLE4)),
      'table4_total_components':60,
      'table4_single_mismatch':{'location':'A1-C3-neutral','reproduced':float(coll[0,2,1]),'printed':0.09,'difference':float(coll[0,2,1]-.09)},
      'q3_max_final_component_diff_vs_locked_6dp':float(np.max(np.abs(fin-LOCKED_FINAL))),
      'q3_max_score_diff_vs_locked':float(np.max(np.abs(sc-LOCKED_SCORES))),
      'raw_first_integer_transition':7,
      'continuous_A2_A5_crossing':float(root),
      'matched_max_score_deviation_q1_10':float(max(np.max(abs(s-base_match)) for s in matched_scores)),
      'nearest_numerically_found_challenger':f'A{nearest["j"]+1}',
      'nearest_rho_lin':float(nearest['rho_lin']),
      'nearest_rho_found_upper_bound':float(nearest['rho_found']),
      'radius_claim':'numerically found feasible boundary distance / upper bound only; not certified global optimum',
      'canonical_metric_transport_check_max_dev':float(max(abs(r[2]-nearest['rho_found']) for r in transport)),
    }
    (OUT/'HE2019_REBUILD_SUMMARY.json').write_text(json.dumps(summary,indent=2))
    print('HE2019 REBUILT EVALUATOR: PASS')
    print(json.dumps(summary,indent=2))

if __name__=='__main__': main()
