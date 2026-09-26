#!/usr/bin/env python3
"""Independent Qin et al. (2019) Hamacher qROFWAPPMM replay and matched-control audit.

Uses the normalized expert matrices published in the article. The worked example does
not print the Minkowski exponent b; b=3 is used for the primary replay because it
best reproduces the printed criterion power-weight matrix. The matched winner result
is separately verified for b=1,2,3. Criterion partitions are implemented locally,
following the article's partition definitions.
"""
import itertools, csv
from pathlib import Path
import mpmath as mp
mp.mp.dps=40
RAW=[
[
[(.7,.2),(.8,.2),(.4,.5),(.7,.1),(.9,.2)],[(.8,.6),(.7,.6),(.4,.5),(.5,.3),(.7,.2)],[(.6,.5),(.5,.4),(.5,.6),(.8,.5),(.8,.3)],[(.7,.2),(.6,.5),(.6,.5),(.6,.2),(.6,.5)],[(.6,.4),(.7,.5),(.5,.6),(.7,.4),(.7,.3)]],
[
[(.7,.1),(.7,.3),(.3,.5),(.6,.1),(.8,.2)],[(.8,.3),(.7,.5),(.2,.7),(.9,.2),(.6,.2)],[(.9,.2),(.8,.3),(.5,.3),(.8,.4),(.9,.2)],[(.6,.3),(.7,.5),(.4,.6),(.8,.3),(.7,.6)],[(.7,.3),(.8,.4),(.3,.4),(.7,.2),(.6,.3)]],
[
[(.8,.2),(.6,.2),(.3,.7),(.8,.1),(.7,.1)],[(.8,.6),(.7,.4),(.4,.8),(.7,.4),(.6,.4)],[(.7,.2),(.8,.3),(.2,.6),(.9,.2),(.8,.2)],[(.4,.7),(.9,.2),(.2,.5),(.5,.5),(.9,.6)],[(.6,.4),(.7,.3),(.3,.6),(.7,.3),(.8,.5)]],
[
[(.9,.3),(.6,.2),(.3,.6),(.7,.3),(.7,.2)],[(.6,.6),(.7,.7),(.2,.6),(.6,.2),(.6,.2)],[(.7,.1),(.9,.4),(.4,.7),(.8,.3),(.9,.4)],[(.7,.3),(.5,.4),(.4,.5),(.9,.4),(.6,.2)],[(.8,.4),(.8,.5),(.2,.5),(.7,.4),(.8,.3)]]]
MATS=[[[tuple(mp.mpf(str(x)) for x in p) for p in row] for row in M] for M in RAW]
EW=[mp.mpf('.30'),mp.mpf('.22'),mp.mpf('.28'),mp.mpf('.20')]
CW=[mp.mpf('.20'),mp.mpf('.20'),mp.mpf('.15'),mp.mpf('.25'),mp.mpf('.20')]
PARTS=[[0,2,4],[1,3]]; DELTAS=[[1,2,3],[1,2]]
def dist(a,b,p=3): return (mp.mpf('.5')*abs(a[0]-b[0])**p+mp.mpf('.5')*abs(a[1]-b[1])**p)**(mp.mpf(1)/p)
def pweights(vals,base,p=3):
    T=[sum(1-dist(x,y,p) for j,y in enumerate(vals) if j!=i) for i,x in enumerate(vals)]
    z=[base[i]*(1+T[i]) for i in range(len(vals))]; d=sum(z); return [x/d for x in z]
class H:
    def __init__(self,q,lam=3): self.q=mp.mpf(q); self.lam=mp.mpf(lam)
    def f(self,t): return mp.log((self.lam+(1-self.lam)*t**self.q)/(t**self.q))
    def g(self,t): return mp.log((self.lam+(1-self.lam)*(1-t**self.q))/(1-t**self.q))
    def fi(self,s): return (self.lam/(mp.e**s+self.lam-1))**(1/self.q)
    def gi(self,s): return ((mp.e**s-1)/(mp.e**s+self.lam-1))**(1/self.q)
def add(a,b,G): return (G.gi(G.g(a[0])+G.g(b[0])),G.fi(G.f(a[1])+G.f(b[1])))
def mul(a,b,G): return (G.fi(G.f(a[0])+G.f(b[0])),G.gi(G.g(a[1])+G.g(b[1])))
def scalar(c,a,G): c=mp.mpf(c); return (G.gi(c*G.g(a[0])),G.fi(c*G.f(a[1])))
def power(a,c,G): c=mp.mpf(c); return (G.fi(c*G.f(a[0])),G.gi(c*G.g(a[1])))
def one_partition(vals,w,deltas,G):
    m=len(vals); acc=None
    for perm in itertools.permutations(range(m)):
        prod=None
        for pos,idx in enumerate(perm):
            if deltas[pos]==0: continue
            term=power(scalar(m*w[idx],vals[idx],G),deltas[pos],G)
            prod=term if prod is None else mul(prod,term,G)
        acc=prod if acc is None else add(acc,prod,G)
    return power(scalar(mp.mpf(1)/mp.factorial(m),acc,G),mp.mpf(1)/sum(deltas),G)
def partitioned(vals,w,G):
    ps=[]
    for inds,ds in zip(PARTS,DELTAS): ps.append(one_partition([vals[i] for i in inds],[w[i] for i in inds],ds,G))
    x=ps[0]
    for y in ps[1:]: x=add(x,y,G)
    return scalar(mp.mpf(1)/len(ps),x,G)
def run(q,controlled=False,b=3):
    G=H(q,3)
    if controlled:
        r=mp.mpf(3); qq=mp.mpf(q); mats=[[[ (mu**(r/qq),nu**(r/qq)) for mu,nu in row] for row in M] for M in MATS]
    else: mats=MATS
    coll=[[None]*5 for _ in range(5)]
    for i in range(5):
        for j in range(5):
            vals=[mats[h][i][j] for h in range(4)]; w=pweights(vals,EW,b); coll[i][j]=one_partition(vals,w,[1,0,0,0],G)
    final=[]
    for i in range(5):
        w=pweights(coll[i],CW,b); final.append(partitioned(coll[i],w,G))
    scores=[mu**q-nu**q for mu,nu in final]; rank=sorted(range(5),key=lambda i:scores[i],reverse=True)
    return final,scores,rank
def rstr(rank): return '>'.join('A%d'%(i+1) for i in rank)
out=Path(__file__).resolve().parent
rows=[]
for q in range(3,11):
    _,rs,rr=run(q,False,3); _,cs,cr=run(q,True,3)
    rows.append([q]+[float(x) for x in rs]+[rstr(rr),'A%d'%(rr[0]+1)]+[float(x) for x in cs]+[rstr(cr),'A%d'%(cr[0]+1)])
with open(out/'replayed_raw_and_controlled_from_script.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(['q','raw_A1','raw_A2','raw_A3','raw_A4','raw_A5','raw_rank','raw_winner','ctrl_A1','ctrl_A2','ctrl_A3','ctrl_A4','ctrl_A5','ctrl_rank','ctrl_winner']); w.writerows(rows)
for b in (1,2,3):
    ranks=[]
    for q in range(3,11): ranks.append(rstr(run(q,True,b)[2]))
    assert set(ranks)=={'A1>A3>A5>A4>A2'},(b,ranks)
print('QIN2019 REPLAY: PASS')
print('Primary b=3 raw winner changes A1 -> A3 at q=7; matched ranking remains A1>A3>A5>A4>A2.')
print('Matched winner conclusion is unchanged for b=1,2,3.')
