#!/usr/bin/env python3
"""Directed-MPFR complete-winner ball covering; all pruning uses validated bounds."""
from __future__ import annotations
import argparse,gzip,json,heapq,time,math,sys
from pathlib import Path
from fractions import Fraction
import numpy as np
import kernel as K
HERE=Path(__file__).resolve().parent
DATA_TEXT=(HERE/'mpfr_interval_kernel.cpp').read_text().split('static const int DATA[200]={')[1].split('}')[0]
DIGITS=np.array(list(map(int,DATA_TEXT.split(','))))
CENTRE=[Fraction(int(x)**3,1000) for x in DIGITS]
def outward(frac,direction):
 a=float(frac)
 if direction<0 and Fraction.from_float(a)>frac:a=np.nextafter(a,-np.inf)
 if direction>0 and Fraction.from_float(a)<frac:a=np.nextafter(a,np.inf)
 return float(a)
def pair_idx(j):return np.array([np.ravel_multi_index((e,a,c,k),(4,5,5,2)) for e in range(4) for a in [0,j] for c in range(5) for k in range(2)],dtype=np.int32)
def rootbox(idx,R):
 r=Fraction.from_float(R)
 return (np.array([outward(max(Fraction(0),CENTRE[i]-r),-1) for i in idx]),np.array([outward(min(Fraction(1),CENTRE[i]+r),1) for i in idx]))
def state_margin(idx,lo,hi,j):
 l=np.full(200,np.nan);h=l.copy();l[idx]=lo;h[idx]=hi
 return K.margin(j,l,h)
def certify(radius,max_evals,out):
 K.lib.kernel_precision(160);out.mkdir(exist_ok=True,parents=True);R=outward(Fraction(radius),1)
 other=[]
 for j in [1,3,4]:
  idx=pair_idx(j);lo,hi=rootbox(idx,R);g=state_margin(idx,lo,hi,j)
  other.append(dict(challenger=f'A{j+1}',margin_lower=float(g[0]),margin_upper=float(g[1]),pass_=bool(g[0]>0)))
 print('OTHER BOXES',other,flush=True)
 idx=pair_idx(2);lo,hi=rootbox(idx,R);x=np.array([float(CENTRE[i]) for i in idx]);grad=np.zeros(len(x));h=1e-7
 for i in range(len(x)):
  xp=x.copy();xm=x.copy();xp[i]+=h;xm[i]-=h
  grad[i]=(np.mean(state_margin(idx,xp,xp,2))-np.mean(state_margin(idx,xm,xm,2)))/(2*h)
 heap=[];counter=evals=leaves=pruned=0;smallest_leaf=math.inf;beg=time.time()
 with gzip.open(out/'COVER_TREE.jsonl.gz','wt',encoding='utf-8') as f:
  def emit(obj):f.write(json.dumps(obj,separators=(',',':'))+'\n')
  emit(dict(kind='header',radius_decimal=radius,radius_dyadic_hex=R.hex(),bits=160,challenger=2,idx=idx.tolist(),library=K.lib.kernel_version().decode()))
  def add(lo,hi,parent,side):
   nonlocal counter,evals,leaves,pruned,smallest_leaf
   node=counter;counter+=1;cc=K.contract(idx,lo,hi,R)
   if cc is None:
    pruned+=1;emit(dict(kind='outside',id=node,parent=parent,side=side));return
   lo,hi=cc;g=state_margin(idx,lo,hi,2);evals+=1
   if g[0]>0:
    leaves+=1;smallest_leaf=min(smallest_leaf,float(g[0]));emit(dict(kind='positive',id=node,parent=parent,side=side,margin_lower=float(g[0]),margin_upper=float(g[1])));return
   emit(dict(kind='pending',id=node,parent=parent,side=side,margin_lower=float(g[0])))
   heapq.heappush(heap,(float(g[0]),node,lo,hi))
  add(lo,hi,None,None)
  while heap and evals<max_evals:
   _,node,lo,hi=heapq.heappop(heap);w=hi-lo;k=int(np.argmax(w*(np.abs(grad)+.05*np.max(np.abs(grad)))))
   mid=float(lo[k]+(hi[k]-lo[k])/2)
   if not lo[k]<mid<hi[k]:raise RuntimeError('unsplittable unresolved box')
   emit(dict(kind='split',id=node,coordinate=k,midpoint_hex=mid.hex()))
   h1=hi.copy();h1[k]=mid;l2=lo.copy();l2[k]=mid;add(lo,h1,node,'L');add(l2,hi,node,'R')
   if evals%1000<2:print('evals',evals,'queue',len(heap),'elapsed',round(time.time()-beg,1),flush=True)
  ok=not heap and all(r['pass_'] for r in other)
  emit(dict(kind='end',pass_=ok,queue=len(heap),evals=evals,nodes=counter,positive_leaves=leaves,outside_leaves=pruned))
 summary=dict(status='PASS' if ok else 'UNRESOLVED_BUDGET',radius_decimal=radius,radius_dyadic_hex=R.hex(),bits=160,interval_evaluations=evals,nodes=counter,queue=len(heap),positive_leaves=leaves,outside_leaves=pruned,smallest_positive_leaf_lower=smallest_leaf,other_challengers=other,elapsed_seconds=time.time()-beg,trust_scope='MPFR-directed-rounding arithmetic plus audited interval transformations and covering algorithm; not a proof-assistant result',radius_upper_unchanged='0.003514796840392008',global_optimum='NOT_CLAIMED')
 (out/'RADIUS_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True);return ok
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--radius',default='0.000995');ap.add_argument('--max',type=int,default=30000);ap.add_argument('--out',type=Path,default=HERE/'radius_000995');a=ap.parse_args();sys.exit(0 if certify(a.radius,a.max,a.out) else 1)
