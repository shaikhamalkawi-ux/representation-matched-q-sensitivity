#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction
import csv, math, re, sys
import numpy as np
import kernel as K
HERE=Path(__file__).resolve().parent
VEC=HERE / 'data' / 'QIN2019_TOP60_HP_TIE_VECTOR.csv'
K.lib.kernel_precision(256)
text=(HERE/'mpfr_interval_kernel.cpp').read_text()
digs=list(map(int,text.split('static const int DATA[200]={')[1].split('}')[0].split(',')))

def enclose_frac(x: Fraction):
    a=float(x)
    lo=a; hi=a
    fa=Fraction.from_float(a)
    if fa>x: lo=math.nextafter(a,-math.inf)
    if fa<x: hi=math.nextafter(a, math.inf)
    return (lo,hi)

base=[enclose_frac(Fraction(d**3,1000)) for d in digs]
w=[None]*200
for i,b in enumerate(base): w[i]=b
rows=list(csv.DictReader(VEC.open()))
changed=[]
for r in rows:
    idx=np.ravel_multi_index((int(r['e']),int(r['alt']),int(r['criterion']),int(r['coord'])),(4,5,5,2))
    v=float(r['tie_z'])
    w[idx]=(v,v)
    changed.append(idx)
lo=np.full(200,np.nan);hi=np.full(200,np.nan)
for idx in changed:
    lo[idx]=w[idx][0]; hi[idx]=w[idx][1]
margin=K.margin(2,lo,hi)
# Directed L2 distance from exact baseline to declared finite-precision witness.
s=(0.0,0.0)
for idx in changed:
    d=K.op(1,w[idx],base[idx])
    a=K.op(7,d)
    sq=K.op(4,a,(2.0,2.0))
    s=K.op(0,s,sq)
dist=K.op(5,s,(2.0,2.0))
# Canonical simplex check for all 100 pairs.
maxsum=(-math.inf,-math.inf); minslack=(math.inf,math.inf)
for p in range(100):
    sm=K.op(0,w[2*p],w[2*p+1])
    if sm[1]>maxsum[1]: maxsum=sm
    slack=K.op(1,(1.0,1.0),sm)
    if slack[0]<minslack[0]: minslack=slack
status = bool(margin[1] < 0 and dist[1] > 0 and maxsum[1] <= 1.0 and minslack[0] >= 0)
print('MPFR',K.lib.kernel_version().decode())
print('A1-A3 margin interval =',repr(float(margin[0])),repr(float(margin[1])))
print('L2 distance interval =',repr(float(dist[0])),repr(float(dist[1])))
print('max pair sum interval =',repr(float(maxsum[0])),repr(float(maxsum[1])))
print('minimum simplex slack interval =',repr(float(minslack[0])),repr(float(minslack[1])))
print('DIRECTED UPPER WITNESS:', 'PASS' if status else 'FAIL')
if not status: sys.exit(1)
