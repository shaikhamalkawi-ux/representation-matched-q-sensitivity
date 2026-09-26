"""Validated compact-parameter cover for all finite real q>=3, plus t=0 limit."""
import argparse
import ctypes as C
from fractions import Fraction as F
import hashlib
import json
import math
import os
from pathlib import Path
import re
import time

HERE=Path(__file__).resolve().parent
P=C.POINTER(C.c_double)

def outward(x, upper):
    value=float(x)
    if (upper and F.from_float(value)<x) or (not upper and F.from_float(value)>x):
        value=math.nextafter(value, math.inf if upper else -math.inf)
    return value

def floor12(x):
    n=x.numerator*10**12//x.denominator
    return ('-' if n<0 else '')+f'{abs(n)//10**12}.{abs(n)%10**12:012}'

class Kernel:
    def __init__(self, path, bits, runtime):
        self.handles=[]
        if os.name=='nt':
            for d in (path.parent, *runtime):
                self.handles.append(os.add_dll_directory(str(Path(d).resolve())))
        self.lib=C.CDLL(str(path.resolve()))
        self.lib.kernel_precision.argtypes=[C.c_uint]
        self.lib.kernel_precision(bits)
        self.lib.kernel_error.restype=C.c_char_p
        self.lib.kernel_version.restype=C.c_char_p
        self.lib.kernel_scores_compact.argtypes=[P,P,C.c_double,C.c_double,C.c_int,P]
        self.lib.kernel_scores_compact.restype=C.c_int
        self.lib.kernel_scores_limit.argtypes=[P,P,C.c_int,P]
        self.lib.kernel_scores_limit.restype=C.c_int
        self.lib.kernel_scores_real.argtypes=[P,P,C.c_double,C.c_double,C.c_double,C.c_double,C.c_int,P]
        self.lib.kernel_scores_real.restype=C.c_int

    def evaluate(self,tlo,thi,lo,hi,merge,mode='compact'):
        a=None if lo is None else (C.c_double*200)(*lo)
        b=None if hi is None else (C.c_double*200)(*hi)
        out=(C.c_double*10)()
        if mode=='limit':
            rc=self.lib.kernel_scores_limit(a,b,merge,out)
        elif mode=='real':
            rc=self.lib.kernel_scores_real(a,b,float(tlo),float(thi),float(tlo),float(thi),merge,out)
        else:
            rc=self.lib.kernel_scores_compact(a,b,outward(tlo,False),outward(thi,True),merge,out)
        if rc: raise RuntimeError(self.lib.kernel_error().decode())
        result=[list(out[2*i:2*i+2]) for i in range(5)]
        assert all(math.isfinite(x) and math.isfinite(y) and x<=y for x,y in result)
        return result

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--library',type=Path,required=True)
    p.add_argument('--runtime-dir',action='append',default=[])
    p.add_argument('--bits',type=int,default=192)
    p.add_argument('--merge',type=int,choices=(0,1),required=True)
    p.add_argument('--radius-power',type=int,default=12)
    p.add_argument('--zero-radius',action='store_true')
    p.add_argument('--max-depth',type=int,default=16)
    p.add_argument('--max-nodes',type=int,default=16383)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    assert 64<=a.bits<=4096 and a.radius_power>0
    kernel=Kernel(a.library,a.bits,a.runtime_dir)
    cpp=(HERE/'mpfr_interval_kernel.cpp').read_text()
    data=[int(x) for x in re.search(r'DATA\[200\]=\{([^}]+)\}',cpp).group(1).split(',')]
    assert len(data)==200
    centre=[F(x**3,1000) for x in data]
    radius=F(0) if a.zero_radius else F(1,2**a.radius_power)
    lo=[outward(c-radius,False) for c in centre] if radius else None
    hi=[outward(c+radius,True) for c in centre] if radius else None
    if radius:
        assert min(lo)>0 and max(hi)<1
        assert all(F.from_float(hi[i])+F.from_float(hi[i+1])<1 for i in range(0,200,2))
    crosschecks=[]
    for q in (3,4,7,10,50,100,1000,1000000):
        old=kernel.evaluate(F(q),F(q),None,None,a.merge,'real')
        new=kernel.evaluate(F(1,q),F(1,q),None,None,a.merge)
        assert all(max(x[0],y[0])<=min(x[1],y[1]) for x,y in zip(old,new))
        crosschecks.append({'q':q,'real_q_intervals':old,'compact_intervals':new})
    at_zero=kernel.evaluate(F(0),F(0),None,None,a.merge)
    limit=kernel.evaluate(F(0),F(0),None,None,a.merge,'limit')
    assert all(max(x[0],y[0])<=min(x[1],y[1]) for x,y in zip(at_zero,limit))
    started=time.perf_counter()
    pending=[(F(0),F(1,3),0,'r')]
    nodes=[]; leaves=[]; unresolved=[]
    while pending:
        tl,th,depth,node_id=pending.pop()
        row={'id':node_id,'t_lower_exact':str(tl),'t_upper_exact':str(th),'depth':depth}
        try:
            scores=kernel.evaluate(tl,th,lo,hi,a.merge)
            margins=[F.from_float(scores[0][0])-F.from_float(scores[j][1]) for j in range(1,5)]
            positive=min(margins)>0
            row.update(score_intervals=scores,margin_lowers_exact=list(map(str,margins)),positive=positive)
        except RuntimeError as e:
            positive=False
            row.update(error=str(e),positive=False)
        nodes.append(row)
        if positive: leaves.append(row)
        elif depth>=a.max_depth or len(nodes)+len(pending)+2>a.max_nodes: unresolved.append(row)
        else:
            mid=(tl+th)/2
            pending.extend([(mid,th,depth+1,node_id+'1'),(tl,mid,depth+1,node_id+'0')])
        if len(nodes)%64==0:
            print(f'nodes={len(nodes)} positive={len(leaves)} unresolved={len(unresolved)} pending={len(pending)}',flush=True)
    leaves.sort(key=lambda r:F(r['t_lower_exact']))
    complete=bool(leaves) and not unresolved
    if complete:
        assert F(leaves[0]['t_lower_exact'])==0 and F(leaves[-1]['t_upper_exact'])==F(1,3)
        assert all(x['t_upper_exact']==y['t_lower_exact'] for x,y in zip(leaves,leaves[1:]))
        assert len(nodes)==2*len(leaves)-1
    minimum=min(F(v) for r in leaves for v in r['margin_lowers_exact']) if leaves else None
    result={'status':'PASS' if complete else 'UNRESOLVED_NOT_CERTIFIED','bits':a.bits,
      'mpfr_version':kernel.lib.kernel_version().decode(),'scope':'all finite real q>=3 via t=1/q in (0,1/3], plus positive-input continuous limit t=0',
      'merge':a.merge,'model':'H' if a.merge==0 else 'A','q_min_exact':'3','t_domain_exact':['0','1/3'],
      'canonical_radius_exact':str(radius),'input_intervals':{'lower':lo,'upper':hi},
      'input_set':'exact centre' if not radius else 'entire 200-coordinate canonical L-infinity cube, strictly inside positive simplexes',
      'winner':'A1' if complete else None,'all_four_challengers':True,'cover_complete':complete,
      'nodes':len(nodes),'positive_leaves':len(leaves),'unresolved_leaves':len(unresolved),
      'uniform_margin_lower_exact':str(minimum),'uniform_margin_floor_12dp':floor12(minimum) if minimum is not None else None,
      'elapsed_seconds':time.perf_counter()-started,'finite_q_crosschecks':crosschecks,
      'limit_scores':limit,'t_zero_scores':at_zero,
      'source_hashes':{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
        ('mpfr_interval_kernel.cpp','continuous_q_kernel.cpp','compact_q_kernel.cpp',Path(__file__).name)},
      'library_sha256':hashlib.sha256(a.library.read_bytes()).hexdigest(),
      'trust_boundary':['positive canonical inputs and collective grades for continuity at t=0',
         'declared H/A local-partition-size source reconstructions; fixed non-input parameters',
         'MPFR/GMP and audited interval primitive inclusion and native build',
         'closed t cover with no unresolved leaves, not a sampled infinite-q inference'],
      'leaf_cover':leaves,'all_nodes':nodes,'unresolved':unresolved}
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('leaf_cover','all_nodes','unresolved','input_intervals','finite_q_crosschecks')},indent=2))
    if not complete: raise SystemExit(2)

if __name__=='__main__': main()
