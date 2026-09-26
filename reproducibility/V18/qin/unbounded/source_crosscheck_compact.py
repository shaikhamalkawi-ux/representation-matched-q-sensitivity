"""Independent high-precision source-value containment checks for compact kernel.

Self-contained formulas in source_qin_models.py; requires mpmath and the native
compact library plus MPFR/GMP runtime. Does not import either interval generator
or replay checker. Point containment, not merely overlap, is explicitly tested.
This nondirected source diagnostic complements the separate interval cover proof.
"""
from __future__ import annotations
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

import mpmath as mp
import source_qin_models as model

HERE=Path(__file__).resolve().parent
P=C.POINTER(C.c_double)


def require(ok,message):
    if not ok: raise ValueError(message)


def mp_fraction(value):
    return mp.mpf(value.numerator)/value.denominator


def outward(value,up):
    f=float(value)
    if (up and F.from_float(f)<value) or (not up and F.from_float(f)>value):
        f=math.nextafter(f,math.inf if up else -math.inf)
    require(F.from_float(f)>=value if up else F.from_float(f)<=value,"outward conversion failed")
    return f


class Native:
    def __init__(self,path,runtimes,bits):
        self.handles=[]
        if os.name=='nt':
            for directory in (path.parent,*runtimes):
                self.handles.append(os.add_dll_directory(str(Path(directory).resolve())))
        self.lib=C.CDLL(str(path.resolve()))
        self.lib.kernel_precision.argtypes=[C.c_uint]
        self.lib.kernel_precision(bits)
        self.lib.kernel_error.restype=C.c_char_p
        self.lib.kernel_version.restype=C.c_char_p
        self.lib.kernel_scores_compact.argtypes=[P,P,C.c_double,C.c_double,C.c_int,P]
        self.lib.kernel_scores_compact.restype=C.c_int
        self.lib.kernel_scores_limit.argtypes=[P,P,C.c_int,P]
        self.lib.kernel_scores_limit.restype=C.c_int

    def evaluate(self,point,tlo,thi,merge,limit=False):
        p=None if point is None else (C.c_double*200)(*point)
        out=(C.c_double*10)()
        if limit:
            rc=self.lib.kernel_scores_limit(p,p,merge,out)
        else:
            rc=self.lib.kernel_scores_compact(p,p,outward(tlo,False),outward(thi,True),merge,out)
        require(rc==0,self.lib.kernel_error().decode() if rc else "")
        bounds=[list(out[2*i:2*i+2]) for i in range(5)]
        require(all(math.isfinite(lo) and math.isfinite(hi) and lo<=hi for lo,hi in bounds),"bad bounds")
        return bounds


def contains(value,bounds):
    # Binary64 export endpoints are reconstructed EXACTLY as rationals in mp.
    return mp_fraction(F.from_float(bounds[0]))<=value<=mp_fraction(F.from_float(bounds[1]))


def calculate(point_exact,q_exact,dps,partition,fixed=False):
    with mp.workdps(dps):
        return model.scores([mp_fraction(x) for x in point_exact],mp_fraction(q_exact),partition,fixed)


def maximum_difference(x,y):
    return max(abs(a-b) for xr,yr in zip(x,y) for a,b in zip(xr,yr))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--library',type=Path,required=True)
    parser.add_argument('--runtime-dir',type=Path,action='append',default=[])
    parser.add_argument('--bits',type=int,default=256)
    parser.add_argument('--output',type=Path,default=HERE/'source_compact_crosscheck.json')
    args=parser.parse_args()
    started=time.perf_counter()
    mp.mp.dps=160
    data=[int(x) for expert in model.RAW_TENTHS for row in expert for x in row.split()]
    require(len(data)==200,'source record count')
    centre=[F(x**3,1000) for x in data]
    cpp=(HERE/'mpfr_interval_kernel.cpp').read_text(encoding='utf-8')
    cppdata=[int(x) for x in re.search(r'DATA\[200\]\s*=\s*\{([^}]+)\}',cpp).group(1).split(',')]
    require(data==cppdata,'independent source/C++ DATA mismatch')
    perturbation=[float(x+F((-1)**i,8192)) for i,x in enumerate(centre)]
    perturbed_exact=[F.from_float(x) for x in perturbation]
    require(all(abs(x-c)<F(1,4096) for x,c in zip(perturbed_exact,centre)),'perturbation leaves cube')
    require(all(x>0 for x in perturbed_exact),'positive input restriction')
    require(all(perturbed_exact[i]+perturbed_exact[i+1]<1 for i in range(0,200,2)),'simplex restriction')
    native=Native(args.library,args.runtime_dir,args.bits)
    rows=[]
    comparisons=0
    for label,point_exact,native_point in [('rational_centre',centre,None),('perturbed_binary64',perturbed_exact,perturbation)]:
        for q in [F(25,8),F(10),F(50),F(1000),F(1000000)]:
            direct=calculate(point_exact,q,160,model.java_direct_one_partition)
            logarithmic=calculate(point_exact,q,160,model.one_partition)
            lower_precision=calculate(point_exact,q,110,model.java_direct_one_partition)
            dlog=maximum_difference(direct,logarithmic)
            dprecision=maximum_difference(direct,lower_precision)
            require(dlog<mp.mpf('1e-140'),'direct/logarithmic mismatch')
            require(dprecision<mp.mpf('1e-95'),'110/160-digit convergence failure')
            t=F(1)/q
            for merge in [0,1]:
                point_bounds=native.evaluate(native_point,t,t,merge)
                # Independent nonzero interval containing the exact reciprocal.
                width=t/1024
                range_bounds=native.evaluate(native_point,t-width,t+width,merge)
                for i,reference in enumerate(direct[merge]):
                    require(contains(reference,point_bounds[i]),f'source value outside point enclosure {label} q={q} model={merge} alt={i+1}')
                    require(contains(reference,range_bounds[i]),f'source value outside t-range enclosure {label} q={q} model={merge} alt={i+1}')
                    comparisons+=2
                rows.append(dict(input=label,q_exact=str(q),t_exact=str(t),model='H' if merge==0 else 'A',
                                 source_scores_100_digits=[mp.nstr(x,100) for x in direct[merge]],
                                 source_direct_log_max_abs_error=mp.nstr(dlog,40),
                                 source_110_160_digit_max_abs_error=mp.nstr(dprecision,40),
                                 point_t_binary64_enclosure=[outward(t,False),outward(t,True)],
                                 source_in_point_intervals=True,point_score_intervals=point_bounds,
                                 nonzero_t_interval_exact=[str(t-width),str(t+width)],
                                 source_in_range_intervals=True,range_score_intervals=range_bounds))
            print(f'Source containment passed: {label}, q={q}',flush=True)
        # At t=0, all supports converge to 1 and power weights to base weights.
        # Evaluate the source aggregation with those weights at operator q=1
        # and q=7; q-natural aggregation makes the final powered scores equal.
        direct=calculate(point_exact,F(1),160,model.java_direct_one_partition,True)
        logarithmic=calculate(point_exact,F(7),160,model.one_partition,True)
        lower_precision=calculate(point_exact,F(1),110,model.java_direct_one_partition,True)
        dlog=maximum_difference(direct,logarithmic)
        dprecision=maximum_difference(direct,lower_precision)
        require(dlog<mp.mpf('1e-140'),'base-weight limit source forms/rungs disagree')
        require(dprecision<mp.mpf('1e-95'),'limit precision check failed')
        for merge in [0,1]:
            compact_bounds=native.evaluate(native_point,F(0),F(0),merge)
            limit_bounds=native.evaluate(native_point,F(0),F(0),merge,True)
            for i,reference in enumerate(direct[merge]):
                require(contains(reference,compact_bounds[i]),'source limit outside compact t=0 bounds')
                require(contains(reference,limit_bounds[i]),'source limit outside base-weight shortcut bounds')
                comparisons+=2
            rows.append(dict(input=label,q_exact='infinite-limit',t_exact='0',model='H' if merge==0 else 'A',
                             source_operator_rungs_exact=['1','7'],source_base_weights_fixed=True,
                             source_scores_100_digits=[mp.nstr(x,100) for x in direct[merge]],
                             source_direct_log_max_abs_error=mp.nstr(dlog,40),
                             source_110_160_digit_max_abs_error=mp.nstr(dprecision,40),
                             source_in_compact_zero_intervals=True,compact_zero_intervals=compact_bounds,
                             source_in_limit_intervals=True,limit_intervals=limit_bounds))
        print(f'Source containment passed: {label}, t=0 base-weight limit',flush=True)
    files=[Path(__file__),HERE/'source_qin_models.py',HERE/'compact_q_kernel.cpp',HERE/'continuous_q_kernel.cpp',HERE/'mpfr_interval_kernel.cpp']
    result=dict(status='PASS',source_decimal_precisions=[110,160],interval_precision_bits=args.bits,
                mpfr_version=native.lib.kernel_version().decode(),mpmath_version=mp.__version__,
                finite_q_exact=['25/8','10','50','1000','1000000'],limit_t_exact='0',
                models=['H','A'],input_cases=['rational_centre','perturbed_binary64'],
                independent_source_scores=120,source_containment_assertions=comparisons,
                independent_source_DATA_matches_cpp=True,
                perturbed_point_binary64=perturbation,
                scope='nondirected independent source-value diagnostics; continuum/cube guarantee supplied by separate directed interval cover',
                endpoint_handling='exact rational q and 1/q; outward binary64 t bounds; exact binary64 source-coordinate and interval-endpoint reconstruction',
                source_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                native_library_sha256=hashlib.sha256(args.library.read_bytes()).hexdigest(),
                elapsed_seconds=time.perf_counter()-started,cases=rows)
    require(comparisons==240,'comparison count mismatch')
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ['cases','perturbed_point_binary64']},indent=2))


if __name__=='__main__': main()
