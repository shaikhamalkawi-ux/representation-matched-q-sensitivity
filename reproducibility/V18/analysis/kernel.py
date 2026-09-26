"""Portable build/ctypes loader; original mathematical wrappers unchanged.

Requires preinstalled C++17, MPFR and GMP. Set QROF_MPFR_TOOLCHAIN to an
optional compiler prefix containing bin; otherwise use CXX or g++ on PATH.
No dependencies are downloaded or installed by this module.
"""
from __future__ import annotations
import ctypes as C
import os
from pathlib import Path
import shutil
import subprocess
import numpy as np
HERE=Path(__file__).resolve().parent
SO=HERE/('mpfr_interval_kernel.dll' if os.name=='nt' else 'mpfr_interval_kernel.so')
prefix=os.environ.get('QROF_MPFR_TOOLCHAIN')
TOOLCHAIN=Path(prefix).resolve() if prefix else None
DLL_DIRECTORY=os.add_dll_directory(str(TOOLCHAIN/'bin')) if os.name=='nt' and TOOLCHAIN else None

def build():
    if not SO.exists() or SO.stat().st_mtime < (HERE/'mpfr_interval_kernel.cpp').stat().st_mtime:
        compiler=(str(TOOLCHAIN/'bin'/('g++.exe' if os.name=='nt' else 'g++'))
                  if TOOLCHAIN else os.environ.get('CXX') or shutil.which('g++'))
        if not compiler:
            raise RuntimeError('Preinstalled C++17 compiler required; set CXX or QROF_MPFR_TOOLCHAIN')
        env=os.environ.copy()
        if TOOLCHAIN:
            env['PATH']=str(TOOLCHAIN/'bin')+os.pathsep+env.get('PATH','')
        flags=['-static-libgcc','-static-libstdc++'] if os.name=='nt' else ['-fPIC']
        command=[compiler,'-O2','-std=c++17','-fno-fast-math','-shared',*flags,
                 str(HERE/'mpfr_interval_kernel.cpp'),'-lmpfr','-lgmp','-o',str(SO)]
        subprocess.run(command,check=True,env=env)
build()
lib=C.CDLL(str(SO));P=C.POINTER(C.c_double); PI=C.POINTER(C.c_int)
lib.kernel_scores.argtypes=[P,P,C.c_uint,C.c_uint,P];lib.kernel_scores.restype=C.c_int
lib.kernel_margin.argtypes=[P,P,C.c_int,P];lib.kernel_margin.restype=C.c_int
lib.kernel_contract.argtypes=[C.c_int,PI,P,P,C.c_double];lib.kernel_contract.restype=C.c_int
lib.kernel_op.argtypes=[C.c_int,C.c_double,C.c_double,C.c_double,C.c_double,P];lib.kernel_op.restype=C.c_int
lib.kernel_version.restype=C.c_char_p;lib.kernel_error.restype=C.c_char_p
lib.kernel_precision.argtypes=[C.c_uint]

def ptr(a): return a.ctypes.data_as(P)
def checked(rc):
    if rc: raise RuntimeError(lib.kernel_error().decode())
def bounds(lo,hi):
    lo=np.ascontiguousarray(lo,dtype=np.float64).ravel()
    hi=np.ascontiguousarray(hi,dtype=np.float64).ravel()
    if lo.size!=200 or hi.size!=200: raise ValueError('200 canonical endpoints required')
    if not np.array_equal(np.isnan(lo),np.isnan(hi)): raise ValueError('Both endpoints must use the same fixed-coordinate NaN sentinel')
    mask=~np.isnan(lo)
    if np.any(~np.isfinite(lo[mask])) or np.any(~np.isfinite(hi[mask])) or np.any(lo[mask]>hi[mask]) or np.any(lo[mask]<0) or np.any(hi[mask]>1):
        raise ValueError('Canonical intervals must lie in [0,1]')
    return lo,hi
def scores(qe=3,qc=3,lo=None,hi=None):
    out=np.zeros(10,dtype=np.float64)
    if lo is not None:
        lo,hi=bounds(lo,hi)
    elif hi is not None: raise ValueError('Both lower and upper arrays are required')
    if not isinstance(qe,int) or not isinstance(qc,int) or qe<1 or qc<1: raise ValueError('Positive integer rungs required')
    checked(lib.kernel_scores(ptr(lo) if lo is not None else None,ptr(hi) if hi is not None else None,qe,qc,ptr(out)))
    return out.reshape(5,2)
def margin(j,lo,hi):
    out=np.zeros(2);lo,hi=bounds(lo,hi)
    checked(lib.kernel_margin(ptr(lo),ptr(hi),j,ptr(out)));return out

def contract(idx,lo,hi,r):
    idx=np.ascontiguousarray(idx,dtype=np.int32);lo=np.array(lo,dtype=np.float64,copy=True);hi=np.array(hi,dtype=np.float64,copy=True)
    if not np.isfinite(r) or r<=0 or len(idx)!=len(lo) or len(idx)!=len(hi) or np.any(idx<0) or np.any(idx>=200) or np.any(~np.isfinite(lo)) or np.any(~np.isfinite(hi)) or np.any(lo>hi):
        raise ValueError('Invalid ball-contractor arguments')
    rc=lib.kernel_contract(len(idx),idx.ctypes.data_as(PI),ptr(lo),ptr(hi),r)
    if rc<0:checked(rc)
    return None if rc else (lo,hi)
def op(code,a,b=(0.,0.)):
    out=np.zeros(2);checked(lib.kernel_op(code,*map(float,a),*map(float,b),ptr(out)));return out
if __name__=='__main__':
    import time
    print('MPFR',lib.kernel_version().decode())
    t=time.time();print(scores());print('elapsed',time.time()-t)
