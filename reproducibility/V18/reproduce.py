"""Run released scientific checks in a fresh working copy, without installation.

Python 3.10+, numpy, scipy, mpmath, C++17, MPFR and GMP must already exist.
The release files themselves are never modified; use a new --output directory.
The slow Qin q=3 radius-cover replay is explicit via --radius-replay.
"""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent

def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def hashes(root): return {p.relative_to(root).as_posix():digest(p) for p in sorted(root.rglob('*')) if p.is_file()}
def require(ok,message):
    if not ok: raise RuntimeError(message)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--toolchain',type=Path,help='Optional preinstalled compiler prefix containing bin')
    parser.add_argument('--python-deps',type=Path,help='Optional preinstalled mpmath dependency directory')
    parser.add_argument('--radius-replay',action='store_true',help='Also validate the slower 18,515-node Qin q=3 cover')
    args=parser.parse_args()
    require(not sys.flags.optimize,'Do not use Python -O; inherited checks include assertions')
    output=args.output.resolve()
    require(not output.exists(),'Refusing to reuse an existing output directory')
    require(not output.is_relative_to(HERE) and not HERE.is_relative_to(output),'Output must be outside the release and not its ancestor')
    before=hashes(HERE)
    output.mkdir(parents=True,exist_ok=False)
    work=output/'work'; shutil.copytree(HERE,work)
    logs=output/'logs'; logs.mkdir()
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1';env.pop('PYTHONOPTIMIZE',None)
    prefix=args.toolchain.resolve() if args.toolchain else None
    if prefix:
        env['QROF_MPFR_TOOLCHAIN']=str(prefix)
        env['PATH']=str(prefix/'bin')+os.pathsep+env.get('PATH','')
    if args.python_deps:
        env['PYTHONPATH']=str(args.python_deps.resolve())+os.pathsep+env.get('PYTHONPATH','')
    compiler=str(prefix/'bin'/('g++.exe' if os.name=='nt' else 'g++')) if prefix else env.get('CXX') or shutil.which('g++')
    require(bool(compiler),'Preinstalled C++17 compiler required')
    suffix='.dll' if os.name=='nt' else '.so'
    runtime=['--runtime-dir',str(prefix/'bin')] if prefix and os.name=='nt' else []
    steps=[];start=time.perf_counter()
    report={'status':'RUNNING','steps':steps,'scope':'Fresh V17 core, seven-case benchmark, Qin/He quick checks, Qin all-q and full-order retained covers, Zhang control.',
            'radius_replay_requested':args.radius_replay,'manifest_before':before,'historical_30_candidate_ledger_available':False}
    def run(label,cmd,cwd=work):
        then=time.perf_counter()
        p=subprocess.run([str(x) for x in cmd],cwd=cwd,env=env,capture_output=True,text=True,
                         encoding='utf-8',errors='replace',timeout=1800)
        (logs/(label+'.log')).write_text(p.stdout+p.stderr,encoding='utf-8')
        steps.append({'step':label,'returncode':p.returncode,'seconds':time.perf_counter()-then})
        print(label,'PASS' if p.returncode==0 else 'FAIL',flush=True)
        require(p.returncode==0,label+' failed; inspect fresh log')
    def py(label,relative,*opts):run(label,[sys.executable,'-B',work/relative,*opts])
    def compile_kernel(stem):
        base=work/'qin/unbounded';target=base/(stem+suffix)
        flags=['-static-libgcc','-static-libstdc++'] if os.name=='nt' else ['-fPIC']
        run('compile_'+stem,[compiler,'-O2','-std=c++17','-fno-fast-math','-shared',*flags,
                            base/(stem+'.cpp'),'-lmpfr','-lgmp','-o',target],base)
        return target
    try:
        run('dependencies',[sys.executable,'-B','-c','import numpy,scipy,mpmath; print(numpy.__version__,scipy.__version__,mpmath.__version__)'])
        opts=['--output',output/'seikh_fresh']
        if prefix:opts+=['--toolchain',prefix]
        if args.python_deps:opts+=['--python-deps',args.python_deps.resolve()]
        py('seikh_full_integer_mpfr','seikh/reproduce_v17.py',*opts)
        py('benchmark_seven_cases','benchmark/run_all.py')
        py('qin_he_quick','analysis/run_public_analysis.py')
        if args.radius_replay:py('qin_q3_retained_radius_cover','analysis/replay_certificate.py')
        compact=compile_kernel('compact_q_kernel')
        continuous=compile_kernel('continuous_q_kernel')
        for model,merge in [('H','0'),('A','1')]:
            folder=work/'qin/unbounded'
            py('qin_unbounded_'+model,'qin/unbounded/replay_unbounded_q.py',folder/(model+'_UNBOUNDED_R12.json'),
               '--library',compact,*runtime,'--expected-merge',merge,'--expected-radius','1/4096',
               '--allow-rebuilt-library','--negative-controls','--output',output/(model+'_UNBOUNDED_REPLAY.json'))
            py('qin_fullorder_'+model,'qin/fullorder/faramondi_replay_allrank.py',
               work/'qin/fullorder'/('faramondi_'+model+'_allrank_R16.json'),
               '--library',continuous,'--kernel-sources',folder,*runtime,'--expected-model',model,
               '--expected-radius','1/65536','--expected-order','A1,A3,A5,A4,A2',
               '--validation-mode','rebuilt-binary','--malformed-tests','--output',output/(model+'_FULLORDER_REPLAY.json'))
        py('qin_source_crosscheck','qin/unbounded/source_crosscheck_compact.py','--library',compact,*runtime,
           '--output',output/'QIN_SOURCE_CROSSCHECK.json')
        py('zhang_control','zhang/zhang_pipeline_audit.py','--output-dir',output/'zhang_fresh')
        report['status']='PASS'
    except Exception as exc:
        report['status']='FAIL';report['error']=str(exc)
        raise
    finally:
        after=hashes(HERE)
        report['release_sources_unchanged']=before==after
        if before!=after:report['status']='FAIL'
        report['seconds']=time.perf_counter()-start
        (output/'REPLAY_SUMMARY.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    require(report['status']=='PASS','Replay failed')
    print('FINAL PASS',flush=True)

if __name__=='__main__': main()
