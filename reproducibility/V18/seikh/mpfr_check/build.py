"""Compile V17 MPFR adapter without installs or modifying received material.

The three inherited files and input JSON are copied unchanged into this folder
after pinned SHA256 verification. The copied vendor files travel with the new
package, so --vendor-only can rebuild it without the historical workspace tree.
"""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,time
HERE=Path(__file__).resolve().parent
WORKSPACE=HERE.parents[2]
ORIGINAL=WORKSPACE/'inputs/QROF_PostV14_Audit_20260926/extracted/QROF_PostV14_IndependentArithmetic_Audit_20260925'
if os.name=='nt':ORIGINAL=Path('\\\\?\\'+str(ORIGINAL))
DEFAULT_TOOLCHAIN=Path(os.environ['QROF_MPFR_TOOLCHAIN']) if os.environ.get('QROF_MPFR_TOOLCHAIN') else None
PINS={
 'seikh_unbounded_kernel.cpp':'66f706f01537295a0edd8e7427b1225b08cf54cd476108fc0c9d61d6fa6a277c',
 'inherited_mpfr_interval_kernel.cpp':'128559782aa752ef39e34d0d18914b1c1033711d6c726dca95847769530053d1',
 'data_values.inc':'405df7e1c3ccc85ff65f4468072930d3e151a194ffbb36067847767fbab61696',
 'model_inputs.json':'71351416bb73b22f3781bad962a023b6ea7362e8020429aab97b47d3f36d1797'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build(toolchain=DEFAULT_TOOLCHAIN,vendor_only=False):
 if os.name=='nt' and toolchain is None:
  raise ValueError('Windows requires --toolchain PATH or QROF_MPFR_TOOLCHAIN; no machine-specific default is assumed')
 if toolchain is not None:toolchain=Path(toolchain).resolve()
 vendor=HERE/'vendor';vendor.mkdir(exist_ok=True)
 for name,expected in PINS.items():
  dest=HERE/name if name=='model_inputs.json' else vendor/name
  if not vendor_only:
   source=ORIGINAL/('independent' if name=='model_inputs.json' else 'inherited/QROF_PostV14_Contribution_Study')/name
   if sha(source)!=expected:raise ValueError('Original pin mismatch: '+name)
   shutil.copy2(source,dest)
  if sha(dest)!=expected:raise ValueError('Vendor pin mismatch: '+name)
 env=os.environ.copy()
 if toolchain is not None:env['PATH']=str(toolchain/'bin')+os.pathsep+env.get('PATH','')
 compiler=toolchain/'bin/g++.exe' if os.name=='nt' else Path('g++')
 library=HERE/('seikh_v17_components.dll' if os.name=='nt' else 'seikh_v17_components.so')
 command=[str(compiler),'-O2','-std=c++17','-fno-fast-math','-shared',
          '-static-libgcc','-static-libstdc++',str(HERE/'seikh_v17_components.cpp'),'-lmpfr','-lgmp','-o',str(library)]
 if os.name!='nt':command.insert(5,'-fPIC')
 start=time.perf_counter();run=subprocess.run(command,env=env,cwd=HERE,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 (HERE/'BUILD.log').write_text(run.stdout,encoding='utf-8')
 if run.returncode:raise RuntimeError('Build failed; see BUILD.log')
 version=subprocess.check_output([str(compiler),'--version'],env=env,text=True).splitlines()[0]
 report={'status':'PASS','compiler':version,'command':command,'toolchain':str(toolchain) if toolchain is not None else 'system PATH',
         'elapsed_seconds':time.perf_counter()-start,'vendor_pins':PINS,
         'adapter_sha256':sha(HERE/'seikh_v17_components.cpp'),'library_sha256':sha(library)}
 (HERE/'BUILD.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
 return report
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--toolchain',type=Path,default=DEFAULT_TOOLCHAIN);p.add_argument('--vendor-only',action='store_true');a=p.parse_args()
 print(json.dumps(build(a.toolchain,a.vendor_only),indent=2))
