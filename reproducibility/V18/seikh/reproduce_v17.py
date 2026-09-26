"""Optimizer-free V17 clean-directory reproduction; no network or installation.

Run this script from the extracted V17 research folder. All scientific inputs
are resolved relative to this file. Python and a preinstalled C++/MPFR/GMP
toolchain are required; mpmath is required for the finite diagnostic tests.
Use --toolchain and --python-deps for explicit dependency locations if needed.
Existing output directories are never reused or overwritten. Every retained
certificate is regenerated, including every internal and terminal tree record.
"""
from __future__ import annotations
import argparse,copy,hashlib,json,os,platform,shutil,subprocess,sys,time
from fractions import Fraction as F
from pathlib import Path

HERE=Path(__file__).resolve().parent
INPUT_SHA='71351416bb73b22f3781bad962a023b6ea7362e8020429aab97b47d3f36d1797'
PINS={
 'vendor/model_inputs.json':INPUT_SHA,
 'vendor/exact_models.py':'8d5941c32cd7cf7140c3273ac6805a370f73bcea627d321cd2ba3c2bfc62f427',
 'vendor/dyadic_interval.py':'38cae911d1feb49562dcfd3693bb56a8e1bb8219c8afa51f89d44aff5998852e',
 'mpfr_check/model_inputs.json':INPUT_SHA,
 'mpfr_check/vendor/seikh_unbounded_kernel.cpp':'66f706f01537295a0edd8e7427b1225b08cf54cd476108fc0c9d61d6fa6a277c',
 'mpfr_check/vendor/inherited_mpfr_interval_kernel.cpp':'128559782aa752ef39e34d0d18914b1c1033711d6c726dca95847769530053d1',
 'mpfr_check/vendor/data_values.inc':'405df7e1c3ccc85ff65f4468072930d3e151a194ffbb36067847767fbab61696'}
EXPECTED={
 'Q4_H_0_019.json':{'radius':'19/1000','q':'4','nodes':1,'leaves':1},
 'Q8_H_0_0215.json':{'radius':'43/2000','q':'8','nodes':1,'leaves':1},
 'Q16_H_0_0235.json':{'radius':'47/2000','q':'16','nodes':1,'leaves':1},
 'UNIFORM_H_0_019.json':{'radius':'19/1000','q':None,'nodes':73,'leaves':37}}
CODE=[
 'seikh_certificates.py','test_entropy_range.py',
 'witness_search/verify_fixed_rung_witnesses.py',
 'mpfr_check/build.py','mpfr_check/verify.py','mpfr_check/replay_v17.py',
 'mpfr_check/test_verify.py','mpfr_check/seikh_v17_components.cpp']
TIMING={'elapsed_seconds','elapsed_s','seconds','runtime_seconds'}

def require(condition,message):
 if not condition:raise ValueError(message)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def read(path):return json.loads(path.read_text(encoding='utf-8'))
def write_new(path,value):
 with path.open('x',encoding='utf-8') as stream:stream.write(json.dumps(value,indent=2)+'\n')
def substantive(value):
 if isinstance(value,dict):return {k:substantive(v) for k,v in value.items() if k not in TIMING}
 if isinstance(value,list):return [substantive(v) for v in value]
 return value
def accept_replay(candidate,fresh):
 require(substantive(candidate)==substantive(fresh),'deterministic replay mismatch beyond timing')
def check_identity(record,spec):
 require(record['status']=='PASS_COMPUTATIONAL' and record['version']==17,'certificate acceptance/version')
 require(record['input_model_sha256']==INPUT_SHA,'certificate input pin')
 require(F(record['normalized_baseline_raw_halfwidth'])==F(spec['radius']),'certificate radius identity')
 require(record['q_exact']==spec['q'],'certificate rung identity')
 require(record['node_count']==spec['nodes'] and record['leaf_count']==spec['leaves'],'unexpected tree counts')
 require(len(record['nodes'])==spec['nodes'] and len(record['leaves'])==spec['leaves'],'unexpected record lengths')
 require(record['unresolved_count']==0 and not record['unresolved'],'unresolved certificate')
 require(record['max_nodes']==2047 and record['max_depth']==18,'unexpected deterministic search limits')

def replay_negative_tests(fresh):
 """Compare deliberate, sometimes internally consistent forgeries with replay.

 `fresh` was generated in the new directory immediately before these tests.
 The acceptance comparison checks all substantive metadata and endpoint data,
 not merely topology or positivity. No bad fixture is saved over a real result.
 """
 def resummarize(d):
  d['node_count']=len(d['nodes']);d['leaf_count']=len(d['leaves'])
  d['margin_floor']=str(min(F(x) for n in d['leaves'] for x in n['margins']))
 def score_shift(d,all_scores=False):
  target=d['leaves'][0]['id'];node=next(n for n in d['nodes'] if n['id']==target)
  for i in (range(4) if all_scores else [1]):
   node['scores'][i]=[str(F(x)+F(1,2**224)) for x in node['scores'][i]]
  sc=[[F(x) for x in pair] for pair in node['scores']]
  node['margins']=[str(sc[1][0]-sc[j][1]) for j in (0,2,3)]
  node['positive']=min(map(F,node['margins']))>0
  d['leaves']=[copy.deepcopy(node) if n['id']==target else n for n in d['leaves']]
  resummarize(d)
 def remove_leaf(d):
  target=d['leaves'].pop()['id'];d['nodes']=[n for n in d['nodes'] if n['id']!=target];resummarize(d)
 def remove_internal(d):
  target=next(n['id'] for n in d['nodes'] if not n['positive'])
  d['nodes']=[n for n in d['nodes'] if n['id']!=target];resummarize(d)
 def mutate_interval(d):
  target=d['leaves'][0]['id']
  for collection in ('nodes','leaves'):
   n=next(n for n in d[collection] if n['id']==target);n['t'][0]=str(F(n['t'][0])+F(1,2**30))
 def mutate_weight(d):
  target=d['leaves'][0]['id']
  for collection in ('nodes','leaves'):
   n=next(n for n in d[collection] if n['id']==target)
   n['weights'][0]=[str(F(x)+F(1,2**224)) for x in n['weights'][0]]
 def duplicate(d):
  d['nodes'].append(copy.deepcopy(d['nodes'][-1]));d['node_count']=len(d['nodes'])
 cases=[('self_consistent_winner_score_margin_floor',score_shift),
        ('self_consistent_all_scores_translation',lambda d:score_shift(d,True)),
        ('removed_leaf_with_updated_counts',remove_leaf),
        ('removed_internal_node_with_updated_count',remove_internal),
        ('mutated_leaf_interval_and_copy',mutate_interval),
        ('mutated_weight_and_copy',mutate_weight),('duplicate_node_with_updated_count',duplicate),
        ('altered_radius_label',lambda d:d.update(normalized_baseline_raw_halfwidth='18/1000'))]
 results=[]
 for label,mutate in cases:
  candidate=copy.deepcopy(fresh);mutate(candidate)
  try:accept_replay(candidate,fresh)
  except ValueError as exc:results.append({'case':label,'rejected':True,'reason':str(exc)})
  else:raise ValueError('corrupted record accepted: '+label)
 timing=copy.deepcopy(fresh);timing['elapsed_seconds']=123456789
 accept_replay(timing,fresh)
 return {'status':'PASS','rejected_count':len(results),'tests':results,
         'timing_only_change_accepted':True,
         'scope':'Equality against newly regenerated integer records, excluding only timing fields.'}

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--output',type=Path,default=HERE/'clean_replay_01')
 parser.add_argument('--toolchain',type=Path,help='Preinstalled MPFR/GMP compiler prefix (containing bin).')
 parser.add_argument('--python-deps',type=Path,help='Optional directory containing preinstalled mpmath.')
 args=parser.parse_args();output=args.output.resolve()
 require(not sys.flags.optimize,'Do not use Python -O for diagnostic tests')
 require(not HERE.is_relative_to(output),'Output must not be the source folder or its ancestor')
 require(not output.exists(),'Refusing to overwrite/reuse an existing replay directory')
 actual={p.name for p in (HERE/'certificates').glob('*.json')}
 require(actual==set(EXPECTED),'Missing or unexpected certificate records: '+repr(actual))
 for relative,expected in PINS.items():require(sha(HERE/relative)==expected,'Immutable pin mismatch: '+relative)
 witness_path=HERE/'witness_search/FIXED_RUNG_WITNESSES.json';witness=read(witness_path)
 require(len(witness['witnesses'])==3 and {w['q_exact'] for w in witness['witnesses']}=={'4','8','16'},'Witness q-set/count mismatch')
 for name,spec in EXPECTED.items():check_identity(read(HERE/'certificates'/name),spec)
 selected=list(PINS)+CODE+['witness_search/FIXED_RUNG_WITNESSES.json']+['certificates/'+name for name in EXPECTED]
 before={name:sha(HERE/name) for name in selected}
 output.mkdir(parents=True,exist_ok=False);(output/'logs').mkdir();(output/'certificates').mkdir();(output/'retained_certificates').mkdir()
 for relative in list(PINS)+CODE+['witness_search/FIXED_RUNG_WITNESSES.json']:
  target=output/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(HERE/relative,target)
 for name in EXPECTED:shutil.copy2(HERE/'certificates'/name,output/'retained_certificates'/name)
 # This copy records the exact launcher used, but is not invoked recursively.
 shutil.copy2(Path(__file__).resolve(),output/'reproduce_v17.py')
 env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1';env.pop('PYTHONOPTIMIZE',None)
 if args.python_deps:env['PYTHONPATH']=str(args.python_deps.resolve())+os.pathsep+env.get('PYTHONPATH','')
 if args.toolchain:
  args.toolchain=args.toolchain.resolve();env['QROF_MPFR_TOOLCHAIN']=str(args.toolchain)
  env['PATH']=str(args.toolchain/'bin')+os.pathsep+env.get('PATH','')
 start=time.perf_counter();steps=[]
 summary={'status':'RUNNING','source':str(HERE),'output':str(output),'steps':steps,
          'python':sys.version,'platform':platform.platform(),'immutable_pins':PINS,
          'source_hashes_before':before,'expected_certificates':EXPECTED,
          'warning':'Finite tests are diagnostics; no optimizer, network, package install, or formal proof assistant is used.'}
 def run(label,command):
  then=time.perf_counter();p=subprocess.run([str(x) for x in command],cwd=output,env=env,
      stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,encoding='utf-8',errors='replace',timeout=600)
  (output/'logs'/(label+'.log')).write_text(p.stdout,encoding='utf-8')
  steps.append({'step':label,'returncode':p.returncode,'elapsed_seconds':time.perf_counter()-then,'command':[str(x) for x in command]})
  print(label,'PASS' if p.returncode==0 else 'FAIL',flush=True)
  require(p.returncode==0,label+' failed; see its fresh log');return p.stdout
 def py(label,script,*options):return run(label,[sys.executable,'-B',output/script,*options])
 try:
  dep=run('dependency_versions',[sys.executable,'-B','-c','import json,sys,mpmath; print(json.dumps({"python":sys.version,"mpmath":mpmath.__version__,"mpmath_file":mpmath.__file__}))'])
  summary['dependencies']=json.loads(dep)
  parity={}
  for name,spec in EXPECTED.items():
   options=['--radius',spec['radius'],'--max-nodes','2047','--output',output/'certificates'/name]
   if spec['q'] is not None:options.extend(['--q',spec['q']])
   py('integer_'+Path(name).stem,'seikh_certificates.py',*options)
   fresh=read(output/'certificates'/name);check_identity(fresh,spec)
   accept_replay(read(output/'retained_certificates'/name),fresh);parity[name]=True
  summary['all_substantive_certificate_fields_match']=parity
  negative=replay_negative_tests(read(output/'certificates/UNIFORM_H_0_019.json'))
  write_new(output/'INTEGER_REPLAY_NEGATIVE_TESTS.json',negative);summary['integer_record_negative_tests']=negative['rejected_count']
  py('exact_point_witnesses','witness_search/verify_fixed_rung_witnesses.py','--cross-rung')
  witness_result=read(output/'witness_search/EXACT_WITNESS_REPLAY.json')
  require(witness_result['status']=='PASS' and len(witness_result['checked_witnesses'])==3 and
          {x['q'] for x in witness_result['checked_witnesses']}=={'4','8','16'},'Exact witness replay incomplete')
  summary['integer_witness_negative_tests']=witness_result['negative_test_count']
  py('integer_entropy_diagnostics','test_entropy_range.py','--output',output/'ENTROPY_RANGE_TESTS.json')
  entropy=read(output/'ENTROPY_RANGE_TESTS.json')
  require(entropy['status']=='PASS' and entropy['tests_run']==21 and entropy['skipped_count']==0 and entropy['mpmath_version'],'Integer diagnostics incomplete/skipped')
  summary['integer_entropy_tests']=entropy['tests_run']
  buildopts=['--vendor-only']
  if args.toolchain:buildopts+=['--toolchain',args.toolchain]
  py('compile_fresh_mpfr','mpfr_check/build.py',*buildopts)
  py('mpfr_diagnostics','mpfr_check/test_verify.py')
  test_report=read(output/'mpfr_check/TESTS.json');require(test_report['status']=='PASS','MPFR diagnostics failed')
  summary['mpfr_negative_tests']=next(x['count'] for x in test_report['checks'] if x['check']=='negative_record_tests')
  require(summary['mpfr_negative_tests']==20,'MPFR negative test count changed')
  py('mpfr_certificates_and_witnesses','mpfr_check/replay_v17.py')
  mpfr=read(output/'mpfr_check/V17_MPFR_RECHECK.json')
  require(mpfr['status']=='PASS' and mpfr['certificate_count']==4 and mpfr['witness_count']==3,'MPFR replay count/status')
  require({x['file'] for x in mpfr['certificates']}==set(EXPECTED),'MPFR missing certificate')
  require({x['q'] for x in mpfr['witnesses']}=={'4','8','16'},'MPFR missing witness')
  require(sum(x['leaf_count'] for x in mpfr['certificates'])==40,'MPFR expected 40 leaf checks')
  summary['mpfr_recheck']={'certificates':4,'leaves':40,'witnesses':3,'bits':mpfr['mpfr_bits'],'version':test_report['mpfr_version']}
  summary['compiler']=read(output/'mpfr_check/BUILD.json')['compiler']
  summary['status']='PASS'
 except Exception as exc:
  summary['status']='FAIL';summary['error']=repr(exc);raise
 finally:
  after={name:sha(HERE/name) for name in selected};summary['source_hashes_after']=after
  summary['original_inputs_and_sources_unchanged']=before==after
  copied=[name for name in list(PINS)+CODE+['witness_search/FIXED_RUNG_WITNESSES.json']]
  summary['copied_sources_match_original']=all(sha(output/name)==before[name] for name in copied)
  if not summary['original_inputs_and_sources_unchanged'] or not summary['copied_sources_match_original']:summary['status']='FAIL'
  summary['elapsed_seconds']=time.perf_counter()-start
  summary['output_hashes']={p.relative_to(output).as_posix():sha(p) for p in sorted(output.rglob('*')) if p.is_file()}
  write_new(output/'CLEAN_REPLAY_SUMMARY.json',summary)
 require(summary['status']=='PASS','Fresh replay did not pass immutability checks')
 print('FINAL PASS',output,flush=True)

if __name__=='__main__':main()
