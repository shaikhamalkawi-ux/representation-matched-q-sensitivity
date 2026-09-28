"""Portable, bounded replay of V22 scoped numerical evidence and certificates.

No downloads, dependency installs, optimization, theorem proving or publication.
Original scripts are used unchanged; outputs always go to a new external folder.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
if os.name == 'nt' and not str(HERE).startswith('\\\\?\\'):
    HERE = Path('\\\\?\\'+str(HERE))
ORIGINAL = HERE/'original'
RADIUS = ORIGINAL/'outputs/FSS_STRENGTHENING_20260926/radius_refinement'
FLOOR = ORIGINAL/'outputs/FSS_EIGHT_HOUR_RESEARCH_20260926'
WITNESS = ORIGINAL/'outputs/FSS_V18_RELEASE_20260926/public_artifact/seikh/witness_search'


def pin(path):
    return {'bytes':path.stat().st_size, 'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}


def check_inputs():
    data = json.loads((HERE/'THEORY_MANIFEST.json').read_text(encoding='utf-8'))
    for name, expected in data['payloads'].items():
        if pin(HERE/name) != expected:
            raise ValueError('Changed package input: ' + name)
    return data['payloads']


def commands(output, profile):
    def py(path, *args, optimize=False):
        return [sys.executable, '-B'] + (['-OO'] if optimize else []) + [str(path), *map(str,args)]
    compact = FLOOR/'novelty_challenge/source_floor/compact'
    decimal = FLOOR/'entropy_family/source_floor/floor16_independent/check_floor16_decimal.py'
    jobs = {
        'floor16_mpfr': py(compact/'check_floor16_mpfr.py','verify','--source',compact/'FLOOR16_SOURCE_AND_INTEGER_NODES_R1.json','--output',output/'floor16_mpfr.json'),
        'floor16_decimal': py(decimal,'--output',output/'floor16_decimal.json'),
        'seikh_fixed_integer': py(RADIUS/'independent_audit/audit_fixed_q.py','--output',output/'seikh_fixed_integer.json'),
        'seikh_witnesses': py(WITNESS/'verify_fixed_rung_witnesses.py','--record',WITNESS/'FIXED_RUNG_WITNESSES.json','--output',output/'seikh_witnesses.json'),
    }
    for q, radius in [(4,'0_0232195'),(8,'0_0250515'),(16,'0_0258973')]:
        jobs[f'seikh_q{q}_mpfr'] = py(HERE/'mpfr_output_adapter.py','--record',RADIUS/f'Q{q}_H_{radius}.json','--output',output/f'seikh_q{q}_mpfr.json')
    if profile == 'full':
        jobs['floor16_decimal_optimized'] = py(decimal,'--output',output/'floor16_decimal_optimized.json',optimize=True)
        jobs['seikh_uniform_integer'] = py(RADIUS/'independent_audit/audit_uniform_joint.py','--input','UNIFORM_JOINT_H_0_0232195.json','--output',output/'seikh_uniform_integer.json')
        jobs['seikh_uniform_mpfr'] = py(HERE/'mpfr_output_adapter.py','--record',RADIUS/'UNIFORM_JOINT_H_0_0232195.json','--output',output/'seikh_uniform_mpfr.json')
        jobs['qol_decimal_cover'] = py(HERE/'qol_cover_adapter.py','--arithmetic','decimal','--output',output/'qol_decimal_cover.json')
        jobs['qol_mpfr_cover'] = py(HERE/'qol_cover_adapter.py','--arithmetic','mpfr','--output',output/'qol_mpfr_cover.json')
        fragile=FLOOR/'external_case/compact_multiwinner'
        jobs['s4_compact_mpfr']=py(fragile/'certify_compact_multiwinner.py','verify','--input',fragile/'COMPACT_SOURCE_AND_NODES.json','--output',output/'s4_compact_mpfr.json')
        jobs['s4_compact_decimal']=py(FLOOR/'root_capacity/replay_compact_multiwinner_decimal.py','--source',fragile/'COMPACT_SOURCE_AND_NODES.json','--output',output/'s4_compact_decimal.json')
        jobs['s4_exact_domain']=py(fragile/'check_compact_domain_exact.py','--input',fragile/'COMPACT_SOURCE_AND_NODES.json','--output',output/'s4_exact_domain.json')
        jobs['s4_rounding']=py(FLOOR/'entropy_family/compact_precision/check_decimal_rounding.py','--output',output/'s4_rounding.json')
        jobs['s6_cached_diagnostic']=py(HERE/'diagnostic_output_adapter.py','--output',output/'s6_cached_diagnostic.json')
        jobs['source2022']=py(HERE/'source2022_output_adapter.py','--output',output/'source2022.json')
    return jobs


def run_job(name, command, output, timeout):
    start = time.perf_counter()
    result = {'job':name, 'timeout_seconds':timeout, 'status':'STARTED'}
    # Public receipt deliberately records script basename/options, not host paths.
    result['script'] = next(Path(x).name for x in command if x.endswith('.py'))
    try:
        p = subprocess.run(command, cwd=output, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           encoding='utf-8', errors='replace', timeout=timeout)
        result.update(returncode=p.returncode, status='PASS_EXIT_ZERO' if p.returncode == 0 else 'FAIL_EXIT')
        stdout, stderr = p.stdout, p.stderr
    except subprocess.TimeoutExpired as exc:
        result.update(status='TIMEOUT', returncode=None)
        stdout, stderr = exc.stdout or '', exc.stderr or ''
        if isinstance(stdout, bytes): stdout = stdout.decode('utf-8', 'replace')
        if isinstance(stderr, bytes): stderr = stderr.decode('utf-8', 'replace')
    except OSError as exc:
        result.update(status='EXECUTION_FAILURE', returncode=None)
        stdout, stderr = '', str(exc)
    for suffix, value in [('stdout.log',stdout),('stderr.log',stderr)]:
        with (output/f'{name}.{suffix}').open('x',encoding='utf-8') as stream:
            stream.write(value)
    result['elapsed_seconds'] = time.perf_counter()-start
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profile', choices=['quick','full'], default='quick')
    p.add_argument('--out', type=Path, required=True, help='New directory outside this package')
    p.add_argument('--workers', type=int, choices=[1,2], default=2)
    p.add_argument('--job-timeout', type=int, default=600, help='Seconds per job; 1..600')
    args = p.parse_args()
    if sys.flags.optimize:
        raise ValueError('Use normal Python for the runner; only its explicit floor16 -OO control is optimized')
    if not 1 <= args.job_timeout <= 600:
        raise ValueError('Bound job timeout to 1..600 seconds')
    output = args.out.resolve()
    if os.name == 'nt' and not str(output).startswith('\\\\?\\'):
        output = Path('\\\\?\\'+str(output))
    if output.is_relative_to(HERE) or HERE.is_relative_to(output):
        raise ValueError('Use a new separate output directory, not inside/above this package')
    output.mkdir(parents=True,exist_ok=False)
    start=time.perf_counter()
    receipt={'schema':'V22_THEORY_FRESH_REPLAY_R1','status':'STARTED','profile':args.profile,
             'started_utc':datetime.now(timezone.utc).isoformat(),'workers':args.workers,
             'scope':'Scoped finite-node/full-source-box, Seikh/QOL certificates and ordinary S6/2022 diagnostics; not proof of the analytic capacity/resource/Frank claims',
             'empirical_module_replayed':False,'original_source_modified':False}
    try:
        before=check_inputs()
        jobs=commands(output,args.profile)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures={n:pool.submit(run_job,n,c,output,args.job_timeout) for n,c in jobs.items()}
            receipt['jobs']={n:f.result() for n,f in futures.items()}
        failed=[n for n,r in receipt['jobs'].items() if r['status']!='PASS_EXIT_ZERO']
        if failed: raise RuntimeError('Retained failed jobs: '+', '.join(failed))
        records={n:json.loads((output/f'{n}.json').read_text(encoding='utf-8')) for n in jobs}
        for n, r in records.items():
            allowed_rounding=(n=='s4_rounding' and r.get('status')=='COMPLETE_DECIMAL100_SPECIFIC_ROUNDING_PROTOCOL_AUDIT')
            if not r.get('status','').startswith('PASS') and not allowed_rounding: raise ValueError('Unexpected result status: '+n)
        expected_nodes={'seikh_q4_mpfr':11,'seikh_q8_mpfr':3,'seikh_q16_mpfr':3,'seikh_uniform_mpfr':811}
        for name,count in expected_nodes.items():
            if name in records and sum(c['node_count'] for c in records[name]['cases'])!=count:
                raise ValueError('Unexpected certificate node count: '+name)
        if records['floor16_decimal']['scores_contained_in_MPFR'] != 256:
            raise ValueError('Decimal score containment incomplete')
        if records['floor16_mpfr']['winner_sequence'] != ['C8','C7','C6','C5']:
            raise ValueError('Wrong finite-node winner sequence')
        if args.profile=='full':
            if (output/'floor16_decimal.json').read_bytes() != (output/'floor16_decimal_optimized.json').read_bytes():
                raise ValueError('Optimized Decimal replay differs')
            if records['seikh_uniform_integer']['nodes'] != 811:
                raise ValueError('Uniform integer cover incomplete')
            for name in ('qol_decimal_cover','qol_mpfr_cover'):
                if [records[name]['groups'][g]['accepted_leaves'] for g in ('primary','secondary')] != [32,61]:
                    raise ValueError('Incomplete QOL cover')
            if records['s4_compact_mpfr']['minimum_genuine_real_q_changes'] != 6 or records['s4_compact_decimal']['at_least_changes'] != 6:
                raise ValueError('Wrong S4 finite-node sequence')
            trials=records['s4_rounding']['trials']
            if [t['digits'] for t in trials] != list(range(6,17)):
                raise ValueError('Missing S4 rounding attempts')
            if any(t['center_all_seven_required_winners'] or t['full_raw_box_all_seven_required_winners'] for t in trials[:4]):
                raise ValueError('S4 six-through-nine-place failures not reproduced')
            if records['s6_cached_diagnostic']['complete_cases'] != 18 or records['s6_cached_diagnostic']['literal_probes'] != 252:
                raise ValueError('Incomplete S6 diagnostic')
            if records['source2022']['finite_node_count'] != 74:
                raise ValueError('Incomplete source2022 screen')
        if check_inputs()!=before:raise ValueError('Input manifest changed')
        receipt.update(status='PASS_SCOPED_THEORY_REPLAY',job_count=len(jobs),all_original_and_wrapper_pins_unchanged=True,
                       outputs={p.name:pin(p) for p in output.iterdir() if p.is_file() and p.suffix=='.json'},
                       uniform_all_finite_q_radius_replayed=args.profile=='full')
    except Exception:
        receipt.update(status='FAIL_PRESERVED',traceback=traceback.format_exc())
        raise
    finally:
        receipt.update(elapsed_seconds=time.perf_counter()-start,finished_utc=datetime.now(timezone.utc).isoformat())
        with (output/'THEORY_REPLAY_RECEIPT.json').open('x',encoding='utf-8') as stream:
            json.dump(receipt,stream,indent=2);stream.write('\n')
        print(json.dumps({k:receipt[k] for k in ('status','profile','elapsed_seconds')},indent=2),flush=True)


if __name__=='__main__':main()
