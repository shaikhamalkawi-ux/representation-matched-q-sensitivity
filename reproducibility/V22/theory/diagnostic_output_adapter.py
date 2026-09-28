"""New-output adapter for the unchanged cached S6 numerical diagnostic.

The literal original whole script is NOT run. Its unchanged 180-bisection body
is replayed through the previously audited identical-argument cache, with
252 literal equality probes and exact comparison to all 18 retained cases.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
if os.name == 'nt' and not str(HERE).startswith('\\\\?\\'):
    HERE = Path('\\\\?\\'+str(HERE))


def write(path, record):
    with path.open('x',encoding='utf-8') as stream:
        json.dump(record,stream,indent=2);stream.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if sys.flags.optimize:
        raise ValueError('Use normal Python: the frozen diagnostic uses assertions')
    path=HERE/'original/outputs/FSS_V19_DRAFT_20260927/diagnostic_optimization/check_joint_diagnostic_cached_R2.py'
    manifest=json.loads((HERE/'THEORY_MANIFEST.json').read_text(encoding='utf-8'))
    expected=manifest['payloads'][path.relative_to(HERE).as_posix()]['sha256']
    actual=hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:raise ValueError('Cached diagnostic source changed')
    spec=importlib.util.spec_from_file_location('frozen_s6_cache',path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    original, body=mod.load_original()
    namespace, body_sha=mod.build_functions(original,body)
    reference=json.loads(mod.EXPECTED.read_text(encoding='utf-8'))
    if [tuple(c[k] for k in ('m','k','n','r')) for c in reference['cases']] != mod.CASE_SPECS:
        raise ValueError('Original case order differs')
    base={'checker_sha256':actual,'original_script_sha256':mod.ORIGINAL_SHA,
          'expected_json_sha256':mod.EXPECTED_SHA,'unchanged_check_core_ast_sha256':body_sha,
          'precision_decimal_digits':90,'bisections_per_target':180,
          'output_adapter':'Public path-free receipt; unchanged cached numerical functions'}
    result={'status':'STARTED','scope':'S6 ordinary Decimal90 numerical diagnostic only',
            'outward_arithmetic':False,'continuum_certificate':False,
            'original_full_script_executed':False}
    try:
        probes=mod.equivalence_probes(namespace,120,base)
        write(args.output.with_name('s6_literal_equivalence.json'),probes)
        if probes['total_cases'] != 18 or probes['total_probes'] != 252:
            raise ValueError('Incomplete literal equivalence probe set')
        receipt, payload=mod.full_compare(namespace,420,reference,base)
        write(args.output.with_name('s6_cached_compare.json'),receipt)
        if payload is not None:write(args.output.with_name('s6_cached_results.json'),payload)
        if receipt['status'] != 'PASS_OPTIMIZED_FULL_READ_COMPARE' or receipt['completed_cases'] != 18:
            raise ValueError('Cached diagnostic did not complete all retained cases')
        result.update(status='PASS_SCOPED_CACHED_DIAGNOSTIC',literal_probes=252,
                      complete_cases=18,all_scientific_values_equal_frozen_reference=True,
                      original_180_bisection_body_preserved=True,
                      original_script_sha256=mod.ORIGINAL_SHA,expected_json_sha256=mod.EXPECTED_SHA,
                      historical_original_full_script_timeouts_not_reclassified=True)
    except Exception as exc:
        result.update(status='FAIL_PRESERVED',failure_type=type(exc).__name__,message=str(exc))
        raise
    finally:write(args.output,result)


if __name__=='__main__':main()
