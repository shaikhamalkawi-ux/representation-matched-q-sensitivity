"""Path/provenance-only wrapper around unchanged 2022 baseline and screen.

No source formula, data, score, grid or precision is edited. All scientific
fields are compared exactly to explicit projections of the frozen records.
"""
import argparse
import contextlib
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
if os.name=='nt' and not str(HERE).startswith('\\\\?\\'):
    HERE=Path('\\\\?\\'+str(HERE))
SOURCE=HERE/'original/outputs/FSS_OPEN_ENDED_RESEARCH_20260927/source2022'


def write(path,record):
    with path.open('x',encoding='utf-8') as stream:
        json.dump(record,stream,indent=2);stream.write('\n')


def load(name):
    path=SOURCE/name
    manifest=json.loads((HERE/'THEORY_MANIFEST.json').read_text(encoding='utf-8'))
    if hashlib.sha256(path.read_bytes()).hexdigest()!=manifest['payloads'][path.relative_to(HERE).as_posix()]['sha256']:
        raise ValueError('Changed original script: '+name)
    spec=importlib.util.spec_from_file_location(name[:-3],path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def compare(actual,name,allowed):
    expected=json.loads((HERE/'expected'/name).read_text(encoding='utf-8'))
    if set(expected['excluded_top_level_fields'])!=set(allowed):
        raise ValueError('Unexpected comparison exclusion schema')
    if any(key not in actual for key in allowed):raise ValueError('Missing source metadata')
    projection={key:value for key,value in actual.items() if key not in allowed}
    if projection!=expected['scientific_result']:
        raise ValueError('Scientific result differs: '+name)
    return len(projection)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    baseline_path=args.output.with_name('source2022_baseline.json')
    screen_path=args.output.with_name('source2022_screen.json')
    result={'status':'STARTED','formula_or_numerical_input_modified':False,
            'output_adapter':'New output paths and public operator-provenance note only',
            'source_is_html_not_authenticated_pdf':True,'continuum_certificate':False}
    try:
        baseline=load('seikh2022_html_diagnostic_r2.py')
        baseline.INPUT=HERE/'data/SOURCE2022_INPUT.json'
        baseline.OPERATOR=HERE/'SOURCE2022_OPERATOR_SCOPE.md'
        baseline.OUT=baseline_path
        # Old mains print host output paths; suppress only printing, not checks.
        with contextlib.redirect_stdout(io.StringIO()):baseline.main()
        data=json.loads(baseline_path.read_text(encoding='utf-8'))
        nbase=compare(data,'source2022_baseline_expected.json',
                      ['utc','input_name','input_sha256','operator_provenance_name','operator_provenance_sha256','formula_provenance'])
        screen=load('seikh2022_matched_screen_r1.py')
        screen.SOURCE=HERE/'data/SOURCE2022_INPUT.json'
        screen.OPERATOR=HERE/'SOURCE2022_OPERATOR_SCOPE.md'
        screen.BASELINE=baseline_path
        screen.OUTPUT=screen_path
        with contextlib.redirect_stdout(io.StringIO()):screen.main()
        screened=json.loads(screen_path.read_text(encoding='utf-8'))
        nscreen=compare(screened,'source2022_screen_expected.json',['utc','input_pins'])
        if screened['grid']['finite_node_count']!=74:raise ValueError('Incomplete frozen grid')
        for kind in ('WA','WG'):
            if screened['summary'][kind]['observed_adjacent_winner_changes']!=[]:
                raise ValueError('Unexpected source-screen winner changes')
        result.update(status='PASS_ORDINARY_SOURCE2022_REPLAY',finite_node_count=74,
                      baseline_exact_scientific_top_level_fields=nbase,
                      screen_exact_scientific_top_level_fields=nscreen,
                      scientific_payloads_equal_frozen_records=True,
                      summary=screened['summary'],source_doi='10.1007/s41066-021-00290-2',
                      failures_preserved='All original baseline rounding comparisons, including failures and post-hoc truncation diagnostics, remain in source2022_baseline.json')
    except Exception as exc:
        result.update(status='FAIL_PRESERVED',failure_type=type(exc).__name__,message=str(exc))
        raise
    finally:write(args.output,result)


if __name__=='__main__':main()
