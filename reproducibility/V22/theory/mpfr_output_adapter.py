"""Output-only adapter for the unchanged directed-MPFR certificate verifier.

The original verifier's CLI writes alongside its source. This adapter imports
its unchanged verify(record) function and writes exclusively to a new caller-
selected output, preserving all original code/input hashes and mathematics.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent/'original'
if os.name == 'nt' and not str(ROOT).startswith('\\\\?\\'):
    ROOT = Path('\\\\?\\'+str(ROOT))
CHECKER = ROOT/'outputs/FSS_STRENGTHENING_20260926/radius_refinement/mpfr_audit/verify_refined_mpfr.py'
EXPECTED = '867b871950cf41960ee1274411b2581cda88377ecdc512617991c1471487592a'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--record', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Output exists')
    if hashlib.sha256(CHECKER.read_bytes()).hexdigest() != EXPECTED:
        raise ValueError('Original verifier pin mismatch')
    spec = importlib.util.spec_from_file_location('unchanged_public_mpfr_verifier', CHECKER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    start = time.perf_counter()
    raw = args.record.read_bytes()
    result = module.verify(json.loads(raw))
    result.update(record_sha256=hashlib.sha256(raw).hexdigest(), checker_sha256=EXPECTED,
                  elapsed_seconds=time.perf_counter()-start, adapter='output_destination_only_no_formula_change')
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')
    print(json.dumps({'status':result['status'], 'case_count':len(result['cases']),
                      'nodes':sum(c['node_count'] for c in result['cases'])}), flush=True)


if __name__ == '__main__':
    main()
