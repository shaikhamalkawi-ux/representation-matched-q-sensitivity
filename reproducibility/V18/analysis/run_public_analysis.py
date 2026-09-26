#!/usr/bin/env python3
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser(description='Run public Qin/He verification workflow.')
    g=ap.add_mutually_exclusive_group()
    g.add_argument('--radius-replay',action='store_true',help='Also replay the retained 18,515-node radius certificate (slow).')
    g.add_argument('--full',action='store_true',help='Rebuild the radius-cover search and replay it (very slow).')
    args=ap.parse_args()
    steps=[
        ['run_stagewise.py'],
        ['verify_new_analysis.py'],
        ['verify_structural_controls.py'],
        ['verify_checker_negative.py'],
        ['verify_upper_witness_directed.py'],
    ]
    if args.full:
        steps.append(['certify_radius.py','--radius','0.000995','--max','30000'])
        steps.append(['replay_certificate.py'])
    elif args.radius_replay:
        steps.append(['replay_certificate.py'])
    for step in steps:
        print(f"\n==> {step[0]}", flush=True)
        subprocess.run([sys.executable,*step],cwd=HERE,check=True)
    print('\nPUBLIC QIN/HE ANALYSIS: PASS')

if __name__=='__main__':
    main()
