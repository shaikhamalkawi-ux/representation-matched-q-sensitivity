#!/usr/bin/env python3
from __future__ import annotations
import argparse, subprocess, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
    ap=argparse.ArgumentParser(description='Run the public QROF reproduction workflow.')
    g=ap.add_mutually_exclusive_group()
    g.add_argument('--radius-replay',action='store_true',help='Replay the retained Qin radius certificate when the Zenodo cover trace is present (slow).')
    g.add_argument('--full',action='store_true',help='Rebuild and replay the Qin radius cover (very slow).')
    args=ap.parse_args()
    print('==> Seven-case benchmark')
    subprocess.run([sys.executable,'run_all.py'],cwd=ROOT/'benchmark',check=True)
    print('\n==> Qin/He robustness analyses')
    cmd=[sys.executable,'run_public_analysis.py']
    if args.radius_replay: cmd.append('--radius-replay')
    if args.full: cmd.append('--full')
    subprocess.run(cmd,cwd=ROOT/'analysis',check=True)
    print('\nQROF PUBLIC REPRODUCTION: PASS')

if __name__=='__main__':
    main()
