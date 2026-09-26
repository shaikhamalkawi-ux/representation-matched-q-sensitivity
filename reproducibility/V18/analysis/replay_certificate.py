#!/usr/bin/env python3
"""Reconstruct a saved binary cover and re-evaluate its leaves at higher precision.

The search heuristic is not executed. Both children of every pending node must
be present; contraction is replayed with its declared arithmetic and every
positive leaf is independently re-evaluated at 256 bits. The same small MPFR
kernel remains in the trusted computing base; this is not a proof assistant.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json,time
from pathlib import Path
from fractions import Fraction
import numpy as np
import kernel as K
from certify_radius import CENTRE, pair_idx, rootbox, state_margin
if not __debug__: raise RuntimeError('Verification requires assertions; do not use Python -O')
HERE=Path(__file__).resolve().parent

def replay(directory:Path):
    started=time.time();trace=directory/'COVER_TREE.jsonl.gz'
    original=json.loads((directory/'RADIUS_SUMMARY.json').read_text())
    nodes={};splits={};children={};leaf_lower=[];outside=0;end=None
    with gzip.open(trace,'rt') as f:
        header=json.loads(next(f));assert header['kind']=='header'
        assert header['challenger']==2, 'The full-winner cover must certify A1 against A3'
        assert header['radius_decimal']==original['radius_decimal']
        assert header['radius_dyadic_hex']==original['radius_dyadic_hex']
        assert header['bits']==original['bits']==160
        assert original['status']=='PASS'
        radius=float.fromhex(header['radius_dyadic_hex'])
        assert Fraction.from_float(radius)>=Fraction(header['radius_decimal'])
        idx=pair_idx(header['challenger']);assert idx.tolist()==header['idx']
        rootlo,roothi=rootbox(idx,radius)
        for line in f:
            row=json.loads(line);kind=row['kind']
            if kind=='end': end=row;continue
            node=row['id']
            if kind=='split':
                assert node in nodes and nodes[node]['kind']=='pending' and node not in splits
                lo,hi=nodes[node]['box'];k=row['coordinate'];mid=float.fromhex(row['midpoint_hex'])
                assert 0<=k<len(idx) and lo[k]<mid<hi[k]
                splits[node]=(k,mid);children[node]={};continue
            assert node not in nodes and kind in ('pending','positive','outside')
            parent,side=row['parent'],row['side']
            if parent is None:
                assert not nodes and node==0 and side is None
                lo,hi=rootlo.copy(),roothi.copy()
            else:
                assert parent in splits and side in ('L','R') and side not in children[parent]
                lo,hi=(x.copy() for x in nodes[parent]['box']);k,mid=splits[parent]
                if side=='L':hi[k]=mid
                else:lo[k]=mid
                children[parent][side]=node
            K.lib.kernel_precision(header['bits']);box=K.contract(idx,lo,hi,radius)
            if kind=='outside':
                assert box is None;outside+=1;nodes[node]={'kind':kind};continue
            assert box is not None
            nodes[node]={'kind':kind,'box':box}
            if kind=='positive':
                K.lib.kernel_precision(256);margin=state_margin(idx,*box,header['challenger'])
                assert margin[0]>0, (node,margin)
                leaf_lower.append(float(margin[0]))
    assert end is not None and end['pass_'] and end['queue']==0
    assert len(nodes)==end['nodes']==original['nodes']
    pending={n for n,v in nodes.items() if v['kind']=='pending'}
    assert pending==set(splits)
    assert all(set(x)=={'L','R'} for x in children.values())
    assert len(nodes)==2*len(splits)+1
    assert len(leaf_lower)==end['positive_leaves'] and outside==end['outside_leaves']
    other=[];K.lib.kernel_precision(256)
    for j in (1,3,4):
        ii=pair_idx(j);lo,hi=rootbox(ii,radius);g=state_margin(ii,lo,hi,j)
        assert g[0]>0;other.append({'challenger':f'A{j+1}','margin_lower':float(g[0]),'margin_upper':float(g[1])})
    result={'status':'PASS','radius_decimal':header['radius_decimal'],'radius_dyadic_hex':header['radius_dyadic_hex'],
      'cover_nodes':len(nodes),'binary_splits':len(splits),'positive_leaves':len(leaf_lower),'outside_leaves':outside,
      'unresolved_nodes':0,'minimum_256_bit_leaf_lower':min(leaf_lower),'leaf_evaluation_bits':256,
      'cover_contraction_bits':header['bits'],'other_challengers':other,'elapsed_seconds':time.time()-started,
      'trace_sha256':hashlib.sha256(trace.read_bytes()).hexdigest(),
      'trust_boundary':'Independent cover reconstruction and higher-precision leaf re-evaluation, sharing the audited MPFR interval kernel; not a proof-assistant check.'}
    (directory/'REPLAY_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('directory',nargs='?',type=Path,default=HERE/'radius_000995');a=ap.parse_args();replay(a.directory)
