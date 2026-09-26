"""Prespecified local search, followed by exact integer-interval witness checks.

Exploration is not a lower-bound proof or a global optimization certificate.
The immutable model implementation is imported read-only from the bundled vendor
directory (or an explicit --model-dir/--input directory with identical hashes).
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.optimize import minimize
import scipy

HERE = Path(__file__).resolve().parent
DEFAULT_MODEL = HERE.parent / 'vendor'
EXPECTED_HASHES = {
    'model_inputs.json':'71351416bb73b22f3781bad962a023b6ea7362e8020429aab97b47d3f36d1797',
    'exact_models.py':'8d5941c32cd7cf7140c3273ac6805a370f73bcea627d321cd2ba3c2bfc62f427',
    'dyadic_interval.py':'38cae911d1feb49562dcfd3693bb56a8e1bb8219c8afa51f89d44aff5998852e'
}

Q_VALUES = (4, 8, 16)
RIVALS = (0, 2, 3)
SEED = 177031
BISECTION_STEPS = 22
MAX_ITERATIONS = 400
INITIAL_UPPER_RADIUS = F(3, 100)
EXACT_RADIUS_GRID = 10**8
WITNESS_SAFETY_INCREMENT = F(1, 10**7)
CANDIDATE_GRID = 10**12

def floating_model(x, q, rival=None):
    """Score and optional analytical gradient of S_Y2 - S_rival."""
    raw = np.asarray(x).reshape(4, 5, 2)
    z = raw**4
    p = 4.0/q
    transported = raw**p
    e = 1 + (transported * np.log(transported)).sum(axis=(0, 2))/4
    w = e/e.sum()
    lp = np.log1p(-z[:, :, 0])
    ln = np.log(z[:, :, 1])
    pp = np.exp((w*lp).sum(axis=1))
    pn = np.exp((w*ln).sum(axis=1))
    scores = 1 - pp - pn
    if rival is None:
        return scores
    weight_derivative = -pp[:, None]*lp - pn[:, None]*ln
    a = weight_derivative[1] - weight_derivative[rival]
    e_derivative = p*raw**(p-1)*(p*np.log(raw)+1)/4
    gradient = e_derivative*((a - np.dot(a, w))/e.sum())[None, :, None]
    direct = np.zeros_like(raw)
    for i, sign in ((1, 1), (rival, -1)):
        direct[i, :, 0] = sign*4*raw[i, :, 0]**3*w*pp[i]/(1-z[i, :, 0])
        direct[i, :, 1] = -sign*4*w*pn[i]/raw[i, :, 1]
    return float(scores[1]-scores[rival]), (gradient+direct).ravel()


def search_one(q, rival, h, starts):
    def objective(u):
        gap, gradient = floating_model(BASELINE+h*u, q, rival)
        return gap, h*gradient
    candidates = []
    for start in starts:
        result = minimize(objective, start, method='L-BFGS-B', jac=True,
                          bounds=[(-1.0, 1.0)]*40,
                          options={'ftol':1e-14, 'gtol':1e-10,
                                   'maxiter':MAX_ITERATIONS, 'maxls':40})
        candidates.append(result)
    best = min(candidates, key=lambda result: result.fun)
    return best, {
        'q':q, 'rival':f'Y{rival+1}', 'halfwidth_float':h,
        'minimum_gap_found_float':float(best.fun),
        'starts':[{'success':bool(r.success), 'iterations':int(r.nit),
                   'function_evaluations':int(r.nfev), 'message':str(r.message),
                   'gap_float':float(r.fun)} for r in candidates]
    }


def exact_witness(q, rival, radius, vector):
    candidate = [max(c-radius, min(c+radius, F(round(float(x)*CANDIDATE_GRID), CANDIDATE_GRID)))
                 for c,x in zip(BASELINE_EXACT, vector)]
    canonical = [x**4 for x in candidate]
    strict_interior = all(0<x<1 for x in candidate)
    pair_max = max(canonical[i]+canonical[i+1] for i in range(0,40,2))
    actual_radius = max(abs(x-c) for x,c in zip(candidate, BASELINE_EXACT))
    scores, weights = seikh(IV.frac(F(1,q)), [IV.frac(x) for x in canonical], tight=True)
    margins = [IV(scores[rival].lo-scores[j].hi,
                  scores[rival].hi-scores[j].lo) for j in range(4) if j != rival]
    gap = IV(scores[1].lo-scores[rival].hi, scores[1].hi-scores[rival].lo)
    accepted = strict_interior and pair_max<=1 and actual_radius<=radius and gap.hi<0
    return {
        'status':'PASS_WITNESS' if accepted else 'REJECTED',
        'q_exact':str(q), 'winner':f'Y{rival+1}',
        'raw_baseline_normalized_inputs':[str(x) for x in candidate],
        'radius_limit':str(radius), 'actual_linf_change':str(actual_radius),
        'admissible':strict_interior and pair_max<=1,
        'maximum_canonical_pair_sum':str(pair_max),
        'scores_exact':[s.endpoints() for s in scores],
        'weights_exact':[w.endpoints() for w in weights],
        'margin_Y2_minus_rival':gap.endpoints(),
        'rival_margins_over_other_alternatives':[m.endpoints() for m in margins],
        'rival_is_unique_winner':min(m.lo for m in margins)>0,
        'bits':PREC,
        'interpretation':'Feasible upper witness only, not a global minimum or lower bound.'
    }


def main():
    global IV, SCALE, PREC, INPUTS, domain, seikh, BASELINE_EXACT, BASELINE
    parser = argparse.ArgumentParser()
    parser.add_argument('--model-dir', '--input', dest='model_dir', type=Path,
                        default=DEFAULT_MODEL,
                        help='Directory containing the three hash-pinned model files.')
    parser.add_argument('--output', type=Path, default=HERE/'FIXED_RUNG_WITNESSES.json',
                        help='Fresh output path; existing files are never overwritten.')
    args = parser.parse_args()
    model=args.model_dir.resolve()
    output=args.output.resolve()
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite existing results: {output}')
    if not output.parent.is_dir():
        raise ValueError(f'Output parent does not exist: {output.parent}')
    hashes={name:hashlib.sha256((model/name).read_bytes()).hexdigest() for name in EXPECTED_HASHES}
    if hashes != EXPECTED_HASHES:
        raise ValueError('Immutable model hash mismatch')
    sys.path.insert(0,str(model))
    from dyadic_interval import IV, SCALE, PREC
    from exact_models import INPUTS, domain, seikh
    BASELINE_EXACT=[F(n, INPUTS['seikh']['denominator']) for n in INPUTS['seikh']['numerators']]
    BASELINE=np.array([float(x) for x in BASELINE_EXACT])
    started = time.perf_counter()
    domain('seikh', 'source', INITIAL_UPPER_RADIUS)
    rng = np.random.default_rng(SEED)
    starts = [np.zeros(40)] + [rng.uniform(-1,1,40) for _ in range(3)]
    plan = {
        'q_values':list(Q_VALUES), 'rivals':[f'Y{i+1}' for i in RIVALS],
        'seed':SEED, 'starts':'Baseline and three seeded uniform points in each normalized box.',
        'initial_upper_radius':str(INITIAL_UPPER_RADIUS),
        'bisection_steps':BISECTION_STEPS, 'max_iterations_per_local_run':MAX_ITERATIONS,
        'method':'L-BFGS-B with analytical gradient and four starts; 22-step heuristic bisection after a negative gap is found.',
        'safety_increment':str(WITNESS_SAFETY_INCREMENT),
        'scope':'A fixed-rung witness search in the unchanged 40-component normalized baseline-raw L-infinity geometry, baseline rung 4.',
        'warning':'An unsuccessful local search is not evidence of nonexistence. Bisection lower endpoints are never certified lower bounds.'
    }
    # Check the hand-derived gradient before using it for search. This is numerical QA,
    # not part of the final exact witness acceptance rule.
    gradient_checks = []
    for q in Q_VALUES:
        for rival in RIVALS:
            _, analytic = floating_model(BASELINE, q, rival)
            finite = np.empty(40)
            for k in range(40):
                step = np.zeros(40);step[k]=1e-6
                finite[k] = (floating_model(BASELINE+step,q,rival)[0]-
                             floating_model(BASELINE-step,q,rival)[0])/(2e-6)
            error = float(np.max(np.abs(finite-analytic)))
            if error>2e-8:
                raise RuntimeError(f'gradient QA failed {q} {rival}: {error}')
            gradient_checks.append({'q':q, 'rival':rival+1, 'max_absolute_error':error})
    records = []
    witnesses = []
    for q in Q_VALUES:
        for rival in RIVALS:
            hi=float(INITIAL_UPPER_RADIUS);lo=0.0
            result, record = search_one(q,rival,hi,starts)
            records.append(record)
            if result.fun>=0:
                print({'q':q,'rival':rival+1,'status':'NO_WITNESS_FOUND_AT_UPPER_SEARCH_RADIUS'},flush=True)
                continue
            for _ in range(BISECTION_STEPS):
                h=(lo+hi)/2
                trial, record=search_one(q,rival,h,starts)
                records.append(record)
                if trial.fun<0:
                    hi=h;result=trial
                else:
                    lo=h
            radius=F(int(np.ceil(hi*EXACT_RADIUS_GRID)),EXACT_RADIUS_GRID)+WITNESS_SAFETY_INCREMENT
            final, record=search_one(q,rival,float(radius),starts)
            records.append(record)
            witness=exact_witness(q,rival,radius,BASELINE+float(radius)*final.x)
            witness['heuristic_transition_interval_float']=[lo,hi]
            witnesses.append(witness)
            print({k:witness[k] for k in ['q_exact','winner','status','actual_linf_change','rival_is_unique_winner']},flush=True)
    result = {
        'plan':plan, 'gradient_checks':gradient_checks, 'records':records, 'witnesses':witnesses,
        'immutable_input_sha256':hashes,
        'versions':{'numpy':np.__version__,'scipy':scipy.__version__,'python':sys.version},
        'elapsed_seconds':time.perf_counter()-started,
        'summary_status':'PASS' if witnesses and all(w['status']=='PASS_WITNESS' for w in witnesses) else 'INCOMPLETE'
    }
    with output.open('x',encoding='utf-8') as stream:
        stream.write(json.dumps(result,indent=2)+'\n')
    print({'output':str(output),'elapsed_seconds':result['elapsed_seconds']},flush=True)


if __name__=='__main__':
    main()
