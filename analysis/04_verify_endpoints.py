"""Independent verification of the profile-likelihood interval endpoints from Script 03.

Purpose
-------
Script 03 builds each profile by warm-started outward walks plus a small fixed grid of
cold starts (seed 12345). A warm-started walk can follow one branch of the likelihood
surface and miss a lower one, which would make an interval look narrower than it is.
This script re-checks every reported endpoint WITHOUT using any Script 03 solution:

  * reference minima are recomputed from random starts only;
  * at each finite endpoint e the profile is re-evaluated at e - d, e and e + d
    (d = 0.005 on the log scale for beta and Q, 0.25 day for the change time),
    each from N_STARTS Latin-hypercube starts spread over the WHOLE search box,
    using a different random seed (20271);
  * the open upper end of Q in Model B is re-checked at Q = 1e7 and in the exact
    Q -> infinity (no-depletion) limit.

A continuous refinement of the change-time minimum is also made (bounded scalar search
around the best cold-start day, warm-started from that cold fit) because Script 03 located
that minimum on a 2-day grid. Change-time endpoints are re-derived against both references.

Pass criteria for a finite endpoint (LR threshold T = chi2_1(0.95)/2 = 1.9207 on the NLL scale):
  1. delta(inside) < T < delta(outside)  - the threshold is crossed between e - d and e + d;
  2. the recomputed reference minimum is not lower than Script 03's by more than 0.01;
  3. the endpoint implied by linear interpolation lies within d of the reported value.

The same likelihood, model and search box as Script 03 are used (imported from it), so this
checks the optimisation, not the modelling assumptions.

Run:  python analysis/04_verify_endpoints.py      (about 10-20 min on 2 cores)
"""
import importlib.util
import json
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import qmc

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location('profile03', HERE / '03_profile_likelihood.py')
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)          # Script 03's main block is guarded; only definitions load

SEED = 20271
N_STARTS = 80
T = m.THRESH
D_LOG, D_TAU = 0.005, 0.25
REPORTED = json.loads((HERE.parent / 'results' / 'profile_likelihood_summary.json').read_text())


def random_starts(free, n, seed):
    """Latin-hypercube starts over the full log-scale search box (no information from Script 03)."""
    lo = np.array([m.box(k)[0] for k in free]); hi = np.array([m.box(k)[1] for k in free])
    u = qmc.LatinHypercube(d=len(free), seed=seed).random(n)
    x = lo + u * (hi - lo)
    return [dict(zip(free, row)) for row in x]


def cold_profile(job):
    """Profile NLL at one point from random starts only. job = (label, model, param, value, tau, seed)."""
    label, model, param, value, tau, seed = job
    t0 = time.time()
    if param == 'tau':
        free, fixed, tau_use = m.keys_for(model), {}, value
    elif param is None:                                  # unconstrained reference minimum
        free, fixed, tau_use = m.keys_for(model), {}, tau
    elif param == 'Q' and not np.isfinite(value):        # exact Q -> infinity limit
        free, fixed, tau_use = m.keys_for(model, q_infinite=True), {'Q': np.inf}, tau
    else:
        free = [k for k in m.keys_for(model) if k != param]
        fixed, tau_use = {param: value}, tau
    r = m.optimise(fixed, free, model, tau_use, random_starts(free, N_STARTS, seed))
    return dict(label=label, model=model, param=param, value=value, tau=tau_use, nll=r['nll'],
                n_starts=r['n_starts'], n_at_best=r['n_at_best'], polished=r['polished'],
                proj_grad_norm=r['proj_grad_norm'], nuisance_at_bound=r['nuisance_at_bound'],
                fit={k: (float(v) if np.isfinite(v) else None) for k, v in r['logp'].items()},
                seconds=time.time() - t0)


def build_jobs():
    jobs, s = [], SEED
    def add(*a):
        nonlocal s
        s += 1
        jobs.append(a + (s,))
    # reference minima, recomputed cold
    add('ref A', 'A', None, None, None)
    add('ref B tau=56', 'B', None, None, m.TAU_PRIMARY)
    for tau in (52.0, 53.0, 54.0, 55.0, 56.0):
        add('ref B tau-free', 'B', 'tau', tau, None)
    # open end of Q in model B
    add('B: Q open end', 'B', 'Q', np.log(1e7), m.TAU_PRIMARY)
    add('B: Q open end', 'B', 'Q', np.inf, m.TAU_PRIMARY)
    # finite endpoints
    spec = {'A: Q': ('A', 'Q', None), 'A: beta0': ('A', 'b0', None),
            'B: Q': ('B', 'Q', m.TAU_PRIMARY), 'B: beta_before': ('B', 'b0', m.TAU_PRIMARY),
            'B: beta_after': ('B', 'b1', m.TAU_PRIMARY), 'B: change time': ('B', 'tau', None)}
    endpoints = []
    for name, (model, param, tau) in spec.items():
        info = REPORTED['profiles'][name]
        for side, val, status in zip(('lower', 'upper'), info['interval_95_LR'], info['interval_status']):
            if not status.startswith('found'):
                continue
            e = val if param == 'tau' else np.log(val)
            d = D_TAU if param == 'tau' else D_LOG
            inside, outside = (e + d, e - d) if side == 'lower' else (e - d, e + d)
            endpoints.append(dict(name=name, side=side, model=model, param=param, e=e, d=d,
                                  inside=inside, outside=outside, reported=val,
                                  reported_ref=info['nll_min_used']))
            for where, v in (('inside', inside), ('at', e), ('outside', outside)):
                add(f'{name}|{side}|{where}', model, param, v, tau)
    return jobs, endpoints


def natural(param, v):
    return v if param == 'tau' else float(np.exp(v))


if __name__ == '__main__':
    t0 = time.time()
    jobs, endpoints = build_jobs()
    print(f'{len(jobs)} cold profile evaluations x {N_STARTS} random starts (seed {SEED})', flush=True)
    with Pool(2) as pool:
        res = pool.map(cold_profile, jobs, chunksize=1)
    by = {}
    for r in res:
        by.setdefault(r['label'], []).append(r)

    # recomputed reference minima; any lower value seen anywhere in the same setting also counts
    refA = min(r['nll'] for r in res if r['model'] == 'A')
    refB56 = min([by['ref B tau=56'][0]['nll']] +
                 [r['nll'] for r in res if r['model'] == 'B' and r['param'] in ('Q', 'b0', 'b1')])
    refBtau = min(r['nll'] for r in res if r['model'] == 'B')        # tau-free includes every fit
    refs = {'A': refA, 'B56': refB56, 'Btau': refBtau}
    rep = REPORTED['profiles']
    ref_check = {
        'Model A': dict(recomputed=refA, script03=rep['A: Q']['nll_min_used']),
        'Model B, change fixed at 1 Jan (day 56)': dict(recomputed=refB56, script03=rep['B: Q']['nll_min_used']),
        'Model B, change time free': dict(recomputed=refBtau, script03=rep['B: change time']['nll_min_used']),
    }
    for v in ref_check.values():
        v['difference'] = v['recomputed'] - v['script03']
        v['pass'] = bool(v['difference'] > -0.01)

    rows, verdicts = [], []
    for ep in endpoints:
        ref = refs['A'] if ep['model'] == 'A' else (refs['Btau'] if ep['param'] == 'tau' else refs['B56'])
        vals = {w: by[f"{ep['name']}|{ep['side']}|{w}"][0] for w in ('inside', 'at', 'outside')}
        dl = {w: vals[w]['nll'] - ref for w in vals}
        # linear interpolation of the T crossing between inside and outside points
        x_in, x_out = ep['inside'], ep['outside']
        if dl['outside'] != dl['inside']:
            x_cross = x_in + (T - dl['inside']) * (x_out - x_in) / (dl['outside'] - dl['inside'])
        else:
            x_cross = np.nan
        c1 = bool(dl['inside'] < T < dl['outside'])
        c3 = bool(np.isfinite(x_cross) and abs(x_cross - ep['e']) <= ep['d'])
        ok = c1 and c3
        for w in vals:
            rows.append(dict(profile=ep['name'], side=ep['side'], point=w,
                             value=natural(ep['param'], {'inside': x_in, 'at': ep['e'], 'outside': x_out}[w]),
                             delta_nll=dl[w], n_starts=vals[w]['n_starts'], n_at_best=vals[w]['n_at_best'],
                             proj_grad_norm=vals[w]['proj_grad_norm'],
                             nuisance_at_bound=';'.join(vals[w]['nuisance_at_bound'])))
        verdicts.append(dict(profile=ep['name'], side=ep['side'], reported=ep['reported'],
                             recomputed=natural(ep['param'], x_cross) if np.isfinite(x_cross) else None,
                             delta_at_reported=dl['at'], threshold=T,
                             crossing_bracketed=c1, within_tolerance=c3, pass_=ok,
                             min_starts_agreeing=int(min(v['n_at_best'] for v in vals.values()))))

    # continuous refinement of the change-time minimum (Script 03 used a 2-day grid there)
    best_tau = min(by['ref B tau-free'], key=lambda r: r['nll'])
    warm = [{k: v for k, v in best_tau['fit'].items()}]
    prof_tau = lambda t: m.profile_point('tau', t, 'B', None, warm, 0)['nll']
    rs = minimize_scalar(prof_tau, bounds=(best_tau['value'] - 1.0, best_tau['value'] + 1.0),
                         method='bounded', options={'xatol': 0.01})
    tau_refined = dict(change_day=float(rs.x), nll=float(rs.fun),
                       date=str(pd.Timestamp('2023-11-06') + pd.Timedelta(days=float(rs.x))))
    ref_tau_final = min(refs['Btau'], float(rs.fun))
    for v in verdicts:
        if v['profile'] != 'B: change time':
            continue
        pts = {r['point']: r for r in rows if r['profile'] == v['profile'] and r['side'] == v['side']}
        d_in = pts['inside']['delta_nll'] + refs['Btau'] - ref_tau_final
        d_out = pts['outside']['delta_nll'] + refs['Btau'] - ref_tau_final
        x = pts['inside']['value'] + (T - d_in) * (pts['outside']['value'] - pts['inside']['value']) / (d_out - d_in)
        v['recomputed_vs_refined_minimum'] = float(x)
        v['date_vs_refined_minimum'] = str(pd.Timestamp('2023-11-06') + pd.Timedelta(days=float(x)))

    open_end = [dict(Q=('infinity' if not np.isfinite(r['value']) else float(np.exp(r['value']))),
                     delta_nll=r['nll'] - refs['B56'], below_threshold=bool(r['nll'] - refs['B56'] < T),
                     n_at_best=r['n_at_best'], n_starts=r['n_starts']) for r in by['B: Q open end']]
    tau_ref = [dict(change_day=r['value'], nll=r['nll'], n_at_best=r['n_at_best']) for r in by['ref B tau-free']]

    endpoints_pass = all(v['pass_'] for v in verdicts)
    refs_pass = all(v['pass'] for v in ref_check.values())
    open_pass = all(o['below_threshold'] for o in open_end)
    all_pass = endpoints_pass and refs_pass and open_pass
    out = dict(purpose='Independent cold-start verification of Script 03 interval endpoints',
               seed=SEED, random_starts_per_point=N_STARTS, start_design='Latin hypercube over full log search box',
               lr_threshold_nll=T, offsets=dict(log_scale=D_LOG, change_time_days=D_TAU),
               reference_minima=ref_check, change_time_reference_scan=tau_ref,
               change_time_refined_minimum=tau_refined,
               model_B_Q_open_end=open_end, endpoints=verdicts, endpoints_pass=bool(endpoints_pass),
               reference_minima_pass=bool(refs_pass), open_end_pass=bool(open_pass),
               all_checks_pass=bool(all_pass),
               runtime_seconds=time.time() - t0)
    R = HERE.parent / 'results'
    (R / 'endpoint_verification.json').write_text(json.dumps(out, indent=2))
    pd.DataFrame(rows).to_csv(R / 'endpoint_verification_points.csv', index=False)
    print(json.dumps({k: out[k] for k in ('reference_minima', 'change_time_refined_minimum', 'model_B_Q_open_end', 'endpoints',
                                     'endpoints_pass', 'reference_minima_pass', 'open_end_pass', 'all_checks_pass')},
                     indent=1, default=str))
    print(f'done in {time.time() - t0:.0f} s')
