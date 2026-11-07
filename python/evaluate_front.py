"""Probe and scanning responses of every set of a front, with the number of responding afferents.

For every set in results/<NAME>/params.csv the probe (6 speeds), the dot pattern (6 spacings) and the single dots
(3 widths) are simulated with the revised implementation (rate threshold theta, the same dt for all stimuli).
Besides FR/TV RMSE and RW, the number of responding afferents (at least one spike) and of spikes is stored for every
scanning condition. Used to compare the equal-weight and relative-error-weight fronts (relative_weights.py).

Usage: python evaluate_front.py [NAME]   (default: unitA; also unitA_relative_weights)
Output: results/<NAME>/front_summary.csv
"""

import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import lif
import scan


def first_spike(s):
    return s[0] if len(s) else 300.0


def rate_hz(s):
    return 1000.0 / np.mean(np.diff(s)) if len(s) > 1 else np.nan


def fss(obs, sim):
    return -(1 - np.sum((obs - sim) ** 2) / np.sum(obs**2))


def evaluate(args):
    n, row, theta, dt = args
    probe_data, actual = lif.load_probe_inputs()
    exp = [c[~np.isnan(c) & (c > 0)] for c in actual.T]
    fsl_obs = np.array([first_spike(s) for s in exp])
    rate_obs = np.array([rate_hz(s) for s in exp])
    u_thr = theta * dt
    rec = {"set": n}
    fsl, rate = np.zeros(6), np.zeros(6)
    for k, ps in enumerate(lif.PROBE_SPEEDS):
        t, sed = lif.resample_sed(probe_data, k, dt=dt)
        s = lif.simulate(t, sed, *row, u_thr=u_thr)
        fsl[k], rate[k] = first_spike(s), rate_hz(s)
        rec[f"fsl_{ps}"], rec[f"rate_{ps}"], rec[f"nspk_{ps}"] = fsl[k], rate[k], len(s)
    rec["dFSL"] = fsl[0] - fsl_obs[0]
    rec["J_fsl"], rec["J_rate"] = fss(fsl_obs, fsl), fss(rate_obs, rate)
    rec["rms_rate_err"] = np.sqrt(np.nanmean((rate_obs - rate) ** 2)) if not np.isnan(rate).any() else np.nan
    rec["rms_fsl_err"] = np.sqrt(np.mean((fsl_obs - fsl) ** 2))

    exp_fr, exp_tv = scan.load_dot_texture_experiment()
    fr, tv = np.zeros(6), np.zeros(6)
    for mlp, (ds, finger) in enumerate(scan.DOT_SPACINGS):
        trains = scan.run(scan.DOT_TEXTURE_DIR, finger, row, u_thr=u_thr, dt=dt)
        fr[5 - mlp], tv[5 - mlp] = scan.dot_texture_metrics(trains, finger)
        rec[f"fr_{ds}"], rec[f"tv_{ds}"] = fr[5 - mlp], tv[5 - mlp]
        rec[f"resp_{ds}"] = sum(len(s) > 0 for s in trains)
        rec[f"nspkscan_{ds}"] = sum(len(s) for s in trains)
    rec["fr_rmse"], rec["tv_rmse"] = scan.rmse(exp_fr, fr), scan.rmse(exp_tv, tv)
    for dw, finger in scan.DOT_WIDTHS:
        trains = scan.run(scan.SINGLE_DOT_DIR, finger, row, u_thr=u_thr, dt=dt)
        try:
            rec[f"rw_{dw}"] = scan.response_width(trains, finger)
        except (IndexError, ValueError):
            rec[f"rw_{dw}"] = np.nan
        rec[f"resp_single_{dw}"] = sum(len(s) > 0 for s in trains)
        rec[f"nspk_single_{dw}"] = sum(len(s) for s in trains)
    return rec


def main(name, theta=lif.THETA, dt=lif.DT):
    folder = lif.RESULTS / name
    p = pd.read_csv(folder / "params.csv")
    jobs = [(int(r.set), tuple(r[lif.PARAM_COLUMNS]), theta, dt) for _, r in p.iterrows()]
    workers = int(os.environ.get("RA_N_WORKERS", os.cpu_count() or 1))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        records = list(pool.map(evaluate, jobs))
    df = pd.DataFrame(records).sort_values("set")
    df.to_csv(folder / "front_summary.csv", index=False)
    print(f"{name}: {len(df)} sets evaluated at dt={dt:g}, theta={theta:g} -> {folder / 'front_summary.csv'}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "unitA")
