"""Time-step sensitivity of the probe responses of the 53 parameter sets of the original submission (v1.0).

The v1.0 front was optimised with the probe at dt = 0.5 ms and an SED-increment threshold of 1.1 J m^-3 per step
(rate threshold theta = 2.2 J m^-3 ms^-1). The 53 sets (param.mat of tag v1.0) are recomputed at dt = 0.5, 0.1,
0.0125 and 0.00625 ms with the rate threshold fixed (u_thr = theta * dt), for theta = 2.2 (probe setting of v1.0)
and 0.44 (rate threshold of the v1.0 scanning simulations, 1.1 J m^-3 per 2.5 ms step); T_ref as in v1.0.
For each set: FSL, 1/mean ISI and spike count per speed, J_FSL and J_Rate (as in plt_main_FST2.m /
plt_main_ISI2.m) and dFSL at 9.5 mm/s.

Requires a checkout of tag v1.0 (see README.md).
Usage: python original_sets_dt.py
Output: results/v1.0_sets/probe_dt_original.csv (per set x speed) and probe_dt_original_sets.csv (per set)
"""

import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import lif

OUT = lif.RESULTS / "v1.0_sets"
DTS = (0.5, 0.1, 0.0125, 0.00625)
THETAS = (2.2, 0.44)


def first_spike(s):
    return s[0] if len(s) else 300.0


def rate_hz(s):
    return 1000.0 / np.mean(np.diff(s)) if len(s) > 1 else np.nan


def job(args):
    n, row, theta, dt = args
    out = []
    for k, speed in enumerate(lif.PROBE_SPEEDS):
        s = np.asarray(lif.run_probe(row, k, dt=dt, u_thr=theta * dt))
        out.append(dict(set=n, theta=theta, dt=dt, speed=speed, n_spikes=len(s), fsl=first_spike(s), rate=rate_hz(s)))
    return out


def main():
    params = lif.load_v1_params()
    _, actual = lif.load_probe_inputs()
    exp = [c[~np.isnan(c) & (c > 0)] for c in actual.T]
    fsl_obs = np.array([first_spike(s) for s in exp])
    rate_obs = np.array([rate_hz(s) for s in exp])
    jobs = [(n, row, th, dt) for th in THETAS for dt in DTS for n, row in enumerate(params, 1)]
    workers = int(os.environ.get("RA_N_WORKERS", min(8, os.cpu_count() or 1)))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        rows = [r for res in pool.map(job, jobs, chunksize=4) for r in res]
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "probe_dt_original.csv", index=False)

    sets = []
    for (th, dt, n), g in df.groupby(["theta", "dt", "set"]):
        g = g.set_index("speed").loc[lif.PROBE_SPEEDS]
        fsl, rate = g.fsl.to_numpy(), g.rate.to_numpy()
        sets.append(dict(theta=th, dt=dt, set=n, dFSL=fsl[0] - fsl_obs[0],
                         J_fsl=-(1 - np.sum((fsl_obs - fsl) ** 2) / np.sum(fsl_obs ** 2)),
                         J_rate=-(1 - np.sum((rate_obs - rate) ** 2) / np.sum(rate_obs ** 2)),
                         rms_rate_err=np.sqrt(np.mean((rate_obs - rate) ** 2)),
                         n_spikes_total=int(g.n_spikes.sum())))
    pd.DataFrame(sets).to_csv(OUT / "probe_dt_original_sets.csv", index=False)
    print(f"done: {len(df)} rows -> {OUT}")


if __name__ == "__main__":
    main()
