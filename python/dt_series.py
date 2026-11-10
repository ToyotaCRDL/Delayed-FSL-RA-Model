"""Probe responses of the revised sets at a series of time steps (choice of dt = 0.0125 ms).

The 150 sets of the revised implementation (theta 0.13 J m^-3 ms^-1, T_ref 2.1 ms; fitted at 0.0125 ms) are
simulated on the probe (6 speeds) at each time step, halving from 0.5 to 0.0125 ms (0.0625 to 0.025 ms: factor 2.5).

Usage: python dt_series.py [NAME]   (default: unitA)
Output: results/<NAME>/dt_series/probe_dt_series.csv (one row per time step, set and speed, with all spike
        times) and summary.csv (each halving: conditions whose spike count changed, largest FSL and rate changes;
        conditions differing from 0.0125 ms)
"""

import sys

import numpy as np
import pandas as pd

import lif

RES = lif.RESULTS
THETA = lif.THETA
STEPS = [0.5, 0.25, 0.125, 0.0625, 0.025, 0.0125]


def run(params, dt):
    data, _ = lif.load_probe_inputs()
    X = params[["beta_mA", "kd", "t_ref_ms", "v_thr_mV"]].to_numpy()
    rows = []
    for k, ps in enumerate(lif.PROBE_SPEEDS):
        t, sed = lif.resample_sed(data, k, dt=dt)
        cur = lif.afferent_current(sed[:, None], t[1] - t[0], X[:, 0], X[:, 1], THETA * dt)
        for p, s in zip(params.set, lif.integrate(t, cur, X[:, 2], X[:, 3])):
            rows.append(dict(dt=dt, set=p, speed=ps, nspk=len(s), fsl=s[0] if len(s) else np.nan,
                             rate=1000 / np.mean(np.diff(s)) if len(s) > 1 else np.nan,
                             spikes=" ".join(f"{x:.6g}" for x in s)))
    return pd.DataFrame(rows)


def main(folder):
    folder = RES / folder
    params = pd.read_csv(folder / "params.csv")
    out = folder / "dt_series"
    out.mkdir(exist_ok=True)
    df = pd.concat([run(params, dt) for dt in STEPS], ignore_index=True)
    df.to_csv(out / "probe_dt_series.csv", index=False)
    piv = {c: df.pivot_table(index=["set", "speed"], columns="dt", values=c) for c in ["nspk", "fsl", "rate"]}
    rows = []
    for coarse, fine in zip(STEPS[:-1], STEPS[1:]):
        changed = piv["nspk"][coarse] != piv["nspk"][fine]
        rows.append(dict(coarse=coarse, fine=fine, changed_conditions=int(changed.sum()),
                         changed_sets=int(changed.groupby("set").any().sum()),
                         changed_speeds=" ".join(sorted(set(changed[changed].index.get_level_values("speed")))),
                         max_fsl_diff=np.nanmax(np.abs(piv["fsl"][coarse] - piv["fsl"][fine])),
                         max_rate_diff=np.nanmax(np.abs(piv["rate"][coarse] - piv["rate"][fine]))))
    summary = pd.DataFrame(rows)
    finest = piv["nspk"][STEPS[-1]]
    summary["differs_from_finest"] = [int((piv["nspk"][c] != finest).sum()) for c in summary.coarse]
    summary.to_csv(out / "summary.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "unitA")
