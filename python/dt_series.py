"""Probe responses of the revised sets at a series of time steps (choice of dt = 0.0125 ms).

The 150 sets of the revised implementation (theta 0.13 J m^-3 ms^-1, T_ref 2.1 ms; fitted at 0.0125 ms) are
simulated on the probe (6 speeds) at each time step, halving from 0.5 ms (0.0625 to 0.025 ms: factor 2.5).

Usage: python dt_series.py [NAME]   (default: unitA)
Output: results/<NAME>/dt_series/probe_dt_series.csv (one row per time step, set and speed, with all spike
        times), summary.csv (each halving: conditions whose spike count changed, largest FSL and rate changes;
        conditions differing from the finest step) and margin_60-6.csv (at 0.0125 ms and 60.6 mm/s: the peak
        membrane potential / V_thr after the last spike, to show how close the sets that gain a spike are to threshold)
"""

import sys

import numpy as np
import pandas as pd

import lif

RES = lif.RESULTS
THETA = lif.THETA
STEPS = [0.5, 0.25, 0.125, 0.0625, 0.025, 0.0125, 0.00625, 0.003125, 0.0015625]


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


def peak_after_last_spike(params, dt, k=5):
    """Spike count and peak V/V_thr after the last spike (same RK4 and refractory test as lif.integrate)."""
    data, _ = lif.load_probe_inputs()
    X = params[["beta_mA", "kd", "t_ref_ms", "v_thr_mV"]].to_numpy()
    t, sed = lif.resample_sed(data, k, dt=dt)
    cur = lif.afferent_current(sed[:, None], dt, X[:, 0], X[:, 1], THETA * dt)
    n = len(X)
    v, last, ns, peak = np.zeros(n), -X[:, 2].copy(), np.zeros(n, int), np.zeros(n)
    f = lambda vn, i: (lif.R * cur[i] - vn) / lif.TAU
    i, tn = 0, 0.0
    while tn < t[-1] - 1:
        k1 = f(v, i); k2 = f(v + dt / 2 * k1, i); k3 = f(v + dt / 2 * k2, i); k4 = f(v + dt * k3, i)
        v = np.where(tn >= last + X[:, 2], v + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4), v)
        tn = t[i + 1]
        peak = np.maximum(peak, v / X[:, 3])
        fired = v > X[:, 3]
        v[fired], last[fired], peak[fired] = 0.0, tn, 0.0
        ns[fired] += 1
        i += 1
    return ns, peak


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
    ns, peak = peak_after_last_spike(params, 0.0125)
    n6 = piv["nspk"].xs("60-6", level="speed")
    gain = (n6[0.00625] != n6[0.0125]).to_numpy()
    m = pd.DataFrame(dict(set=params.set, nspk_0_0125=ns, nspk_0_00625=n6[0.00625].to_numpy().astype(int),
                          peak_v_over_vthr_after_last=peak))
    m.to_csv(out / "margin_60-6.csv", index=False)
    print(f"60.6 mm/s at 0.0125 ms: {gain.sum()} sets gain a spike at 0.00625 ms; their peak V/V_thr after the last "
          f"spike is {peak[gain].min():.4f}..{peak[gain].max():.4f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "unitA")
