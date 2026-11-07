"""Spikes during the hold phase of the probe indentation, with and without the SED threshold.

A spike counts as a hold-phase spike when it occurs more than HOLD_MARGIN (10 ms) after the end of the indentation
ramp, taken as 1 mm / speed. The probe (six speeds) is simulated at dt = 0.0125 ms with the spline-resampled SED.

1. The 53 parameter sets of the original submission (param.mat of tag v1.0; their own T_ref) without a threshold
   (theta = 0): the residual SED fluctuations after the indentation produce hold-phase spikes.
2. The 150 sets of the revised version (results/unitA/params.csv; T_ref 2.1 ms) with theta = 0.13: no hold-phase spike.

Requires a checkout of tag v1.0 for part 1 (see README.md).
Usage: python hold_spikes.py
Output: results/v1.0_sets/hold_spikes_theta0.csv and results/unitA/hold_spikes.csv (per set and speed) and a
        printed summary
"""

import numpy as np
import pandas as pd

import lif

SPEEDS = [9.5, 14.8, 20.6, 34.0, 46.5, 60.6]  # mm/s, columns of f61-1.xlsx
RAMP_END_MS = [1000.0 / v for v in SPEEDS]  # 1 mm indentation
HOLD_MARGIN = 10.0  # ms


def count(X, theta, dt=lif.DT):
    """Hold-phase spikes per set (rows of X: beta, K_d, T_ref, V_thr) and speed."""
    data, _ = lif.load_probe_inputs()
    out = np.zeros((len(X), 6), dtype=int)
    for k in range(6):
        t, sed = lif.resample_sed(data, k, dt=dt)
        cur = lif.afferent_current(sed[:, None], dt, X[:, 0], X[:, 1], lif.u_thr(theta, dt))
        for p, s in enumerate(lif.integrate(t, cur, X[:, 2], X[:, 3])):
            out[p, k] = np.sum(s > RAMP_END_MS[k] + HOLD_MARGIN)
    return out


def table(n, label):
    df = pd.DataFrame(n, columns=[f"hold_spikes_{ps}" for ps in lif.PROBE_SPEEDS])
    df.insert(0, "set", np.arange(1, len(n) + 1))
    df["hold_spikes_total"] = n.sum(axis=1)
    print(f"{label}: {n.sum()} hold-phase spikes in {(n.sum(axis=1) > 0).sum()} of {len(n)} sets "
          f"(per speed: {', '.join(f'{ps} {c}' for ps, c in zip(lif.PROBE_SPEEDS, n.sum(axis=0)))})")
    return df


def main():
    v1 = lif.load_v1_params()
    out = lif.RESULTS / "v1.0_sets"
    out.mkdir(parents=True, exist_ok=True)
    table(count(v1, 0.0), "53 v1.0 sets, theta 0, dt 0.0125 ms").to_csv(out / "hold_spikes_theta0.csv", index=False)

    revised = pd.read_csv(lif.RESULTS / "unitA" / "params.csv")[lif.PARAM_COLUMNS].to_numpy()
    table(count(revised, lif.THETA), f"150 revised sets, theta {lif.THETA:g}, dt 0.0125 ms").to_csv(
        lif.RESULTS / "unitA" / "hold_spikes.csv", index=False)


if __name__ == "__main__":
    main()
