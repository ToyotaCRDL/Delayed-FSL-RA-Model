"""Basis of the SED-rate threshold theta = 0.13 J m^-3 ms^-1 and a check that no set fires during the hold.

1. Residual SED fluctuation after the probe indentation. For each of the six speeds, the FE SED (output every 2.5 ms,
   Sheet1 of f61-1.xlsx) is differenced over the period after the probe displacement (Sheet2) reaches its maximum.
   A change counts as a fluctuation when its sign differs from the preceding change (the SED oscillates rather than
   relaxing in one direction). The largest |dU/dt| of these fluctuations over all speeds (0.124) sets theta: the
   smallest value above it at two decimals (0.13). The same maximum is found when the indentation is taken to end at
   the nominal time 1 mm / speed instead of the plateau in Sheet2.
2. Hold-period firing. With the parameter sets of FOLDER, the probe is simulated at DT (spline-resampled SED, as in
   the optimisation). The afferent current is zero once |dU/dt| stays below theta, so a spike after that time would be
   firing during the hold. For each speed the script reports the displacement plateau, the last time |dU/dt| exceeds
   theta, the latest spike over all sets, and how many sets fire after the plateau / after the last suprathreshold step.

Usage: python theta_basis.py [NAME] [THETA] [DT]   (default: unitA 0.13 0.0125)
Output: results/<NAME>/theta_basis.csv (one row per speed) and a printed summary
"""

import sys

import numpy as np
import pandas as pd

import lif

RES = lif.RESULTS
FE_STEP = 2.5  # ms, FE output interval


def residual_fluctuation(sed, t, t_end):
    """Largest |dU/dt| [J m^-3 ms^-1] among sign-alternating changes of the FE SED that start at or after t_end."""
    rate = np.diff(sed) / FE_STEP
    after = t[1:] > t_end  # change over (t[i-1], t[i]] starts at or after t_end
    r, tr = rate[after], t[1:][after]
    alternating = np.r_[False, np.sign(r[1:]) != np.sign(r[:-1])]
    if not alternating.any():
        return 0.0, np.nan
    i = np.flatnonzero(alternating)[np.argmax(np.abs(r[alternating]))]
    return abs(r[i]), tr[i]


def main(folder, theta, dt):
    data, _ = lif.load_probe_inputs()
    disp = lif._read_sheet(lif.PROBE_DIR / "finger61-1" / "f61-1.xlsx", "Sheet2")
    t_fe, t_disp = data[:, 6] * 1000.0, disp[:, 0] * 1000.0
    p = pd.read_csv(RES / folder / "params.csv")
    X = p[["beta_mA", "kd", "t_ref_ms", "v_thr_mV"]].to_numpy()

    rows = []
    for k, ps in enumerate(lif.PROBE_SPEEDS):
        d = disp[:, k + 1]
        plateau = t_disp[np.argmax(d >= d.max() - 1e-9)]
        nominal = d.max() / float(ps.replace("-", ".")) * 1000.0
        fluct, at = residual_fluctuation(data[:, k], t_fe, plateau)
        fluct_nom, _ = residual_fluctuation(data[:, k], t_fe, nominal)

        t, sed = lif.resample_sed(data, k, dt=dt)
        cur = lif.afferent_current(sed[:, None], dt, X[:, 0], X[:, 1], theta * dt)
        last_active = t[np.flatnonzero(cur[:, 0] != 0)].max()
        last_spike = np.array([s[-1] if len(s) else np.nan for s in lif.integrate(t, cur, X[:, 2], X[:, 3])])
        rows.append(dict(speed=ps, displacement_plateau_ms=plateau, nominal_end_ms=nominal,
                         max_fluctuation=fluct, at_ms=at, max_fluctuation_nominal_end=fluct_nom,
                         last_above_theta_ms=last_active, latest_spike_ms=np.nanmax(last_spike),
                         sets_spiking_after_plateau=int(np.sum(last_spike > plateau)),
                         sets_spiking_after_last_above_theta=int(np.sum(last_spike > last_active + dt))))
    df = pd.DataFrame(rows)
    df.to_csv(RES / folder / "theta_basis.csv", index=False)
    print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    m = df.max_fluctuation.max()
    print(f"\nlargest residual fluctuation: {m:.3f} J m^-3 ms^-1 ({df.speed[df.max_fluctuation.idxmax()]} mm/s); "
          f"smallest theta above it at two decimals: {np.floor(m * 100 + 1) / 100:.2f}")
    print(f"{len(p)} sets at theta {theta:g}, dt {dt:g} ms: spikes after the last suprathreshold step: "
          f"{df.sets_spiking_after_last_above_theta.sum()}")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0] if a else "unitA",
         float(a[1]) if len(a) > 1 else 0.13, float(a[2]) if len(a) > 2 else 0.0125)
