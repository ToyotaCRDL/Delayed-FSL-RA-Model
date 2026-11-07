"""Fronts optimised with equal and with relative-error weights across the probe speeds, compared.

Both fronts use the revised implementation (theta 0.13, T_ref 2.1 ms, dt 0.0125 ms, NSGA-II, same seed).
Equal weights: J = -(1 - sum (obs - sim)^2 / sum obs^2) (manuscript). Relative weights: each speed's squared error
divided by obs^2, i.e. J = -(1 - mean squared relative error) (optimize.py --relative).

Prints, for each front: ranges of dFSL and RMS rate error, per-speed shares of each objective (under the weighting
that front was optimised with), Spearman rho of the SEP metrics against dFSL, the dFSL at which the 1.3 mm response
disappears, and the representative sets (knee, dFSL closest to 10 and 30 ms). Also the difference in SEP metrics
between the fronts at matched dFSL (equal-weight front interpolated at the dFSL of each relative-weight set).

Requires results/unitA/front_summary.csv and results/unitA_relative_weights/front_summary.csv (evaluate_front.py).
Usage: python relative_weights.py
Output: results/unitA_relative_weights/comparison.csv (per set, both fronts) and stdout
"""

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import lif
from optimize import ProbeObjectives

FRONTS = {"equal": "unitA", "relative": "unitA_relative_weights"}
SPEEDS = ["9-5", "14-8", "20-6", "34-0", "46-5", "60-6"]
SPACINGS = ["2-4", "3-2", "4-3", "5-2", "6-2"]
METRICS = ["fr_rmse", "tv_rmse", "rw_4-1", "rw_2-6", "rw_1-6"]


def observed():
    p = ProbeObjectives(0.13, 0.0125)
    return p.fsl_obs, 1000 * p.rate_obs


def knee(df):
    f = df[["J_fsl", "J_rate"]].to_numpy()
    f = (f - f.min(0)) / (f.max(0) - f.min(0))
    return int(df.set.iloc[np.argmin(np.hypot(f[:, 0], f[:, 1]))])


def shares(df, obs, prefix, relative):
    sim = df[[f"{prefix}_{s}" for s in SPEEDS]].to_numpy()
    sq = (sim - obs) ** 2 / (obs**2 if relative else 1)
    return sq / sq.sum(axis=1, keepdims=True)


def main():
    fsl_obs, rate_obs = observed()
    out = []
    for name, folder in FRONTS.items():
        df = pd.read_csv(lif.RESULTS / folder / "front_summary.csv")
        df["front"] = name
        rel = name == "relative"
        print(f"\n=== {name} weights ({folder}) ===")
        print(f"dFSL {df.dFSL.min():.1f}..{df.dFSL.max():.1f} ms; sets >=20 / >=30 ms: "
              f"{(df.dFSL >= 20).sum()} / {(df.dFSL >= 30).sum()}; "
              f"RMS rate error {df.rms_rate_err.min():.1f}..{df.rms_rate_err.max():.1f} spikes/s")
        sf, sr = shares(df, fsl_obs, "fsl", rel), shares(df, rate_obs, "rate", rel)
        late = df.dFSL.to_numpy() >= 5
        top = np.argmax(sf, axis=1)
        print(f"J_FSL: sets with dFSL >= 5 ms: {late.sum()}; largest share at 9.5 mm/s in {(top[late] == 0).sum()} "
              f"(share {sf[late, 0].min():.0%}..{sf[late, 0].max():.0%})")
        print("J_FSL median share by speed: " + ", ".join(f"{s} {m:.0%}" for s, m in zip(SPEEDS, np.median(sf, 0))))
        print("J_Rate median share by speed: " + ", ".join(f"{s} {m:.0%}" for s, m in zip(SPEEDS, np.median(sr, 0))))
        big = sr.max(axis=1) > 0.5
        print(f"J_Rate: sets with one speed > 50%: {big.sum()} "
              f"({dict(zip(*np.unique(np.array(SPEEDS)[np.argmax(sr[big], 1)], return_counts=True)))})")
        err = df[[f"rate_{s}" for s in SPEEDS]].to_numpy() - rate_obs
        print("rate error median by speed (spikes/s): " + ", ".join(f"{s} {m:.0f}" for s, m in zip(SPEEDS, np.median(err, 0))))
        print("Spearman rho vs dFSL: " + ", ".join(f"{m} {spearmanr(df.dFSL, df[m], nan_policy='omit')[0]:.2f}" for m in METRICS))
        on, off = df[df["fr_1-3"] > 0], df[df["fr_1-3"] == 0]
        print(f"1.3 mm response: present up to dFSL {on.dFSL.max():.1f} ms, absent from {off.dFSL.min():.1f} ms")
        reps = {"knee": knee(df),
                "10 ms": int(df.set.iloc[np.argmin(np.abs(df.dFSL - 10))]),
                "30 ms": int(df.set.iloc[np.argmin(np.abs(df.dFSL - 30))])}
        for label, s in reps.items():
            r = df[df.set == s].iloc[0]
            print(f"  {label}: #{s} dFSL {r.dFSL:.1f} ms, rate err {r.rms_rate_err:.1f}, FR {r.fr_rmse:.1f}, "
                  f"TV {r.tv_rmse:.1f}, RW {r['rw_4-1']:.2f}/{r['rw_2-6']:.2f}/{r['rw_1-6']:.2f}, "
                  f"responders 2.4-6.2 mm {'/'.join(str(int(r['resp_' + d])) for d in SPACINGS)}, "
                  f"single dot {'/'.join(str(int(r['resp_single_' + d])) for d in ['4-1', '2-6', '1-6'])}")
        out.append(df)

    eq, rw = out
    cols = METRICS + ["rms_rate_err"] + [f"resp_{d}" for d in SPACINGS + ["1-3"]] + ["resp_single_1-6"]
    inside = rw[(rw.dFSL >= eq.dFSL.min()) & (rw.dFSL <= eq.dFSL.max())]
    print(f"\n=== relative minus equal at matched dFSL ({len(inside)} relative-weight sets inside the equal range) ===")
    for m in cols:
        ref = eq[["dFSL", m]].dropna().groupby("dFSL", as_index=False)[m].mean().sort_values("dFSL")
        x = inside[["dFSL", m]].dropna()
        d = x[m].to_numpy() - np.interp(x.dFSL, ref.dFSL, ref[m])
        print(f"{m}: median {np.median(d):+.2f}, 5-95% {np.percentile(d, 5):+.2f}..{np.percentile(d, 95):+.2f} (n={len(d)})")

    dest = lif.RESULTS / "unitA_relative_weights"
    pd.concat(out).to_csv(dest / "comparison.csv", index=False)
    print(f"-> {dest / 'comparison.csv'}")


if __name__ == "__main__":
    main()
