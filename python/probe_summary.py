"""Probe responses of every set of a front: FSL and firing rate at each of the six speeds, dFSL, RMS rate error.

Usage: python probe_summary.py [NAME] [--data WORKBOOK]   (default NAME: unitA)
  Unit B (Fig. S1):  python probe_summary.py unitB --data ../LIF_model/stimulus1_probe/finger62-1/f62-1_new.xlsx
Output: results/<NAME>/probe_summary.csv
  dFSL: simulated minus observed FSL at the slowest speed (first column of the workbook), ms
  dFSL_max_abs: largest |simulated - observed FSL| over the six speeds, ms
  rms_rate_err: RMS over the six speeds of the firing-rate error, spikes/s
  fsl_k, rate_k: FSL [ms] and 1/mean ISI [spikes/s] at the k-th speed
"""

import sys

import numpy as np
import pandas as pd

import lif
from optimize import ProbeObjectives


def main(name, theta=lif.THETA, dt=lif.DT, data_file=None):
    folder = lif.RESULTS / name
    p = pd.read_csv(folder / "params.csv")
    obj = ProbeObjectives(theta, dt, data_file)
    X = p[lif.PARAM_COLUMNS].to_numpy()
    fsl, rate = obj.responses(X)
    err = fsl - obj.fsl_obs
    out = p.copy()
    out["dFSL"] = err[:, 0]
    out["dFSL_max_abs"] = np.abs(err).max(axis=1)
    out["rms_rate_err"] = 1000 * np.sqrt(np.nanmean((rate - obj.rate_obs) ** 2, axis=1))
    for k in range(6):
        out[f"fsl_{k + 1}"] = fsl[:, k]
        out[f"rate_{k + 1}"] = 1000 * rate[:, k]
    out.to_csv(folder / "probe_summary.csv", index=False)
    print(f"observed FSL (ms): {np.round(obj.fsl_obs, 1)}; rate (spikes/s): {np.round(1000 * obj.rate_obs, 1)}")
    print(f"{len(out)} sets: dFSL {out.dFSL.min():.1f}..{out.dFSL.max():.1f} ms, "
          f"max |FSL error| {out.dFSL_max_abs.min():.1f}..{out.dFSL_max_abs.max():.1f} ms, "
          f"RMS rate error {out.rms_rate_err.min():.1f}..{out.rms_rate_err.max():.1f} spikes/s")


if __name__ == "__main__":
    a = sys.argv[1:]
    data = None
    if "--data" in a:
        i = a.index("--data")
        data = a[i + 1]
        del a[i:i + 2]
    main(a[0] if a else "unitA", data_file=data)
