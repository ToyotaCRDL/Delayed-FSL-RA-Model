"""Check scan.py against the v1.0 scanning results: SEP spike trains, FR/TV RMSE (Stim. 3), RW and S6 width (Stim. 4).

Runs the 53 v1.0 sets with the v1.0 settings (u_thr 1.1 per step, FE time step 2.5 ms). Requires a checkout of tag
v1.0; its MATLAB outputs (result1/*.csv, set*.mat, spikeWidth.mat) are the reference.
"""

import numpy as np
import scipy.io as sio

import lif
import scan

V1 = lif.V1_LIF_DIR


def load_released_sep(path, drop_placeholder):
    rows = np.loadtxt(path, delimiter=",", ndmin=2)
    if drop_placeholder:
        rows = rows[:, 1:]
    return [r[r > 0] for r in rows]


def same_trains(a, b):
    return len(a) == len(b) and all(len(x) == len(y) and np.allclose(x, y, atol=1e-6) for x, y in zip(a, b))


def main():
    params = lif.load_v1_params()
    exp_fr, exp_tv = scan.load_dot_texture_experiment()
    rel_fr = sio.loadmat(V1 / scan.DOT_TEXTURE_DIR.name / "setFiringRateRMSE.mat")["setFiringRateRMSE"].ravel()
    rel_tv = sio.loadmat(V1 / scan.DOT_TEXTURE_DIR.name / "setTempoVariRMSE.mat")["setTempoVariRMSE"].ravel()
    rel_rw = sio.loadmat(V1 / scan.SINGLE_DOT_DIR.name / "setRespWidth.mat")["setRespWidth"]
    rel_sw = sio.loadmat(V1 / "stimulus4_single_dot_modi_RW" / "spikeWidth.mat")["spikeWidth"]

    sep_ok = sep_bad = 0
    fr_err = []
    tv_err = []
    rw_err = []
    sw_err = []
    for n, row in enumerate(params, start=1):
        fr = np.zeros(6)
        tv = np.zeros(6)
        for mlp, (ds, finger) in enumerate(scan.DOT_SPACINGS):
            trains = scan.run(scan.DOT_TEXTURE_DIR, finger, row)
            ref = load_released_sep(V1 / scan.DOT_TEXTURE_DIR.name / "result1" / f"dot_texture_sep_param_{n}_ds_{ds}.csv", True)
            if same_trains(trains, ref):
                sep_ok += 1
            else:
                sep_bad += 1
                print(f"SEP mismatch: #{n} dot spacing {ds}")
            fr[5 - mlp], tv[5 - mlp] = scan.dot_texture_metrics(trains, finger, histogram=scan.firing_rate_histogram_v1)
        fr_err.append(scan.rmse(exp_fr, fr) - rel_fr[n - 1])
        tv_err.append(scan.rmse(exp_tv, tv) - rel_tv[n - 1])

        for mlp, (dw, finger) in enumerate(scan.DOT_WIDTHS):
            trains = scan.run(scan.SINGLE_DOT_DIR, finger, row)
            ref = load_released_sep(V1 / scan.SINGLE_DOT_DIR.name / "result1" / f"single_dot_sep_param_{n}_dw_{dw}.csv", False)
            if same_trains(trains, ref):
                sep_ok += 1
            else:
                sep_bad += 1
                print(f"SEP mismatch: #{n} dot width {dw}")
            if not (n > 47 and mlp == 0):
                rw_err.append(scan.response_width(trains, finger) - rel_rw[n - 1, mlp])
            sw_err.append(scan.spike_width(trains) - rel_sw[n - 1, mlp])

    print(f"SEP spike trains: match={sep_ok} mismatch={sep_bad}")
    for name, err in [("FR RMSE", fr_err), ("TV RMSE", tv_err), ("RW", rw_err), ("spike width (S6)", sw_err)]:
        err = np.abs(np.array(err))
        print(f"{name}: max |python - released| = {np.nanmax(err):.3g} over {len(err)} values")


if __name__ == "__main__":
    main()
