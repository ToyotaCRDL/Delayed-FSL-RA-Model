"""Check the Braille, imposed-delay and noise code of scan.py against the v1.0 results (v1.0 settings and sets).

Requires a checkout of tag v1.0.

Noise uses MATLAB's randn stream, which cannot be reproduced; its RMSEs are compared only in distribution.
"""

import numpy as np
import scipy.io as sio

import lif
import scan
from verify_v1_scan import load_released_sep, same_trains

DELAYS = [0, 5, 10, 15, 20, 25, 30]
V1 = lif.V1_LIF_DIR
DOT_DELAY_DIR = V1 / "stimulus3_dot_texture_delay"
SINGLE_DELAY_DIR = V1 / "stimulus4_single_dot_delay"
DOT_NOISE_DIR = V1 / "stimulus3_dot_texture_noise"


def check_braille(params):
    ok = bad = 0
    for n, row in enumerate(params, start=1):
        trains = scan.run(scan.BRAILLE_DIR, scan.BRAILLE_FINGER, row)
        ref = load_released_sep(V1 / scan.BRAILLE_DIR.name / "result1" / f"braille_sep_param_{n}.csv", True)
        if same_trains(trains, ref):
            ok += 1
        else:
            bad += 1
            print(f"  Braille mismatch: #{n}")
    print(f"Braille SEP: match={ok} mismatch={bad}")


def check_delay(params):
    row = params[12]  # the control analysis fixes Parameter #13
    exp_fr, exp_tv = scan.load_dot_texture_experiment()
    rel_fr = sio.loadmat(DOT_DELAY_DIR / "setFiringRateRMSE.mat")["setFiringRateRMSE"].ravel()
    rel_tv = sio.loadmat(DOT_DELAY_DIR / "setTempoVariRMSE.mat")["setTempoVariRMSE"].ravel()
    rel_rw = sio.loadmat(SINGLE_DELAY_DIR / "setRespWidth.mat")["setRespWidth"]
    ok = bad = 0
    for k, d in enumerate(DELAYS):
        fr = np.zeros(6)
        tv = np.zeros(6)
        for mlp, (ds, finger) in enumerate(scan.DOT_SPACINGS):
            trains = scan.run(scan.DOT_TEXTURE_DIR, finger, row, delay=d)
            ref = load_released_sep(DOT_DELAY_DIR / "result1" / f"dot_texture_sep_param_13_ds_{ds}_delay_{d}.csv", True)
            ok, bad = (ok + 1, bad) if same_trains(trains, ref) else (ok, bad + 1)
            fr[5 - mlp], tv[5 - mlp] = scan.dot_texture_metrics(trains, finger)
        rw = []
        for mlp, (dw, finger) in enumerate(scan.DOT_WIDTHS):
            trains = scan.run(scan.SINGLE_DOT_DIR, finger, row, delay=d)
            ref = load_released_sep(SINGLE_DELAY_DIR / "result1" / f"single_dot_sep_param_13_dw_{dw}_delay_{d}.csv", False)
            ok, bad = (ok + 1, bad) if same_trains(trains, ref) else (ok, bad + 1)
            rw.append(scan.response_width(trains, finger))
        print(f"  delay {d:2d} ms: FR RMSE {scan.rmse(exp_fr, fr):.4f} (released {rel_fr[k]:.4f}), "
              f"TV RMSE {scan.rmse(exp_tv, tv):.4f} (released {rel_tv[k]:.4f}), "
              f"max |RW - released| {np.max(np.abs(np.array(rw) - rel_rw[k])):.2e}")
    print(f"Delay SEP: match={ok} mismatch={bad}")


def check_noise(params, gamma=2.5e-6, sets=(1, 13, 30, 53)):
    exp_fr, exp_tv = scan.load_dot_texture_experiment()
    rel_fr = sio.loadmat(DOT_NOISE_DIR / "setFiringRateRMSE_noise_2_5.mat")["setFiringRateRMSE"].ravel()
    rel_tv = sio.loadmat(DOT_NOISE_DIR / "setTempoVariRMSE_noise_2_5.mat")["setTempoVariRMSE"].ravel()
    rng = np.random.default_rng(0)
    for n in sets:
        vals = []
        for _ in range(5):
            fr = np.zeros(6)
            tv = np.zeros(6)
            for mlp, (ds, finger) in enumerate(scan.DOT_SPACINGS):
                trains = scan.run(scan.DOT_TEXTURE_DIR, finger, params[n - 1], noise=gamma, rng=rng)
                fr[5 - mlp], tv[5 - mlp] = scan.dot_texture_metrics(trains, finger)
            vals.append((scan.rmse(exp_fr, fr), scan.rmse(exp_tv, tv)))
        v = np.array(vals)
        print(f"  noise #{n}: FR RMSE {v[:, 0].mean():.2f}±{v[:, 0].std():.2f} (released {rel_fr[n - 1]:.2f}), "
              f"TV RMSE {v[:, 1].mean():.2f}±{v[:, 1].std():.2f} (released {rel_tv[n - 1]:.2f})")


def main():
    params = lif.load_v1_params()
    check_braille(params)
    check_delay(params)
    check_noise(params)


if __name__ == "__main__":
    main()
