"""Check that lif.py reproduces the v1.0 probe spike trains (result2/*.csv of tag v1.0) for every set and speed.

Runs the 53 v1.0 sets with the v1.0 settings (u_thr 1.1 per step, dt 0.5 ms). Requires a checkout of tag v1.0.
"""

import numpy as np

import lif

RESULT2 = lif.V1_LIF_DIR / "stimulus1_probe" / "result2"


def load_released(n, ps):
    path = RESULT2 / f"probe_spiketrain_param_{n}_ps_{ps}.csv"
    if not path.exists():
        return None
    vals = np.loadtxt(path, delimiter=",", ndmin=1)
    return vals[vals > 0]


def main():
    data, _ = lif.load_probe_inputs()
    params = lif.load_v1_params()
    n_ok = n_bad = n_missing = 0
    for n, row in enumerate(params, start=1):
        for k, ps in enumerate(lif.PROBE_SPEEDS):
            ref = load_released(n, ps)
            if ref is None:
                n_missing += 1
                continue
            t, sed = lif.resample_sed(data, k)
            sim = lif.simulate(t, sed, *row)
            if len(sim) == len(ref) and np.allclose(sim, ref, atol=1e-6):
                n_ok += 1
            else:
                n_bad += 1
                print(f"#{n} {ps} mm/s: python={sim.tolist()} released={ref.tolist()}")
    print(f"match={n_ok} mismatch={n_bad} missing={n_missing}")


if __name__ == "__main__":
    main()
