"""All analyses of the manuscript for the parameter sets of results/<NAME>/params.csv (revised version).

Same implementation for every stimulus: SED-rate threshold theta (|dU/dt| > theta), time step dt, T_ref from
params.csv; the scanning SED is spline-resampled in time onto dt.

Representative sets (representative.json; keys name the sets of the original submission they replace):
  #1 = best FSL (set 1); #13 = best balance of FSL and firing rate, i.e. the set closest to the ideal point after
  scaling both objectives to [0, 1]; #20, #53 = the sets with dFSL closest to 10 and 30 ms.
The control analysis (imposed unresponsive period of 0-30 ms) uses the set replacing #13.

Noise: the v1.0 scripts add N(0, gamma^2) current noise per 2.5 ms step. To keep the same effect on the membrane
potential at a finer step the standard deviation is scaled by sqrt(2.5 / dt); file names keep the nominal gamma
(noise_2_5, noise_5_0). The random streams differ from MATLAB's.

Usage: python run_all.py [NAME]   (default: unitA)
Output:
  results/<NAME>/summary.csv          per set: dFSL, objectives, probe FSL/rate/spike count, FR/TV per spacing and
                                      their RMSEs, RW and S6 spike width per dot, the same under noise
  results/<NAME>/representative.json
  results/<NAME>/delay_control.csv    control analysis: FR/TV RMSE and RW difference from the experiment per period
  ../LIF_model/*/param.mat, popCost.mat, result1/, result2/ and "statistical analysis"/*.csv in the layout and
  formats of v1.0 (files of the SA1 model and experimental values are left unchanged)
"""

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import scipy.io as sio

import lif
import scan

RW_EXP = np.array([5.13, 4.84, 3.88])  # LaMotte et al. (1998), mm, for the 4.1, 2.6 and 1.6 mm dots
DELAYS = [0, 5, 10, 15, 20, 25, 30]
NOISES = {"2_5": 2.5e-6, "5_0": 5.0e-6}
SPACINGS_ASC = ["1-3", "2-4", "3-2", "4-3", "5-2", "6-2"]  # row order of the FR/TV tables (expFireRate.mat)
RA_FOLDERS = ["stimulus1_probe", "stimulus2_braille", "stimulus3_dot_texture", "stimulus3_dot_texture_delay",
              "stimulus3_dot_texture_noise", "stimulus4_single_dot", "stimulus4_single_dot_delay",
              "stimulus4_single_dot_modi_RW", "stimulus4_single_dot_noise"]
POPCOST_FOLDERS = ["stimulus1_probe", "stimulus3_dot_texture", "stimulus3_dot_texture_delay",
                   "stimulus3_dot_texture_noise"]


def probe_responses(row, theta, dt):
    data, _ = lif.load_probe_inputs()
    out = []
    for k in range(6):
        t, sed = lif.resample_sed(data, k, dt=dt)
        out.append(lif.simulate(t, sed, *row, u_thr=lif.u_thr(theta, dt)))
    return out


def safe_rw(trains, finger):
    """RW, or NaN when it cannot be determined (too few spikes for the 25% crossings)."""
    try:
        return scan.response_width(trains, finger)
    except (IndexError, ValueError):
        return np.nan


def run_set(args):
    n, row, theta, dt = args
    u = lif.u_thr(theta, dt)
    scale = np.sqrt(2.5 / dt)
    res = {"set": n, "probe": probe_responses(row, theta, dt)}
    res["braille"] = scan.run(scan.BRAILLE_DIR, scan.BRAILLE_FINGER, row, u_thr=u, dt=dt)
    res["dot"] = {ds: scan.run(scan.DOT_TEXTURE_DIR, f, row, u_thr=u, dt=dt) for ds, f in scan.DOT_SPACINGS}
    res["single"] = {dw: scan.run(scan.SINGLE_DOT_DIR, f, row, u_thr=u, dt=dt) for dw, f in scan.DOT_WIDTHS}
    res["noise"] = {}
    for label, gamma in NOISES.items():
        rng = np.random.default_rng([n, int(gamma * 1e7)])
        res["noise"][label] = {
            "dot": {ds: scan.run(scan.DOT_TEXTURE_DIR, f, row, u_thr=u, dt=dt, noise=gamma * scale, rng=rng)
                    for ds, f in scan.DOT_SPACINGS},
            "single": {dw: scan.run(scan.SINGLE_DOT_DIR, f, row, u_thr=u, dt=dt, noise=gamma * scale, rng=rng)
                       for dw, f in scan.DOT_WIDTHS},
        }
    return res


def run_delay(args):
    d, row, theta, dt = args
    u = lif.u_thr(theta, dt)
    return d, {
        "dot": {ds: scan.run(scan.DOT_TEXTURE_DIR, f, row, u_thr=u, dt=dt, delay=d) for ds, f in scan.DOT_SPACINGS},
        "single": {dw: scan.run(scan.SINGLE_DOT_DIR, f, row, u_thr=u, dt=dt, delay=d) for dw, f in scan.DOT_WIDTHS},
    }


def dot_metrics(dot):
    """FR and TV per spacing, ordered 1.3 ... 6.2 mm like expFireRate."""
    fr = np.zeros(6)
    tv = np.zeros(6)
    for mlp, (ds, finger) in enumerate(scan.DOT_SPACINGS):
        fr[5 - mlp], tv[5 - mlp] = scan.dot_texture_metrics(dot[ds], finger)
    return fr, tv


def single_metrics(single):
    rw = np.array([safe_rw(single[dw], f) for dw, f in scan.DOT_WIDTHS])
    sw = np.array([scan.spike_width(single[dw]) if sum(map(len, single[dw])) else np.nan for dw, _ in scan.DOT_WIDTHS])
    return rw, sw


def representative_sets(params, summary):
    f = params[["J_fsl", "J_rate"]].to_numpy()
    z = (f - f.min(0)) / (f.max(0) - f.min(0))
    knee = int(params.set.iloc[int(np.argmin(np.hypot(z[:, 0], z[:, 1])))])
    return {"published_#1": 1, "published_#13": knee,
            "published_#20": int((summary.dFSL - 10.0).abs().idxmin()),
            "published_#53": int((summary.dFSL - 30.0).abs().idxmin()),
            "rule": "#1 = best FSL (set 1); #13 = best balance of FSL and rate (closest to the ideal point with both "
                    "objectives scaled to [0, 1]); #20, #53 = sets with ΔFSL closest to 10 and 30 ms"}


# ------------------------------------------------------------------ writers (v1.0 CSV formats)

def write_sep(path, trains, placeholder):
    """One row per afferent, zero-padded; `placeholder`: leading 0 column (spike.m's t_spike(1)) as in v1.0."""
    width = max([len(s) for s in trains] + [1]) + (1 if placeholder else 0)
    mat = np.zeros((len(trains), width))
    for i, s in enumerate(trains):
        start = 1 if placeholder else 0
        mat[i, start:start + len(s)] = s
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, mat, delimiter=",", fmt="%.10g")


def write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, np.atleast_2d(rows), delimiter=",", fmt="%.10g")


def write_ycoord(path, stim_dir, finger):
    _, _, y = scan.centre_line_sed(stim_dir, finger)
    write_rows(path, y[:, None])


def clear_ra_outputs(folder):
    """Remove the RA output files of a previous parameter-set series (e.g. the 53 sets of v1.0)."""
    for sub in ("result1", "result2"):
        for f in (folder / sub).glob("*_param_*.csv"):
            f.unlink()


def write_lif_model(params, summary, results, delay_res, rep13):
    L = lif.LIF_DIR
    n_sets = len(params)
    X = params[lif.PARAM_COLUMNS].to_numpy()
    J = params[["J_fsl", "J_rate"]].to_numpy()
    for name in RA_FOLDERS:
        clear_ra_outputs(L / name)
        sio.savemat(L / name / "param.mat", {"param": X})
    for name in POPCOST_FOLDERS:
        sio.savemat(L / name / "popCost.mat", {"popCost": J})

    order = summary.sort_values("dFSL").index
    d_sorted = summary.dFSL.loc[order].to_numpy()
    exp_fr, exp_tv = scan.load_dot_texture_experiment()

    # Stimulus 1: probe
    probe = L / "stimulus1_probe"
    write_rows(probe / "result1" / "probe_pareto_sol.csv", J)
    for res in results:
        for k, ps in enumerate(lif.PROBE_SPEEDS):
            write_rows(probe / "result2" / f"probe_spiketrain_param_{res['set']}_ps_{ps}.csv",
                       np.concatenate([[0.0], res["probe"][k]])[:, None])

    # Stimulus 2: Braille
    braille = L / "stimulus2_braille" / "result1"
    for res in results:
        write_sep(braille / f"braille_sep_param_{res['set']}.csv", res["braille"], True)
    write_ycoord(braille / "braille_sep_ycoord.csv", scan.BRAILLE_DIR, scan.BRAILLE_FINGER)

    # Stimulus 3: dot pattern (and noise)
    dot = L / "stimulus3_dot_texture"
    dot_noise = L / "stimulus3_dot_texture_noise"
    for res in results:
        for ds, _ in scan.DOT_SPACINGS:
            write_sep(dot / "result1" / f"dot_texture_sep_param_{res['set']}_ds_{ds}.csv", res["dot"][ds], True)
            for label in NOISES:
                write_sep(dot_noise / "result1" / f"dot_texture_sep_param_{res['set']}_ds_{ds}_noise_{label}.csv",
                          res["noise"][label]["dot"][ds], True)
    for folder in (dot, dot_noise, L / "stimulus3_dot_texture_delay"):
        write_ycoord(folder / "result1" / "dot_texture_sep_ycoord.csv", scan.DOT_TEXTURE_DIR, scan.DOT_SPACINGS[0][1])
    for name, prefix, exp in (("firerate", "fr", exp_fr), ("tempvar", "tv", exp_tv)):
        sim = summary[[f"{prefix}_{s}" for s in SPACINGS_ASC]].to_numpy().T  # 6 x sets
        diff = exp[:, None] - sim
        rmse = summary[f"{prefix}_rmse"]
        write_rows(dot / "result2" / f"dot_texture_{name}_sim.csv", sim)
        write_rows(dot / "result2" / f"dot_texture_{name}_diff.csv", diff)
        write_rows(dot / "result2" / f"dot_texture_{name}_diff2.csv", np.vstack([d_sorted, diff[:, order - 1]]))
        write_rows(dot / "result2" / f"dot_texture_{name}_rmse.csv", rmse.to_numpy()[:, None])
        rmse2 = np.column_stack([d_sorted, rmse.loc[order]])
        write_rows(dot / "result2" / f"dot_texture_{name}_rmse2.csv", rmse2)
        write_rows(L / "statistical analysis" / f"dot_texture_{name}_rmse2.csv", rmse2)
        for label in NOISES:
            col = summary[f"{prefix}_rmse_noise_{label}"]
            write_rows(dot_noise / "result2" / f"dot_texture_{name}_rmse_noise_{label}.csv", col.to_numpy()[:, None])
            write_rows(dot_noise / "result2" / f"dot_texture_{name}_rmse2_noise_{label}.csv",
                       np.column_stack([d_sorted, col.loc[order]]))

    # Stimulus 4: single dot (and noise, S6 spike width)
    single = L / "stimulus4_single_dot"
    single_noise = L / "stimulus4_single_dot_noise"
    for res in results:
        for dw, _ in scan.DOT_WIDTHS:
            write_sep(single / "result1" / f"single_dot_sep_param_{res['set']}_dw_{dw}.csv", res["single"][dw], False)
            for label in NOISES:
                write_sep(single_noise / "result1" / f"single_dot_sep_param_{res['set']}_dw_{dw}_noise_{label}.csv",
                          res["noise"][label]["single"][dw], False)
    for folder in (single, single_noise, L / "stimulus4_single_dot_delay"):
        write_ycoord(folder / "result1" / "single_dot_sep_ycoord.csv", scan.SINGLE_DOT_DIR, scan.DOT_WIDTHS[0][1])
    dws = [dw for dw, _ in scan.DOT_WIDTHS]
    rw = summary[[f"rw_{dw}" for dw in dws]].to_numpy()
    write_rows(single / "result2" / "single_dot_respwidth.csv", rw)
    write_rows(single / "result2" / "single_dot_respwidth2.csv", np.column_stack([d_sorted, (rw - RW_EXP)[order - 1]]))
    for j, dw in enumerate(dws):
        col = np.column_stack([d_sorted, (rw - RW_EXP)[order - 1, j]])
        write_rows(L / "statistical analysis" / f"single_dot_respwidth2_dw_{dw.replace('-', '_')}.csv",
                   col[~np.isnan(col[:, 1])])
    sw = summary[[f"sw_{dw}" for dw in dws]].to_numpy()
    modi = L / "stimulus4_single_dot_modi_RW" / "result2"
    write_rows(modi / "single_dot_respwidth.csv", sw)
    write_rows(modi / "single_dot_respwidth2.csv", np.column_stack([d_sorted, sw[order - 1]]))
    for label in NOISES:
        nrw = summary[[f"rw_{dw}_noise_{label}" for dw in dws]].to_numpy()
        write_rows(single_noise / "result2" / f"single_dot_respwidth_noise_{label}.csv", nrw)
        write_rows(single_noise / "result2" / f"single_dot_respwidth2_noise_{label}.csv",
                   np.column_stack([d_sorted, (nrw - RW_EXP)[order - 1]]))

    # control analysis (imposed unresponsive period), set rep13
    dot_delay = L / "stimulus3_dot_texture_delay"
    single_delay = L / "stimulus4_single_dot_delay"
    rows = []
    for d in DELAYS:
        r = delay_res[d]
        for ds, _ in scan.DOT_SPACINGS:
            write_sep(dot_delay / "result1" / f"dot_texture_sep_param_{rep13}_ds_{ds}_delay_{d}.csv", r["dot"][ds], True)
        for dw in dws:
            write_sep(single_delay / "result1" / f"single_dot_sep_param_{rep13}_dw_{dw}_delay_{d}.csv",
                      r["single"][dw], False)
        fr, tv = dot_metrics(r["dot"])
        rw_d, _ = single_metrics(r["single"])
        rows.append((d, scan.rmse(exp_fr, fr), scan.rmse(exp_tv, tv), *rw_d))
    rows = np.array(rows)
    write_rows(dot_delay / "result2" / "dot_texture_firerate_rmse.csv", rows[:, [0, 1]])
    write_rows(dot_delay / "result2" / "dot_texture_tempvar_rmse.csv", rows[:, [0, 2]])
    write_rows(single_delay / "result2" / "single_dot_respwidth.csv", rows[:, [0, 3, 4, 5]])
    write_rows(single_delay / "result2" / "single_dot_respwidth2.csv",
               np.column_stack([rows[:, 0], rows[:, 3:6] - RW_EXP]))
    return rows


def main(name="unitA", theta=lif.THETA, dt=lif.DT):
    out = lif.RESULTS / name
    params = pd.read_csv(out / "params.csv")
    rows = params[lif.PARAM_COLUMNS].to_numpy()

    with ProcessPoolExecutor(max_workers=os.cpu_count()) as pool:
        results = list(pool.map(run_set, [(n, r, theta, dt) for n, r in enumerate(rows, start=1)]))

    _, actual = lif.load_probe_inputs()
    exp_trains = [c[~np.isnan(c) & (c > 0)] for c in actual.T]
    fsl_obs = np.array([s[0] for s in exp_trains])
    exp_fr, exp_tv = scan.load_dot_texture_experiment()

    records = []
    for res in results:
        n = res["set"]
        fsl = np.array([s[0] if len(s) else 300.0 for s in res["probe"]])
        rate = np.array([1000 / np.mean(np.diff(s)) if len(s) > 1 else np.nan for s in res["probe"]])
        fr, tv = dot_metrics(res["dot"])
        rw, sw = single_metrics(res["single"])
        rec = {"set": n, "dFSL": fsl[0] - fsl_obs[0], "J_fsl": params.J_fsl[n - 1], "J_rate": params.J_rate[n - 1],
               "fr_rmse": scan.rmse(exp_fr, fr), "tv_rmse": scan.rmse(exp_tv, tv)}
        for k, ps in enumerate(lif.PROBE_SPEEDS):
            rec[f"fsl_{ps}"], rec[f"rate_{ps}"], rec[f"nspk_{ps}"] = fsl[k], rate[k], len(res["probe"][k])
        for i, ds in enumerate(SPACINGS_ASC):
            rec[f"fr_{ds}"], rec[f"tv_{ds}"] = fr[i], tv[i]
        for j, (dw, _) in enumerate(scan.DOT_WIDTHS):
            rec[f"rw_{dw}"], rec[f"sw_{dw}"] = rw[j], sw[j]
        for label in NOISES:
            nfr, ntv = dot_metrics(res["noise"][label]["dot"])
            nrw, _ = single_metrics(res["noise"][label]["single"])
            rec[f"fr_rmse_noise_{label}"], rec[f"tv_rmse_noise_{label}"] = scan.rmse(exp_fr, nfr), scan.rmse(exp_tv, ntv)
            for j, (dw, _) in enumerate(scan.DOT_WIDTHS):
                rec[f"rw_{dw}_noise_{label}"] = nrw[j]
        records.append(rec)
    summary = pd.DataFrame(records).sort_values("set").set_index("set")
    summary.to_csv(out / "summary.csv")

    reps = representative_sets(params, summary)
    (out / "representative.json").write_text(json.dumps(reps, indent=2, ensure_ascii=False))
    rep13 = reps["published_#13"]

    with ProcessPoolExecutor(max_workers=len(DELAYS)) as pool:
        delay_res = dict(pool.map(run_delay, [(d, rows[rep13 - 1], theta, dt) for d in DELAYS]))

    delay_rows = write_lif_model(params, summary, results, delay_res, rep13)
    pd.DataFrame(np.column_stack([delay_rows[:, :3], delay_rows[:, 3:6] - RW_EXP]),
                 columns=["delay_ms", "fr_rmse", "tv_rmse", "dRW_4-1", "dRW_2-6", "dRW_1-6"]).to_csv(
        out / "delay_control.csv", index=False)
    print(f"{name}: {len(summary)} sets; representative sets {reps}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "unitA")
