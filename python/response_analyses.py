"""Analyses for the response to reviewers, for the revised parameter sets.

Reads results/<NAME>/params.csv and representative.json (run_all.py); same implementation for every stimulus.

Spike-count-matched random thinning: the spike trains of the set replacing #13 are thinned to the total spike count
  of (a) the control analysis (input set to zero for 5-30 ms from each onset, same set) and (b) the sets replacing
  #20 and #53, per stimulus condition (6 dot spacings, 3 dot widths), N_REP times, by keeping that many spikes chosen
  uniformly at random from all spikes of the 51 afferents. In the CSVs, `target` holds the number of the original set
  that a revised set replaces (13, 20, 53) and `set` the revised set.
Spike counts in the 0-100 ms window from stimulus onset: experiment and the four representative sets.
Per-speed error and share of each objective (squared error at that speed / sum over speeds) for every set.

Usage: python response_analyses.py [NAME]   (default: unitA)
Output: results/<NAME>/thinning/{thinning_all,thinning_pareto_all}.csv, spike_counts/spike_counts_0_100ms.csv,
        speed_contribution/per_set_per_speed.csv
"""

import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import lif
import scan

SPEEDS = [9.5, 14.8, 20.6, 34.0, 46.5, 60.6]
DELAYS = (0, 5, 10, 15, 20, 25, 30)
N_REP = 200


def thin(trains, n_keep, rng):
    """Keep n_keep spikes chosen uniformly at random from all spikes of the 51 afferents."""
    owner = np.concatenate([np.full(len(s), a) for a, s in enumerate(trains)])
    times = np.concatenate(trains) if len(owner) else np.array([])
    keep = np.zeros(len(owner), dtype=bool)
    keep[rng.choice(len(owner), size=n_keep, replace=False)] = True
    return [np.sort(times[keep & (owner == a)]) for a in range(len(trains))]


def evaluate(dot_trains, single_trains, exp_fr, exp_tv):
    fr = np.zeros(6)
    tv = np.zeros(6)
    rec = {}
    for mlp, (ds, finger) in enumerate(scan.DOT_SPACINGS):
        tr = dot_trains[ds]
        fr[5 - mlp], tv[5 - mlp] = scan.dot_texture_metrics(tr, finger)
        rec[f"resp_{ds}"] = sum(len(s) > 0 for s in tr)
        rec[f"nspk_{ds}"] = sum(len(s) for s in tr)
    rec["fr_rmse"], rec["tv_rmse"] = scan.rmse(exp_fr, fr), scan.rmse(exp_tv, tv)
    for dw, finger in scan.DOT_WIDTHS:
        tr = single_trains[dw]
        try:
            rec[f"rw_{dw}"] = scan.response_width(tr, finger)
        except (IndexError, ValueError):
            rec[f"rw_{dw}"] = np.nan
        rec[f"resp_single_{dw}"] = sum(len(s) > 0 for s in tr)
        rec[f"nspk_single_{dw}"] = sum(len(s) for s in tr)
    return rec


def scan_set(args):
    row, theta, dt, delay = args
    u = theta * dt
    dot = {ds: scan.run(scan.DOT_TEXTURE_DIR, f, row, u_thr=u, dt=dt, delay=delay) for ds, f in scan.DOT_SPACINGS}
    single = {dw: scan.run(scan.SINGLE_DOT_DIR, f, row, u_thr=u, dt=dt, delay=delay) for dw, f in scan.DOT_WIDTHS}
    return dot, single


def thin_job(args):
    key, rep, base, counts, seed = args
    rng = np.random.default_rng(seed)
    exp_fr, exp_tv = scan.load_dot_texture_experiment()
    dot = {ds: thin(base[0][ds], counts[f"nspk_{ds}"], rng) for ds, _ in scan.DOT_SPACINGS}
    single = {dw: thin(base[1][dw], counts[f"nspk_single_{dw}"], rng) for dw, _ in scan.DOT_WIDTHS}
    return key, rep, evaluate(dot, single, exp_fr, exp_tv)


def probe_trains(row, theta, dt):
    data, _ = lif.load_probe_inputs()
    out = []
    for k in range(6):
        t, sed = lif.resample_sed(data, k, dt=dt)
        out.append(lif.simulate(t, sed, *row, u_thr=theta * dt))
    return out


def probe_job(args):
    n, row, theta, dt = args
    return n, probe_trains(row, theta, dt)


def main(name="unitA", theta=lif.THETA, dt=lif.DT):
    out = lif.RESULTS / name
    params = pd.read_csv(out / "params.csv").set_index("set")
    reps = json.loads((out / "representative.json").read_text())
    rows = {n: tuple(params.loc[n, lif.PARAM_COLUMNS]) for n in params.index}
    stand = {13: reps["published_#13"], 20: reps["published_#20"], 53: reps["published_#53"]}
    exp_fr, exp_tv = scan.load_dot_texture_experiment()
    workers = os.cpu_count()

    # ---------------------------------------------------------------- spike-count-matched thinning
    (out / "thinning").mkdir(exist_ok=True)
    jobs = [(rows[stand[13]], theta, dt, d) for d in DELAYS] + [(rows[stand[k]], theta, dt, 0.0) for k in (20, 53)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        sims = list(pool.map(scan_set, jobs))
    delay_sims, pareto_sims = dict(zip(DELAYS, sims[:len(DELAYS)])), dict(zip((20, 53), sims[len(DELAYS):]))
    base = delay_sims[0]
    delay_rows = {d: dict(evaluate(*delay_sims[d], exp_fr, exp_tv), condition="delay", delay=d, rep=-1)
                  for d in DELAYS}
    pareto_rows = {13: dict(evaluate(*base, exp_fr, exp_tv), condition="pareto", target=13, set=stand[13], rep=-1)}
    for k in (20, 53):
        pareto_rows[k] = dict(evaluate(*pareto_sims[k], exp_fr, exp_tv), condition="pareto", target=k, set=stand[k],
                              rep=-1)
    tjobs = [(("delay", d), r, base, delay_rows[d], 1000 * d + r) for d in DELAYS if d > 0 for r in range(N_REP)]
    tjobs += [(("pareto", k), r, base, pareto_rows[k], 100000 * k + r) for k in (20, 53) for r in range(N_REP)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        thinned = list(pool.map(thin_job, tjobs, chunksize=20))
    d_out = list(delay_rows.values()) + [dict(rec, condition="thinning", delay=key[1], rep=rep)
                                         for key, rep, rec in thinned if key[0] == "delay"]
    p_out = list(pareto_rows.values()) + [dict(rec, condition="thinning", target=key[1], set=stand[13], rep=rep)
                                          for key, rep, rec in thinned if key[0] == "pareto"]
    pd.DataFrame(d_out).to_csv(out / "thinning" / "thinning_all.csv", index=False)
    pd.DataFrame(p_out).to_csv(out / "thinning" / "thinning_pareto_all.csv", index=False)
    print(f"thinning: base set {stand[13]} (as #13), targets {stand[20]} (as #20), {stand[53]} (as #53)", flush=True)

    # ---------------------------------------------------------------- probe: spike counts and per-speed errors
    with ProcessPoolExecutor(max_workers=workers) as pool:
        probe = dict(pool.map(probe_job, [(n, rows[n], theta, dt) for n in params.index]))
    _, actual = lif.load_probe_inputs()
    exp = [c[~np.isnan(c) & (c > 0)] for c in actual.T]
    fsl_obs = np.array([s[0] for s in exp])
    rate_obs = np.array([1000 / np.mean(np.diff(s)) for s in exp])

    # spike counts in 0-100 ms
    (out / "spike_counts").mkdir(exist_ok=True)
    count_rows = []
    sources = [("experiment", None, exp)] + [(f"stand-in for #{k}", reps[f"published_#{k}"],
                                             probe[reps[f"published_#{k}"]]) for k in (1, 13, 20, 53)]
    for label, n, trains in sources:
        for v, s in zip(SPEEDS, trains):
            c = int(np.sum((s > 0) & (s <= 100.0)))
            count_rows.append({"source": label, "set": n, "speed_mm_s": v, "count_0_100_ms": c,
                               "window_rate_spikes_s": c / 0.1, "total_spikes": len(s)})
    counts = pd.DataFrame(count_rows)
    counts.to_csv(out / "spike_counts" / "spike_counts_0_100ms.csv", index=False)
    print(counts.pivot(index="source", columns="speed_mm_s", values="count_0_100_ms").to_string(), flush=True)

    # per-speed errors and objective shares
    (out / "speed_contribution").mkdir(exist_ok=True)
    rrows = []
    for n, trains in probe.items():
        fsl = np.array([s[0] if len(s) else 300.0 for s in trains])
        rate = np.array([1000 / np.mean(np.diff(s)) if len(s) > 1 else np.nan for s in trains])
        se_fsl, se_rate = (fsl_obs - fsl) ** 2, (rate_obs - rate) ** 2
        for k, v in enumerate(SPEEDS):
            rrows.append({"set": n, "speed_mm_s": v, "fsl_sim_ms": fsl[k], "fsl_obs_ms": fsl_obs[k],
                          "fsl_err_ms": fsl[k] - fsl_obs[k], "fsl_share": se_fsl[k] / se_fsl.sum(),
                          "rate_sim_hz": rate[k], "rate_obs_hz": rate_obs[k], "rate_err_hz": rate[k] - rate_obs[k],
                          "rate_share": se_rate[k] / np.nansum(se_rate),
                          "J_fsl": params.J_fsl[n], "J_rate": params.J_rate[n]})
    df = pd.DataFrame(rrows)
    df["delta_fsl_9_5_ms"] = df["set"].map(df[df.speed_mm_s == 9.5].set_index("set")["fsl_err_ms"])
    df.to_csv(out / "speed_contribution" / "per_set_per_speed.csv", index=False)
    j = df.groupby("set").first()
    jf = -(1 - df.assign(e=(df.fsl_obs_ms - df.fsl_sim_ms) ** 2).groupby("set").e.sum() / np.sum(fsl_obs ** 2))
    print(f"per-speed errors: objective check max |J_fsl - params| = {np.abs(jf - j.J_fsl).max():.2e}", flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "unitA")
