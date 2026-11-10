# Python implementation

The results of the revised manuscript were computed with these scripts. The model is the one implemented in the
MATLAB scripts of `../LIF_model` (`spike*.m`, `runge_kutta4.m`, `func.m` and the input preparation of the `main*.m`
scripts), ported line by line; the `verify_v1_*.py` scripts reproduce the results of the original submission (tag
`v1.0`) from the v1.0 settings.

Settings of the revised version (`lif.py`), the same for every stimulus:

| | revised version | original submission (v1.0) |
|---|---|---|
| SED threshold | rate of change, \|dU/dt\| > θ = 0.13 J m⁻³ ms⁻¹ | increment per time step > 1.1 J m⁻³ |
| time step | 0.0125 ms; the FE SED (every 2.5 ms) is spline-resampled in time | 0.5 ms (probe), 2.5 ms (scanning) |
| refractory period | fixed at 2.1 ms | optimised (0–1 ms) |
| optimisation | NSGA-II (pymoo), population 150, 100 generations, seed 1; β, K_d, V_thr | gamultiobj, population 150, 100 generations; β, K_d, T_ref, V_thr |
| parameter sets | all 150 non-dominated sets of the final population | 53 sets returned by gamultiobj |

The SA1 model (Stimulus 2–4, comparison with RA) is unchanged.


## Setup

```bash
python3.14 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Input data are read from `../LIF_model`: `stimulus1_probe/finger61-1/f61-1.xlsx` (probe, Unit A),
`stimulus1_probe/finger62-1/f62-1_new.xlsx` (probe, Unit B), and `dataSED.mat`, `dataTime.mat`, `data_coord.xlsx`
and `expFireRate.mat` / `expTempVar.mat` of `stimulus2_braille/`, `stimulus3_dot_texture/` and
`stimulus4_single_dot/`. The `.rpt` and `.xlsx` FE outputs in the `finger*/` folders are not needed.
The scripts use all CPU cores (`RA_N_WORKERS` limits `evaluate_front.py` and `original_sets_dt.py`).


## Workflow

Run in this folder, in this order. Results are written to `results/<name>/` (`unitA` unless stated otherwise).

| step | script | output | time (18 cores) |
|---|---|---|---|
| 1 | `optimize.py` | `results/unitA/params.csv`: 150 sets (β, K_d, T_ref, V_thr, J_FSL, J_Rate) | 1 min |
| 2 | `run_all.py` | `summary.csv`, `representative.json`, `delay_control.csv`; all SEPs and metrics in `../LIF_model/*/result1`, `result2`, `param.mat`, `popCost.mat` and `../LIF_model/statistical analysis/*.csv` (v1.0 formats) | 5 min |
| 3 | `probe_summary.py` | `probe_summary.csv`: FSL and firing rate at each speed, ΔFSL, RMS rate error | < 1 min |
| 4 | `theta_basis.py` | `theta_basis.csv`: residual SED fluctuation after the probe indentation (basis of θ) and a check that no set fires during the hold | < 1 min |
| 4b | `hold_spikes.py` | hold-phase spikes (more than 10 ms after the end of the 1 mm indentation): `hold_spikes.csv` for the 150 revised sets (θ = 0.13) and `results/v1.0_sets/hold_spikes_theta0.csv` for the 53 v1.0 sets without the threshold (θ = 0); needs the v1.0 checkout | < 1 min |
| 5 | `dt_series.py` | `dt_series/`: probe responses of the 150 sets at time steps from 0.5 to 0.0125 ms | |
| 6 | `response_analyses.py` | `thinning/`: spike-count-matched random thinning; `spike_counts/`: spike counts in 0–100 ms; `speed_contribution/`: per-speed errors and objective shares | |
| 7 | `evaluate_front.py` | `front_summary.csv`: per-set metrics with the numbers of responding afferents | |
| 8 | `optimize.py --relative --out unitA_relative_weights`, `evaluate_front.py unitA_relative_weights`, `relative_weights.py` | `results/unitA_relative_weights/`: front optimised with relative-error weights across speeds, and `comparison.csv` | |
| 9 | `optimize.py --data ../LIF_model/stimulus1_probe/finger62-1/f62-1_new.xlsx --out unitB`, `probe_summary.py unitB --data ...` | `results/unitB/`: Unit B (Fig. S1) | 4 min |

Representative sets (`representative.json`; the keys name the sets of the original submission they replace):
set 1 (best FSL), set 87 (best balance of FSL and firing rate: closest to the ideal point after scaling both
objectives to [0, 1]), and the sets with ΔFSL closest to 10 and 30 ms (92 and 118). The control analysis uses set 87.

Noise analysis: the v1.0 scripts add N(0, γ²) current noise per 2.5 ms step; the standard deviation is scaled by
√(2.5 / dt) to keep its effect on the membrane potential. The random streams differ from MATLAB's.

RW is NaN where it cannot be determined (too few spikes for the 25% crossings), e.g. the 4.1 mm dot for sets with
large ΔFSL.


## Scripts for the original submission (v1.0)

These read files of the original submission from a checkout of tag `v1.0`, located at `../v1.0` or at the path in
the environment variable `RA_V1_DIR`, e.g.

```bash
git worktree add v1.0 v1.0
```

in the repository root.

| script | purpose |
|---|---|
| `verify_v1_probe.py` | reproduces the v1.0 probe spike trains of the 53 sets (318 conditions) |
| `verify_v1_scan.py` | reproduces the v1.0 SEPs of the dot pattern and single dots, FR/TV RMSE, RW and spike width |
| `verify_v1_extra.py` | reproduces the v1.0 Braille SEPs and the control analysis; noise compared in distribution |
| `original_sets_dt.py` | probe responses of the 53 v1.0 sets at finer time steps (`results/v1.0_sets/`) |
| `hold_spikes.py` | hold-phase spikes of the 53 v1.0 sets without the threshold, at dt = 0.0125 ms (`results/v1.0_sets/hold_spikes_theta0.csv`) |


## Files

| file | contents |
|---|---|
| `lif.py` | LIF model, probe input, settings |
| `scan.py` | scanning stimuli (Braille, dot pattern, single dot), imposed unresponsive period, noise, FR/TV/RW and spike width |
| `optimize.py` | objectives and NSGA-II |
| `run_all.py` | all analyses of the manuscript for a parameter-set series |
| other scripts | see the tables above |
