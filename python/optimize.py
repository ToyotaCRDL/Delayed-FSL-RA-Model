"""Multi-objective optimisation of the RA model on the probe indentation (counterpart of multiObj_Pareto.m).

Objectives as in plt_main_FST2.m / plt_main_ISI2.m: J_FSL = -(1 - sum (FSL_obs - FSL_sim)^2 / sum FSL_obs^2) and the
same for the firing rate (1 / mean ISI), over the six probe speeds. beta, K_d and V_thr are optimised within the
bounds of the manuscript; T_ref is fixed. NSGA-II (pymoo; population 150, 100 generations, seed 1) is used instead of
MATLAB's gamultiobj. The output is every non-dominated member of the final population, numbered by ascending J_FSL
(from the upper left of the front, as in the manuscript).

Usage: python optimize.py [--theta 0.13] [--tref 2.1] [--dt 0.0125] [--data WORKBOOK] [--relative] [--out NAME]
  default output: results/unitA/params.csv (the 150 sets of the revised version)
  Unit B (Fig. S1):  --data ../LIF_model/stimulus1_probe/finger62-1/f62-1_new.xlsx --out unitB
  relative-error weights (each speed's squared error divided by obs^2):  --relative --out unitA_relative_weights
Output columns: set, beta_mA, kd, t_ref_ms, v_thr_mV (the order of param.mat), J_fsl, J_rate
"""

import argparse

import numpy as np
import pandas as pd
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting

import lif

LB = np.array([-2.5e-6, 0.1e-7, 0.0, 30.0])  # beta, K_d, T_ref, V_thr (multiObj_Pareto.m)
UB = np.array([2.5e-6, 2.5e-7, 1.0, 35.0])
FREE = [0, 1, 3]  # beta, K_d, V_thr (index 2 = T_ref is fixed)
POP, GEN, SEED = 150, 100, 1
FAILED = 1e6  # objective for sets with fewer than two spikes at some speed (MATLAB gives NaN)


class ProbeObjectives(Problem):
    """J_FSL and J_Rate for rows of (beta, K_d, T_ref, V_thr)."""

    def __init__(self, theta=lif.THETA, dt=lif.DT, data_file=None, relative=False):
        """data_file: probe workbook (default Unit A, f61-1.xlsx). The simulated period follows the FE time axis
        of the workbook (200 ms for Unit A, 600 ms for Unit B).
        relative: weight each speed's squared error by 1/observed^2, i.e. J = -(1 - mean squared relative error).
        Default False: equal weights, as in the manuscript."""
        super().__init__(n_var=4, n_obj=2, xl=LB, xu=UB)
        data, actual = lif.load_probe_inputs(data_file)
        finish = float(np.nanmax(data[:, 6])) * 1000.0
        exp = [c[~np.isnan(c) & (c > 0)] for c in actual.T]
        self.fsl_obs = np.array([s[0] for s in exp])
        self.rate_obs = np.array([1 / np.mean(np.diff(s)) for s in exp])
        self.w_fsl = 1 / self.fsl_obs**2 if relative else np.ones(6)
        self.w_rate = 1 / self.rate_obs**2 if relative else np.ones(6)
        self.inputs = [lif.resample_sed(data, k, dt=dt, finish=finish) for k in range(6)]
        self.u_thr = theta * dt

    def responses(self, X):
        """FSL [ms] and 1/mean ISI [1/ms] at each speed for every row of X."""
        fsl = np.zeros((len(X), 6))
        rate = np.zeros((len(X), 6))
        for k, (t, sed) in enumerate(self.inputs):
            cur = lif.afferent_current(sed[:, None], t[1] - t[0], X[:, 0], X[:, 1], self.u_thr)
            for p, s in enumerate(lif.integrate(t, cur, X[:, 2], X[:, 3])):
                fsl[p, k] = s[0] if len(s) else 300.0
                rate[p, k] = 1 / np.mean(np.diff(s)) if len(s) > 1 else np.nan
        return fsl, rate

    def _evaluate(self, X, out, *args, **kwargs):
        fsl, rate = self.responses(X)
        wf, wr = self.w_fsl, self.w_rate
        j_fsl = -(1 - np.sum(wf * (self.fsl_obs - fsl) ** 2, axis=1) / np.sum(wf * self.fsl_obs**2))
        j_rate = -(1 - np.sum(wr * (self.rate_obs - rate) ** 2, axis=1) / np.sum(wr * self.rate_obs**2))
        out["F"] = np.column_stack([j_fsl, np.where(np.isnan(j_rate), FAILED, j_rate)])


class FixedTref(Problem):
    """ProbeObjectives with T_ref held fixed; the decision variables are beta, K_d and V_thr."""

    def __init__(self, theta, dt, t_ref, data_file=None, relative=False):
        super().__init__(n_var=3, n_obj=2, xl=LB[FREE], xu=UB[FREE])
        self.inner, self.t_ref = ProbeObjectives(theta, dt, data_file, relative), t_ref

    def full(self, X):
        return np.insert(X, 2, self.t_ref, axis=1)

    def _evaluate(self, X, out, *args, **kwargs):
        self.inner._evaluate(self.full(X), out)


def main(theta, t_ref, dt, data_file=None, relative=False, out_name="unitA"):
    problem = FixedTref(theta, dt, t_ref, data_file, relative)
    res = minimize(problem, NSGA2(pop_size=POP), ("n_gen", GEN), seed=SEED, verbose=False)
    X, F = problem.full(res.pop.get("X")), res.pop.get("F")
    front = NonDominatedSorting().do(F, only_non_dominated_front=True)
    X, F = X[front], F[front]
    order = np.lexsort((F[:, 1], F[:, 0]))
    X, F = X[order], F[order]
    df = pd.DataFrame(X, columns=lif.PARAM_COLUMNS)
    df.insert(0, "set", np.arange(1, len(df) + 1))
    df["J_fsl"], df["J_rate"] = F[:, 0], F[:, 1]
    out = lif.RESULTS / out_name
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "params.csv", index=False)
    print(f"theta={theta:g} dt={dt:g} T_ref={t_ref:g} ms{' (relative weights)' if relative else ''}: {len(df)} "
          f"non-dominated sets, J_fsl {F[:, 0].min():.3f}..{F[:, 0].max():.3f}, "
          f"J_rate {F[:, 1].min():.3f}..{F[:, 1].max():.3f} -> {out / 'params.csv'}", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--theta", type=float, default=lif.THETA)
    ap.add_argument("--tref", type=float, default=lif.T_REF)
    ap.add_argument("--dt", type=float, default=lif.DT)
    ap.add_argument("--data", default=None)
    ap.add_argument("--relative", action="store_true")
    ap.add_argument("--out", default="unitA")
    a = ap.parse_args()
    main(a.theta, a.tref, a.dt, a.data, a.relative, a.out)
