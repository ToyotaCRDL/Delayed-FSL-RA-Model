"""LIF model of RA afferents (Python implementation of LIF_model/*/spike*.m, runge_kutta4.m and func.m).

Afferent current (Eq. 1 of the manuscript):
    I(t) = beta + K_d |dU/dt|   if |dU/dt| > theta,   I(t) = 0 otherwise,
where U is the strain energy density (SED) from the FE model. As in the MATLAB code, the threshold is applied to
the SED increment per time step, u_thr = theta * dt, and the membrane equation is integrated with RK4, the current
held constant over each step. After a spike the potential is reset to 0 and held for the refractory period T_ref.

Settings of the revised version (all stimuli): THETA, DT and T_REF below. The original submission (tag v1.0) used a
per-step threshold of 1.1 J m^-3 (U_THR) with dt = 0.5 ms for the probe and 2.5 ms for scanning; U_THR is kept as the
default of the functions so that the verify_v1_*.py scripts reproduce the v1.0 results.

Input data are read from ../LIF_model. Files of the original submission (parameter sets and results of v1.0) are
read from a checkout of tag v1.0 (V1_DIR; see README.md), only by verify_v1_*.py and original_sets_dt.py.
"""

import os
from pathlib import Path

import numpy as np
import openpyxl
import scipy.io as sio
from scipy.interpolate import CubicSpline

REPO = Path(__file__).resolve().parents[1]
LIF_DIR = REPO / "LIF_model"
PROBE_DIR = LIF_DIR / "stimulus1_probe"
RESULTS = Path(__file__).resolve().parent / "results"
V1_DIR = Path(os.environ.get("RA_V1_DIR", REPO / "v1.0"))
V1_LIF_DIR = V1_DIR / "LIF_model"

# Constants shared by every MATLAB script
R = 66.11e6
C = 1.01e-6
TAU = 66.77  # ms, R*C
M = 1.0

# Revised version
THETA = 0.13  # SED-rate threshold [J m^-3 ms^-1]
DT = 0.0125  # time step for every stimulus [ms]
T_REF = 2.1  # refractory period [ms] (Brink and Mackel, 1993)

# Original submission (v1.0)
U_THR = 1.1  # SED-increment threshold per time step [J m^-3]

PROBE_SPEEDS = ["9-5", "14-8", "20-6", "34-0", "46-5", "60-6"]  # mm/s, column order of f61-1.xlsx
PARAM_COLUMNS = ["beta_mA", "kd", "t_ref_ms", "v_thr_mV"]  # column order of param.mat


def u_thr(theta=THETA, dt=DT):
    """Per-step SED threshold equivalent to the rate threshold theta."""
    return theta * dt


def load_v1_params():
    """The 53 parameter sets of the original submission (rows of beta, K_d, T_ref, V_thr), from tag v1.0."""
    return sio.loadmat(V1_LIF_DIR / "stimulus1_probe" / "param.mat")["param"]


def _read_sheet(path, sheet):
    ws = openpyxl.load_workbook(path, data_only=True)[sheet]
    rows = [[np.nan if v is None else float(v) for v in row] for row in ws.iter_rows(values_only=True)]
    return np.array(rows)


def load_probe_inputs(path=None):
    """SED histories (6 speeds + time [s]) and experimental spike times [ms] from f61-1.xlsx (Unit A).

    path: another workbook in the same layout (Unit B: LIF_model/stimulus1_probe/finger62-1/f62-1_new.xlsx).
    """
    path = Path(path) if path is not None else PROBE_DIR / "finger61-1" / "f61-1.xlsx"
    return _read_sheet(path, "Sheet1"), _read_sheet(path, "Sheet3")


def resample_sed(data, speed_idx, dt=0.5, start=0.0, finish=200.0):
    """Spline-resample the FE SED (sampled every 2.5 ms) onto a uniform grid, as in mainProbe.m.

    MATLAB griddedInterpolant(...,'spline') is a not-a-knot cubic spline, like CubicSpline's default.
    """
    x = data[:, 6] * 1000.0
    v = data[:, speed_idx]
    point = int(round((finish - start) / dt)) + 1
    t = np.linspace(start, finish, point)
    return t, CubicSpline(x, v)(t)


def afferent_current(sed, dt, beta, kd, u_thr=U_THR):
    """spike.m: I = beta + M*K_d*|dU|/dt only when |dU| per step exceeds u_thr, otherwise 0.

    `sed` is (time,) or (time, afferent); differences are taken along time. `beta` and `kd` may be arrays
    over afferents (e.g. one parameter set per column when a whole GA population is evaluated at once).
    """
    du = np.diff(sed, axis=0)
    active = np.abs(du) > u_thr
    rest = np.where(active, beta + M * kd * np.abs(du) / dt, 0.0)
    return np.concatenate([np.zeros_like(rest[:1]), rest], axis=0)


def integrate(t, cur, t_ref, v_thr):
    """spike.m main loop (RK4 with the current held constant over each step) for every column of `cur`.

    `t` is the common time grid [ms]; `cur` is (time, afferent); `t_ref` and `v_thr` may be per-afferent arrays.
    Returns one array of spike times [ms] per afferent, without the leading 0 placeholder of t_spike(1).
    """
    dt = t[1] - t[0]
    n = cur.shape[1]

    def f(vn, i):
        return (1.0 / TAU) * (R * cur[i] - vn)

    spikes = [[] for _ in range(n)]
    spike_time_last = -np.broadcast_to(np.asarray(t_ref, dtype=float), (n,)).copy()
    tn = 0.0
    vn = np.zeros(n)
    i = 0
    while tn < t[-1] - 1:
        k1 = f(vn, i)
        k2 = f(vn + dt / 2.0 * k1, i)
        k3 = f(vn + dt / 2.0 * k2, i)
        k4 = f(vn + dt * k3, i)
        dv = dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        integrating = tn >= spike_time_last + t_ref
        vn = np.where(integrating, vn + dv, vn)
        tn = t[i + 1]
        fired = vn > v_thr
        if fired.any():
            vn[fired] = 0.0
            spike_time_last[fired] = tn
            for a in np.flatnonzero(fired):
                spikes[a].append(tn)
        i += 1
    return [np.array(s) for s in spikes]


def simulate_population(t, sed, beta, kd, t_ref, v_thr, u_thr=U_THR):
    """Current (afferent_current) followed by the LIF loop (integrate); `sed` is (time, afferent)."""
    dt = t[1] - t[0]
    return integrate(t, afferent_current(sed, dt, beta, kd, u_thr), t_ref, v_thr)


def simulate(t, sed, beta, kd, t_ref, v_thr, u_thr=U_THR):
    """Single-afferent version of simulate_population. Returns spike times [ms]."""
    return simulate_population(t, sed[:, None], beta, kd, t_ref, v_thr, u_thr)[0]


def run_probe(params_row, speed_idx, dt=0.5, u_thr=U_THR):
    data, _ = load_probe_inputs()
    t, sed = resample_sed(data, speed_idx, dt=dt)
    beta, kd, t_ref, v_thr = params_row
    return simulate(t, sed, beta, kd, t_ref, v_thr, u_thr)
