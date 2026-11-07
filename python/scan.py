"""Scanning stimuli (Python implementation of mainBrailleRA.m, mainDotTexRA.m and mainSingleDotRA.m, incl. the
_delay, _noise and _modi_RW variants).

The SED on the FE surface grid is linearly interpolated (griddata) onto the 51 receptors on the finger centre line
and fed to the LIF model of lif.py (spikeRA.m is identical to spike.m). In the revised version (dt given) the
centre-line SED, output by the FE model every 2.5 ms, is spline-resampled in time onto the common time step; with
dt=None the FE step is used as in v1.0. Metric code (FR, TV, RW) follows the MATLAB loops line by line, because
bin-edge comparisons and the ill-conditioned 8th-degree polyfit are sensitive to the exact arithmetic.
"""

from functools import lru_cache

import numpy as np
import openpyxl
import scipy.io as sio
import scipy.linalg
import scipy.signal
from scipy.interpolate import CubicSpline, LinearNDInterpolator

import lif

LIF_DIR = lif.LIF_DIR
BRAILLE_DIR = LIF_DIR / "stimulus2_braille"
DOT_TEXTURE_DIR = LIF_DIR / "stimulus3_dot_texture"
SINGLE_DOT_DIR = LIF_DIR / "stimulus4_single_dot"

# (label used in the released file names, FE model folder), in the loop order of the MATLAB scripts
BRAILLE_FINGER = "57-7-3"
DOT_SPACINGS = [("6-2", "58-28"), ("5-2", "58-27"), ("4-3", "58-26"), ("3-2", "58-25"), ("2-4", "58-24"), ("1-3", "58-23")]
DOT_WIDTHS = [("4-1", "70-2"), ("2-6", "70-3"), ("1-6", "70-4")]

# How each main*RA.m reads data_coord.xlsx: (sheet for X, sign of X, sheet for Y, which Y entries replace the
# first/last query rows). The dot FE models swap X/Y and flip X; Braille uses the sheets as they are.
COORD_LAYOUT = {
    BRAILLE_DIR: (0, 1.0, 1, ((0, 0), (-1, 0))),
    DOT_TEXTURE_DIR: (1, -1.0, 0, ((0, 0), (-1, 0))),
    SINGLE_DOT_DIR: (1, -1.0, 0, ((0, -1), (0, 0))),
}

CENTRE_COLUMN = 25  # 0-based column of the 51x51 resampled grid; MATLAB nodes 1276-1326 are its 51 rows
N_RECEPTORS = 51


def matlab_colon(a, d, b):
    """a:d:b with MATLAB's colon arithmetic (symmetric about the midpoint), so bin edges match bit for bit."""
    tol = 2.0 * np.finfo(float).eps * max(abs(a), abs(b))
    sig = np.sign(d)
    n = int(np.floor((b - a) / d))
    if sig * (b - (a + (n + 1) * d)) >= -tol:
        n += 1
    last = a + n * d
    if sig * (last - b) > -tol:
        last = b
    v = np.zeros(n + 1)
    k = np.arange(n // 2 + 1)
    v[k] = a + k * d
    v[n - k] = last - k * d
    if n % 2 == 0:
        v[n // 2] = (a + last) / 2
    return v


def _sheet(path, index):
    ws = openpyxl.load_workbook(path, data_only=True).worksheets[index]
    return np.array([[float(v) for v in row] for row in ws.iter_rows(values_only=True)])


@lru_cache(maxsize=None)
def centre_line_sed(stim_dir, finger):
    """SED at the 51 centre-line receptors: returns t [ms], sed (time, 51), receptor y [mm]."""
    x_sheet, x_sign, y_sheet, (first, last) = COORD_LAYOUT[stim_dir]
    coord_x = x_sign * _sheet(stim_dir / "data_coord.xlsx", x_sheet)
    coord_y = _sheet(stim_dir / "data_coord.xlsx", y_sheet)
    fe = stim_dir / f"finger{finger}"
    data_sed = sio.loadmat(fe / "dataSED.mat")["dataSED"]
    data_time = sio.loadmat(fe / "dataTime.mat")["dataTime"]

    resample = matlab_colon(5.0, -0.2, -5.0)
    y_query = resample.copy()
    y_query[0], y_query[-1] = coord_y[first], coord_y[last]
    x_query = np.full(N_RECEPTORS, -1.0 * resample[CENTRE_COLUMN])

    # reshape(dataSED(i,:),[65,65]) pairs node k with coord(:)(k) in column-major order
    points = np.column_stack([coord_x.ravel(order="F"), coord_y.ravel(order="F")])
    interp = LinearNDInterpolator(points, data_sed.T)
    sed = interp(np.column_stack([x_query, y_query])).T
    t = data_time[:, 0] * 1000.0
    return t, sed, y_query


def resample_time(t, sed, dt):
    """Spline-resample the 2.5 ms FE output onto a finer uniform grid, as mainProbe.m does for the probe."""
    n = int(round((t[-1] - t[0]) / dt)) + 1
    t_new = np.linspace(t[0], t[-1], n)
    return t_new, CubicSpline(t, sed, axis=0)(t_new)


def run(stim_dir, finger, params_row, u_thr=lif.U_THR, dt=None, delay=0.0, noise=0.0, rng=None):
    """Simulate the 51 centre-line afferents. dt=None keeps the FE time step (2.5 ms) as in v1.0.

    delay > 0 applies the imposed unresponsive period of the *_delay scripts; noise > 0 adds the Gaussian
    current noise of the *_noise scripts (standard deviation in mA per time step).
    """
    t, sed, _ = centre_line_sed(stim_dir, finger)
    if dt is not None:
        t, sed = resample_time(t, sed, dt)
    beta, kd, t_ref, v_thr = params_row
    cur = lif.afferent_current(sed, t[1] - t[0], beta, kd, u_thr)
    if delay > 0:
        cur = suppress_onsets(t, cur, delay)
    if noise > 0:
        cur = cur + noise * (rng if rng is not None else np.random.default_rng()).standard_normal(cur.shape)
    return lif.integrate(t, cur, t_ref, v_thr)


def islocalmin(y):
    """MATLAB islocalmin with default options: interior minima, a flat minimum marked at its centre."""
    flags = np.zeros(len(y), dtype=bool)
    i = 1
    while i < len(y) - 1:
        if y[i] < y[i - 1]:
            j = i
            while j + 1 < len(y) and y[j + 1] == y[i]:
                j += 1
            if j + 1 < len(y) and y[j + 1] > y[i]:
                flags[(i + j) // 2] = True
            i = j + 1
        else:
            i += 1
    return flags


def suppress_onsets(t, cur, delay):
    """spikeRA.m of the *_delay scripts: zero the current for `delay` ms from each local minimum of the
    25 Hz low-pass filtered current, i.e. from the onset of each loading/unloading response."""
    fs = 1.0 / ((t[2] - t[1]) / 1000.0)
    b, a = scipy.signal.butter(4, 25.0 / (fs / 2.0), "low")
    pad = 3 * (max(len(a), len(b)) - 1)  # MATLAB filtfilt's padding length
    filtered = scipy.signal.filtfilt(b, a, cur, axis=0, padtype="odd", padlen=pad)
    out = cur.copy()
    for col in range(cur.shape[1]):
        for m in t[islocalmin(filtered[:, col])]:
            out[(t >= m) & (t <= m + delay), col] = 0.0
    return out


# ---------------------------------------------------------------- Stimulus 3: FR and TV (mainDotTexRA.m)

def firing_rate_histogram(spikes, edges):
    """allFiringRateHistogram for one afferent; `spikes` [ms] without the 0 placeholder."""
    sp = np.concatenate([[0.0], spikes])  # allSpikeTime keeps spike.m's leading 0
    hist = np.zeros(len(edges) - 1)
    for i in range(len(edges) - 1):
        counts = 0.0
        index = 0
        for j in range(len(sp)):
            tj = sp[j] / 1000
            if edges[i] <= tj and tj <= edges[i + 1] and index == 0:
                counts = (tj - edges[i]) / (tj - sp[j - 1] / 1000)
                index = 1
            elif edges[i] <= tj and tj <= edges[i + 1] and index == 1:
                counts = counts + 1
            elif edges[i + 1] <= tj and index == 1:
                counts = counts + (edges[i + 1] - sp[j - 1] / 1000) / (tj - sp[j - 1] / 1000)
                index = 0
        hist[i] = counts / 0.01
    return hist


def dot_texture_metrics(spike_trains, finger):
    """Mean firing rate (FR) and temporal variation (TV) over the 51 afferents for one dot spacing."""
    edges = matlab_colon(1.15, 0.01, 1.65) if finger == "58-27" else matlab_colon(1.2, 0.01, 1.7)
    hist = np.column_stack([firing_rate_histogram(s, edges) for s in spike_trains])
    nb = hist.shape[0]
    tempo_vari = 0.0
    for i in range(hist.shape[1]):
        pre = 0.0
        for j in range(nb - 13):
            pre = pre + abs(hist[j + 13, i] - hist[j, i]) / (nb - 13)
        tempo_vari = tempo_vari + pre / hist.shape[1]
    return hist.mean(), tempo_vari


def load_dot_texture_experiment():
    """Experimental FR and TV, ordered 1.3, 2.4, 3.2, 4.3, 5.2, 6.2 mm (vec(7-mlp) in MATLAB)."""
    fr = sio.loadmat(DOT_TEXTURE_DIR / "expFireRate.mat")["expFireRate"].ravel()
    tv = sio.loadmat(DOT_TEXTURE_DIR / "expTempVar.mat")["expTempVar"].ravel()
    return fr, tv


def rmse(actual, predicted):
    return np.sqrt(np.mean((actual - predicted) ** 2))


# ---------------------------------------------------------------- Stimulus 4: RW (mainSingleDotRA.m)

STIMULUS_PERIOD = {"70-2": (1245.0, 1795.0), "70-3": (1295.0, 1750.0), "70-4": (1317.5, 1727.5)}
RW_RECEPTORS = range(1298 - 1276, 1304 - 1276 + 1)  # MATLAB nodes 1298-1304


def averaged_rates(spikes):
    """allFiringRate2: mean of the instantaneous rates on either side of each spike (0 for the first)."""
    rate = np.zeros(len(spikes))
    rate[1:] = 1.0 / (np.diff(spikes) * 0.001)
    out = np.zeros(len(spikes))
    for i in range(len(spikes)):
        if i == 0:
            out[i] = 0
        elif i == len(spikes) - 1:
            out[i] = rate[i]
        else:
            out[i] = (rate[i] + rate[i + 1]) * 0.5
    return out


def matlab_polyfit(x, y, n):
    """polyfit as [Q,R] = qr(V,0); p = R\\(Q'*y) (no pivoting).

    The Vandermonde matrix for t = 1000-2000 ms and degree 8 is numerically rank deficient; this form
    reproduces the released RW values, whereas pivoted/rank-truncated solvers give no valid roots.
    """
    V = np.ones((len(x), n + 1))
    for j in range(n - 1, -1, -1):
        V[:, j] = x * V[:, j + 1]
    Q, Rm = np.linalg.qr(V)
    return scipy.linalg.solve_triangular(Rm, Q.T @ y)


def response_width(spike_trains, finger, polyfit=matlab_polyfit):
    """RW [mm]: 25% width of the 8th-degree fit to the binned firing rate of 7 centre afferents."""
    lo, hi = STIMULUS_PERIOD[finger]
    rates = {r: averaged_rates(spike_trains[r]) for r in RW_RECEPTORS}
    bin_rate = np.zeros(51)
    bin_time = np.zeros(51)
    for i in range(51):
        total = 0.0
        count = 0
        for r in RW_RECEPTORS:
            for s, fr in zip(spike_trains[r], rates[r]):
                if i * 20 + 1000 <= s < (i + 1) * 20 + 1000:
                    if s <= lo or hi <= s:
                        continue
                    total = total + fr
                    count = count + 1
        bin_rate[i] = 0 if count == 0 else total / count
        bin_time[i] = i * 20 + 1000
    nz = np.flatnonzero(bin_rate != 0)
    st, en = nz[0] - 1, nz[-1] + 1
    x = bin_time[st:en + 1]
    p = polyfit(x, bin_rate[st:en + 1], 8)
    y_target = np.max(np.polyval(p, x)) * 0.25
    p_adj = p.copy()
    p_adj[-1] -= y_target
    roots = np.roots(p_adj)
    real = np.sort(roots[roots.imag == 0].real)[::-1]
    real = real[(bin_time[st] < real) & (real < bin_time[en])]
    return (real[0] - real[-1]) / 1000.0 * 10.0


def spike_width(spike_trains):
    """Supplementary Fig. S6 measure (stimulus4_single_dot_modi_RW): earliest-to-latest spike over all 51 afferents [mm]."""
    s = np.concatenate(spike_trains)
    return (s.max() - s.min()) / 1000.0 * 10.0
