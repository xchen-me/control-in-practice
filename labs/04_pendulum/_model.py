# Author: Xu Chen, University of Washington
"""Phone pendulum: phyphox gyroscope loader, pendulum physics, and a phyphox-format log simulator."""
import io
import os
import numpy as np
import pandas as pd
from scipy import integrate

g = 9.81                                        # m/s^2


def L_equivalent(d, h, w):
    """Equivalent simple-pendulum length of a phone (h x w plate, centre at distance d from the pivot)
    swinging in the plane of its screen: L_eq = I_pivot/(m d) = d + (h^2 + w^2)/(12 d)."""
    return d + (h ** 2 + w ** 2) / (12 * d)


def period_large_angle(T0, amp):
    """Period of a pendulum at amplitude amp (rad): T0 (1 + amp^2/16 + 11 amp^4/3072 + ...)."""
    amp = np.asarray(amp, float)
    return T0 * (1 + amp ** 2 / 16 + 11 * amp ** 4 / 3072)


# ---------------------------------------------------------------- loading phyphox exports
def load_phyphox(path, axis=None):
    """Read a phyphox gyroscope CSV export in any of its CSV formats: comma, tab or semicolon separated,
    with a decimal point or (tab and semicolon only) a decimal comma. Both are detected from the file.

    Columns like "Time (s)", "Gyroscope x (rad/s)", ... and optionally "Absolute (rad/s)" (ignored for the
    axis choice, it is not signed). The axis with the most energy (variance) is picked unless axis='x'/'y'/'z'.
    Samples are resampled onto a uniform grid at the median sample time (phone time stamps jitter).
    Returns dict(t, w, axis, fs, raw) with t starting at 0 and w in rad/s."""
    with open(path, encoding='utf-8-sig') as fh:
        text = fh.read()
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith('#')]
    head = lines[0]
    sep = '\t' if '\t' in head else (';' if ';' in head else ',')
    # with tab or semicolon separators a comma inside the numbers is a decimal comma
    decimal = ',' if sep != ',' and any(',' in ln for ln in lines[1:20]) else '.'
    raw = pd.read_csv(io.StringIO('\n'.join(lines)), sep=sep, decimal=decimal)
    raw.columns = [c.strip().strip('"') for c in raw.columns]
    tcol = next(c for c in raw.columns if 'time' in c.lower())
    gyro = {}
    for c in raw.columns:
        lc = c.lower()
        if c == tcol or 'abs' in lc:
            continue
        for ax_ in 'xyz':
            if (' %s ' % ax_) in ' %s ' % lc.replace('(', ' (') or lc.endswith(' ' + ax_):
                gyro[ax_] = c
    if not gyro:
        raise ValueError('no gyroscope x/y/z columns found in %s: %s' % (path, list(raw.columns)))
    raw = raw.dropna(subset=[tcol] + list(gyro.values()))
    t_raw = raw[tcol].to_numpy(float)
    order = np.argsort(t_raw)
    t_raw = t_raw[order] - t_raw[order][0]
    energy = {a: np.var(raw[c].to_numpy(float)) for a, c in gyro.items()}
    axis = axis or max(energy, key=energy.get)
    dt = float(np.median(np.diff(t_raw)))
    t = np.arange(0, t_raw[-1], dt)
    w = np.interp(t, t_raw, raw[gyro[axis]].to_numpy(float)[order])
    return dict(t=t, w=w, axis=axis, fs=1 / dt, energy=energy, raw=raw, columns=gyro)


def find_release(t, w, quiet=0.05):
    """Index of the release: the last quiet sample (|w| < quiet rad/s) before the first large swing."""
    big = np.nonzero(np.abs(w) > 0.3 * np.max(np.abs(w)))[0][0]
    q = np.nonzero(np.abs(w[:big]) < quiet)[0]
    return int(q[-1]) if len(q) else 0


def angle_from_gyro(t, w, k_rel, hold_min=0.5):
    """Angle (rad) by integrating the gyro after removing its bias; angle is 0 at the start of the record.

    The bias is first estimated while the phone is held still before the release (at least hold_min s).
    A small residual bias makes the integrated angle drift linearly; the equilibrium of a pendulum does
    not move, so the slope of a straight line fitted to the angle after the release is removed as well.
    Returns (angle, total bias in rad/s)."""
    if t[k_rel] - 0.2 >= hold_min:
        bias = float(np.mean(w[t < t[k_rel] - 0.2]))
    else:
        bias = 0.0
        print('warning: no still segment before the release; bias from the drift fit only')
    th = integrate.cumulative_trapezoid(w - bias, t, initial=0.0)
    slope = np.polyfit(t[k_rel:], th[k_rel:], 1)[0]
    bias += slope
    th = integrate.cumulative_trapezoid(w - bias, t, initial=0.0)
    return th, bias


def cycles(t, th, k_rel, min_amp=0.005):
    """Cycle-by-cycle period (s) and amplitude (rad) of the angle th after the release.

    Periods: upward zero crossings of th about its equilibrium (the mean of the last 20 % of the record),
    linearly interpolated; the angle is fast there and smooth (integrated), so noise barely moves them.
    Amplitude: half the peak-to-peak angle within the cycle. The list stops at the first cycle below
    min_amp (rad) or with a period more than 10 % off the median.
    Returns (cycle mid-times, periods, amplitudes)."""
    x = th[k_rel:] - np.mean(th[-max(len(th) // 5, 10):])
    tt = t[k_rel:]
    i = np.nonzero((x[:-1] < 0) & (x[1:] >= 0))[0]
    tc = tt[i] - x[i] * (tt[i + 1] - tt[i]) / (x[i + 1] - x[i])
    T = np.diff(tc)
    amp = np.array([(np.max(x[p:q]) - np.min(x[p:q])) / 2 for p, q in zip(i[:-1], i[1:])])
    bad = np.nonzero((amp < min_amp) | (np.abs(T / np.median(T) - 1) > 0.1))[0]
    n = bad[0] if len(bad) else len(T)
    return (tc[:-1] + T / 2)[:n], T[:n], amp[:n]


# ---------------------------------------------------------------- simulated phyphox log
def simulate_phyphox(d, theta0, T_rec, seed, h=0.15, w_ph=0.072, zeta_lin=0.006, c_quad=0.08,
                     axis='z', t_hold=3.0, fs=100.0, bias=(-0.005, 0.0015, -0.0025), noise=0.0015,
                     quant=0.00106, misalign_deg=4.0, twist=0.06):
    """Free decay of a phone pendulum, as the phyphox 'Gyroscope (rotation rate)' experiment exports it.

    Nonlinear model: theta'' = -(g/L_eq) sin(theta) - 2 zeta_lin wn theta' - c_quad |theta'| theta'
    (linear pivot/string friction plus quadratic air drag, c_quad in 1/rad). The phone is held still at
    theta0 (hand tremor) for t_hold s, then released. Gyro: bias, white noise, 16-bit quantization, a few
    degrees of misalignment between the swing axis and the phone axis, a small torsional (twist) mode,
    time stamps with jitter around 1/fs."""
    rng = np.random.default_rng(seed)
    Leq = L_equivalent(d, h, w_ph)
    wn = np.sqrt(g / Leq)

    def f(_, x):
        th, om = x
        return [om, -wn ** 2 * np.sin(th) - 2 * zeta_lin * wn * om - c_quad * abs(om) * om]

    n = int(T_rec * fs)
    t = np.sort(np.arange(n) / fs + rng.normal(0, 0.15e-3, n))
    t -= t[0]
    t_after = t[t >= t_hold] - t_hold
    sol = integrate.solve_ivp(f, (0, t_after[-1]), [theta0, 0.0], t_eval=t_after, rtol=1e-9, atol=1e-11)
    om = np.zeros(n)
    om[t >= t_hold] = sol.y[1]
    hold = t < t_hold                            # hand tremor: small zero-mean wobble
    om[hold] = 0.012 * np.sin(2 * np.pi * 7.3 * t[hold]) + 0.008 * rng.standard_normal(hold.sum())
    a = np.radians(misalign_deg)
    # torsional mode of the string (about the vertical axis), excited by the release
    wt = 2 * np.pi * 0.55
    tw = np.where(t >= t_hold, twist * np.exp(-0.03 * (t - t_hold)) * np.sin(wt * (t - t_hold)), 0.0)
    main, cross1, cross2 = om * np.cos(a), om * np.sin(a) + tw, 0.3 * tw
    others = [c for c in 'xyz' if c != axis]
    sig = {axis: main, others[0]: cross1, others[1]: cross2}
    out = {}
    for k, c in enumerate('xyz'):
        v = sig[c] + bias[k] + noise * rng.standard_normal(n)
        out[c] = quant * np.round(v / quant)
    df = pd.DataFrame({'Time (s)': t, 'Gyroscope x (rad/s)': out['x'], 'Gyroscope y (rad/s)': out['y'],
                       'Gyroscope z (rad/s)': out['z']})
    df['Absolute (rad/s)'] = np.sqrt(out['x'] ** 2 + out['y'] ** 2 + out['z'] ** 2)
    return df, dict(L_eq=Leq, wn_small=wn, T0=2 * np.pi / wn)


def data_path(name):
    """Path of data/<name>.csv in the repository."""
    from ctrlprac import data
    return os.path.join(data.DATA, name + '.csv')
