# Author: Xu Chen, University of Washington
"""Lab 01 helpers: the 'true' heater boards, root-locus gain search, and loop margins with an exact delay."""
import numpy as np
import control as ct

TS = 1.0          # s, sample time of the board firmware
QUANT = 0.1       # degC, sensor resolution
NOISE = 0.04      # degC, sensor noise (std)
U_LIM = (0.0, 100.0)

# Representative values in the style of the TCLab kit (not a calibrated datasheet).
# Heater block (tau1) and thermistor lump (tau2) in series, plus a short transport delay.
BOARD_A = dict(name='board A (demo)', K=0.62, tau1=150.0, tau2=15.0, delay=4.0, amb=21.0)


def true_plant(b):
    """Two-lump plant K / ((tau1 s + 1)(tau2 s + 1)) in degC per % heater; the delay is separate."""
    return ct.tf([b['K']], np.polymul([b['tau1'], 1], [b['tau2'], 1]))


def fopdt_pade(fit, order=1):
    """Identified FOPDT with the delay replaced by a Pade approximation (rational, for root locus)."""
    num, den = ct.pade(fit['theta'], order)
    return fit['model'] * ct.tf(num, den)


def dominant(poles):
    """Closed-loop pole with the largest real part (the slowest one)."""
    poles = np.asarray(poles)
    return poles[np.argmax(poles.real)]


def gain_for_zeta(L, zeta, kmax=1e3):
    """Smallest root-locus gain k for which the slowest complex closed-loop pair of 1 + k L(s) reaches damping zeta.

    Returns (k, closed-loop poles)."""
    num, den = np.ravel(L.num[0][0]), np.ravel(L.den[0][0])

    def z_of(k):                             # damping of the slowest complex pair (1 if none)
        p = np.roots(np.polyadd(den, k * num))
        p = p[np.abs(p.imag) > 1e-9]
        if not len(p):
            return 1.0
        p = dominant(p)
        return -p.real / abs(p)

    ks = np.logspace(-4, np.log10(kmax), 4000)
    zs = np.array([z_of(k) for k in ks])
    i = np.argmax(zs <= zeta)
    if zs[i] > zeta:
        raise ValueError('damping %.2f not reached' % zeta)
    lo, hi = ks[max(i - 1, 0)], ks[i]
    for _ in range(60):                      # bisection on the damping of the dominant pair
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if z_of(mid) > zeta else (lo, mid)
    k = hi
    return k, np.roots(np.polyadd(den, k * num))


def freq(L, w, delay=0.0):
    """L(jw) e^{-jw delay} as a complex array."""
    return np.asarray(L(1j * w)).ravel() * np.exp(-1j * w * delay)


def margins(L, delay=0.0, w=None):
    """Gain crossover wc, PM (deg), phase crossover wpc, GM (abs), delay margin (s) of L(s) e^{-s delay}."""
    w = np.logspace(-5, 1, 20000) if w is None else w
    H = freq(L, w, delay)
    mag, ph = np.abs(H), np.degrees(np.unwrap(np.angle(H)))
    out = dict(wc=np.nan, PM=np.nan, wpc=np.nan, GM=np.inf, DM=np.nan)
    i = np.where((mag[:-1] >= 1) & (mag[1:] < 1))[0]
    if len(i):
        i = i[0]
        out['wc'] = np.interp(0, -np.log(mag[i:i + 2]), w[i:i + 2])
        out['PM'] = 180 + np.interp(out['wc'], w[i:i + 2], ph[i:i + 2])
        out['DM'] = np.radians(out['PM']) / out['wc']
    j = np.where((ph[:-1] > -180) & (ph[1:] <= -180))[0]
    if len(j):
        j = j[0]
        out['wpc'] = np.interp(-180, ph[j:j + 2][::-1], w[j:j + 2][::-1])
        out['GM'] = 1 / np.interp(out['wpc'], w[j:j + 2], mag[j:j + 2])
    return out


def step_metrics(t, y, r_final, y_start, t0=0.0):
    """Percent overshoot and 2 % settling time (s, from t0) of a set-point step from y_start to r_final."""
    m = t >= t0
    tt, yy = t[m] - t0, y[m]
    span = r_final - y_start
    OS = max(0.0, 100 * (np.max((yy - y_start) / span) - 1))
    out = np.abs(yy - r_final) > 0.02 * abs(span)
    ts = tt[np.where(out)[0][-1] + 1] if out.any() and np.where(out)[0][-1] + 1 < len(tt) else (
        0.0 if not out.any() else np.nan)
    return OS, ts
