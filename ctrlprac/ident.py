# Author: Xu Chen, University of Washington
"""System identification from measured data: step tests, free decay, frequency response."""
import numpy as np
import control as ct
from scipy import optimize, signal


def fit_fopdt(t, y, du, t_step=0.0):
    """Fit a first-order-plus-dead-time step response y0 + K du (1 - exp(-(t - t_step - theta)/tau)).

    Returns dict(K, tau, theta, y0, rmse, model) where model is K e^{-theta s}/(tau s + 1)
    without the delay (use ct.pade for the delay)."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    y0 = np.mean(y[t <= t_step]) if np.any(t <= t_step) else y[0]
    span = t[-1] - t_step

    def model(p):
        K, tau, theta = p
        tt = np.clip(t - t_step - theta, 0, None)
        return y0 + K * du * (1 - np.exp(-tt / tau))

    K0 = (np.mean(y[-max(5, len(y) // 20):]) - y0) / du
    res = optimize.least_squares(lambda p: model(p) - y, [K0, span / 5, span / 50],
                                 bounds=([-np.inf, 1e-6, 0], [np.inf, np.inf, span / 2]))
    K, tau, theta = res.x
    return dict(K=K, tau=tau, theta=theta, y0=y0,
                rmse=float(np.sqrt(np.mean(res.fun ** 2))), model=ct.tf([K], [tau, 1]))


def fit_second_order_step(t, y, du, t_step=0.0):
    """Fit y0 + du * step of K wn^2/(s^2 + 2 zeta wn s + wn^2). Returns dict(K, zeta, wn, rmse, model)."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    y0 = np.mean(y[t <= t_step]) if np.any(t <= t_step) else y[0]
    tt = t - t_step
    m = tt >= 0

    def model(p):
        K, z, wn = p
        out = np.full_like(t, y0)
        _, s = ct.step_response(ct.tf([K * wn ** 2], [1, 2 * z * wn, wn ** 2]), tt[m])
        out[m] = y0 + du * s
        return out

    yss = np.mean(y[-max(5, len(y) // 20):])
    K0 = (yss - y0) / du
    ipk = np.argmax(np.abs(y[m] - y0))
    os_ = min(max(abs(y[m][ipk] - y0) / max(abs(yss - y0), 1e-12) - 1, 0.01), 0.95)
    z0 = -np.log(os_) / np.sqrt(np.pi ** 2 + np.log(os_) ** 2)
    tp = max(tt[m][ipk], tt[m][1])
    wn0 = np.pi / (tp * np.sqrt(1 - z0 ** 2))
    res = optimize.least_squares(lambda p: model(p) - y, [K0, z0, wn0],
                                 bounds=([-np.inf, 1e-3, 1e-6], [np.inf, 5, np.inf]))
    K, z, wn = res.x
    return dict(K=K, zeta=z, wn=wn, rmse=float(np.sqrt(np.mean(res.fun ** 2))),
                model=ct.tf([K * wn ** 2], [1, 2 * z * wn, wn ** 2]))


def log_decrement(t, y, n_peaks=None):
    """Damping from a free decay: returns dict(zeta, wd, wn, peaks_t, peaks_y).

    Uses the positive peaks of y - mean(y); n_peaks limits how many are used. On fast, quantized data a
    crest can be flat or split into several samples of nearly equal height, so the peaks are kept at least
    half an oscillation period apart (the period comes from the spectrum of the record), and each peak time
    is the centre of the samples at the crest value."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    x = y - np.mean(y[len(y) // 2:])
    dt = float(np.median(np.diff(t)))
    nfft = 4 * 2 ** int(np.ceil(np.log2(len(x))))      # zero padding refines the spectral peak
    spec = np.abs(np.fft.rfft(x - np.mean(x), nfft))
    period = 1 / np.fft.rfftfreq(nfft, dt)[1 + np.argmax(spec[1:])]
    half = max(1, int(0.5 * period / dt))
    # prominence rejects noise wiggles on top of the oscillation peaks; distance merges split crests
    idx, _ = signal.find_peaks(x, prominence=0.05 * np.max(x), distance=half)
    if n_peaks:
        idx = idx[:n_peaks]
    if len(idx) < 3:
        raise ValueError('need at least 3 peaks for a log decrement')
    tp = t[idx].copy()
    for i, j in enumerate(idx):                          # centre of a flat or split crest
        lo, hi = max(0, j - half // 2), min(len(x), j + half // 2 + 1)
        top = lo + np.nonzero(x[lo:hi] == x[j])[0]
        tp[i] = 0.5 * (t[top[0]] + t[top[-1]])
    # slope of log amplitude per cycle gives delta = 2 pi zeta / sqrt(1 - zeta^2)
    n = np.arange(len(idx))
    delta = -np.polyfit(n, np.log(x[idx]), 1)[0]
    zeta = delta / np.sqrt(4 * np.pi ** 2 + delta ** 2)
    wd = 2 * np.pi / np.mean(np.diff(tp))
    return dict(zeta=zeta, wd=wd, wn=wd / np.sqrt(1 - zeta ** 2), peaks_t=tp, peaks_y=y[idx])


def etfe(t, u, y, nperseg=None):
    """Empirical frequency response H(jw) = S_uy / S_uu (Welch averages) and coherence.

    Returns w (rad/s), H (complex), coh. Uniform sampling is assumed."""
    t, u, y = (np.asarray(a, float) for a in (t, u, y))
    fs = 1.0 / np.median(np.diff(t))
    nperseg = nperseg or min(len(t), 2 ** int(np.log2(len(t) / 4)))
    f, Suu = signal.welch(u - u.mean(), fs, nperseg=nperseg)
    _, Suy = signal.csd(u - u.mean(), y - y.mean(), fs, nperseg=nperseg)
    _, coh = signal.coherence(u, y, fs, nperseg=nperseg)
    k = f > 0
    return 2 * np.pi * f[k], Suy[k] / Suu[k], coh[k]


def sine_response(t, u, y, w):
    """Gain and phase (deg) of y relative to u at one frequency w (rad/s), by projection on sin/cos.

    Use a window that holds an integer number of periods in steady state."""
    t, u, y = (np.asarray(a, float) for a in (t, u, y))
    basis = np.exp(-1j * w * t)
    U, Y = np.sum((u - u.mean()) * basis), np.sum((y - y.mean()) * basis)
    H = Y / U
    return abs(H), np.degrees(np.angle(H))


def fit_tf(w, H, nb, na, weight=None, iters=10, delay=0.0):
    """Continuous-time transfer function b(s)/a(s) (deg b = nb, deg a = na, monic a) fitted to H(jw).

    Sanathanan-Koerner iterations of a linear least-squares (Levy) fit. A known delay is removed first."""
    w, H = np.asarray(w, float), np.asarray(H, complex) * np.exp(1j * np.asarray(w) * delay)
    s = 1j * w
    W = np.ones_like(w) if weight is None else np.asarray(weight, float)
    a_prev = np.ones_like(s)
    for _ in range(iters):
        # unknowns: b_nb..b_0, a_{na-1}..a_0 ; equation b(s) - H a'(s) = H s^na
        A = np.hstack([s[:, None] ** np.arange(nb, -1, -1),
                       -H[:, None] * s[:, None] ** np.arange(na - 1, -1, -1)])
        rhs = H * s ** na
        Wk = (W / np.abs(a_prev))[:, None]
        M = np.vstack([(Wk * A).real, (Wk * A).imag])
        v = np.concatenate([(Wk[:, 0] * rhs).real, (Wk[:, 0] * rhs).imag])
        x = np.linalg.lstsq(M, v, rcond=None)[0]
        b, a = x[:nb + 1], np.concatenate([[1.0], x[nb + 1:]])
        a_prev = np.polyval(a, s)
    return ct.tf(b, a)
