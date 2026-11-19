# Author: Xu Chen, University of Washington
"""The identification tools recover known models from simulated data."""
import numpy as np
import control as ct
from ctrlprac import ident, sim


def test_fopdt():
    t, y, _ = sim.open_loop(ct.tf([0.8], [120, 1]), 1.0, np.r_[np.zeros(20), 50 * np.ones(800)],
                            noise_std=0.1, delay=8.0, y0=21.0)
    f = ident.fit_fopdt(t, y, 50, t_step=20)
    assert abs(f['K'] - 0.8) < 0.02 and abs(f['tau'] - 120) < 6 and abs(f['theta'] - 8) < 2


def test_second_order_and_log_decrement():
    G = ct.tf([2 * 25], [1, 2 * 0.15 * 5, 25])
    t, y = ct.step_response(G, np.linspace(0, 10, 2001))
    f = ident.fit_second_order_step(t, y, 1.0)
    assert abs(f['zeta'] - 0.15) < 0.01 and abs(f['wn'] - 5) < 0.05 and abs(f['K'] - 2) < 0.02
    t, y = ct.initial_response(ct.ss(G), np.linspace(0, 10, 4001), X0=[1, 0])
    ld = ident.log_decrement(t, y)
    assert abs(ld['zeta'] - 0.15) < 0.01 and abs(ld['wn'] - 5) < 0.1


def test_etfe_and_fit_tf():
    G = ct.tf([400], [1, 4, 400]) * ct.tf([1], [0.1, 1])
    rng = np.random.default_rng(1)
    Ts, n = 0.002, 2 ** 16
    u = rng.standard_normal(n)
    t = np.arange(n) * Ts
    _, y = ct.forced_response(ct.c2d(G, Ts, 'zoh'), t, u)
    w, H, coh = ident.etfe(t, u, y, nperseg=4096)
    k = (w > 1) & (w < 200)
    Gf = ident.fit_tf(w[k], H[k], 0, 3)
    p = np.sort_complex(ct.poles(Gf))
    assert np.allclose(np.sort_complex(ct.poles(G)), p, rtol=0.05)


def test_pid_antiwindup_tracks():
    out = sim.simulate_loop(ct.tf([1], [10, 1]), sim.PID(5, 1, 0, Ts=0.05, u_lim=(0, 3)), 0.05, 60, 2.0,
                            u_lim=(0, 3))
    assert abs(out['y'][-1] - 2) < 0.02 and out['u'].max() <= 3


def test_log_decrement_flat_crests():
    """Regression: on fast, quantized data a crest spans several samples of (nearly) equal height and used
    to be counted as two or more peaks, which halved the period and spoiled zeta."""
    fs, wn, zeta = 500.0, 3.6, 0.01
    t = np.arange(0, 60, 1 / fs)
    wd = wn * np.sqrt(1 - zeta ** 2)
    y = 2.0 * np.exp(-zeta * wn * t) * np.cos(wd * t) + 0.0015 * np.random.default_rng(4).standard_normal(len(t))
    y = 1.06e-3 * np.round(y / 1.06e-3)              # 16-bit gyro resolution
    ld = ident.log_decrement(t, y, n_peaks=30)
    assert len(ld['peaks_t']) == 30
    assert np.diff(ld['peaks_t']).min() > 0.9 * 2 * np.pi / wd
    assert abs(ld['zeta'] - zeta) < 0.001 and abs(ld['wn'] - wn) < 0.01
