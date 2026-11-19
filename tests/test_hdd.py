# Author: Xu Chen, University of Washington
"""Lab 07 helpers: the sampled plant with delay, the margins from frequency data, the identification of the
demo drive, and the PES prediction from a log."""
import os
import sys
import numpy as np
import control as ct

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'labs', '07_hdd'))
sys.modules.pop('_model', None)
import _model as M  # noqa: E402

DR = M.DRIVE_A
TS = M.sample_time(DR)


def test_plant_delay_and_frf():
    f = np.array([100.0, 1000.0, 6200.0, 15000.0])
    P = M.plant_dt(DR)
    ref = np.array([complex(np.asarray(P(np.exp(2j * np.pi * x * TS))).ravel()[0]) for x in f])
    assert np.allclose(M.frf(P, f), ref, rtol=1e-10)
    # one more sample of delay is one more factor z^-1
    P2 = M.plant_dt(DR, tau=DR['tau_s'] + TS)
    assert np.allclose(M.frf(P2, f), ref * np.exp(-2j * np.pi * f * TS), rtol=1e-10)
    # at low frequency: Kp/s^2 with the delay Ts/2 + tau
    G = M.frf(P, [100.0])[0]
    w = 2 * np.pi * 100.0
    assert abs(G / (-M.rigid_gain(DR) / w ** 2 * np.exp(-1j * w * (TS / 2 + DR['tau_s']))) - 1) < 1e-3


def test_margins_match_python_control():
    C, _ = M.controller(TS, 0.0032, 315.0, 4570.0, 120.0, (6200.0, 20.0, 0.2))
    f = np.logspace(1, np.log10(0.4999 / TS), 20000)
    L = M.frf(C, f) * M.frf(M.plant_dt(DR), f)
    m = M.frf_margins(f, L)
    gm, pm, sm, wpc, wgc, _ = ct.stability_margins(ct.frd(L, 2 * np.pi * f))
    assert abs(m['PM'][0] - pm) < 0.3 and abs(m['fc'][0] - wgc / 2 / np.pi) < 5
    assert abs(20 * np.log10(m['GM_up'][1] / gm)) < 0.1 and abs(m['Ms'] - 1 / sm) < 0.01
    assert np.max(np.abs(M.cl_poles(M.plant_dt(DR), C))) < 1


def test_identify_demo_drive():
    f, G, coh, meta = M.load_frf('07_hdd_frf')
    fit = M.identify(f, G, coh, TS)
    assert abs(fit['Kp'] / M.rigid_gain(DR) - 1) < 0.01
    assert abs(fit['Td'] - (TS / 2 + DR['tau_s'])) < 0.5e-6
    for m, mt in zip(fit['modes'], DR['modes']):
        assert abs(m['f'] / mt['f'] - 1) < 0.01 and abs(m['zeta'] / mt['zeta'] - 1) < 0.1 and abs(m['r'] - mt['r']) < 0.03


def test_predicted_pes_matches_simulation():
    """For the linear loop (no quantization) the log-based prediction equals the simulation."""
    C0, _ = M.controller(TS, 0.001426, 210.0, 3049.0, 80.0, (6200.0, 20.0, 0.2))
    C1, _ = M.controller(TS, 0.0032, 315.0, 4570.0, 120.0, (6200.0, 20.0, 0.2))
    P = M.plant_dt(DR)
    S = lambda C: (lambda ff: 1 / (1 + M.frf(C, ff) * M.frf(P, ff)))
    p0 = M.simulate_pes(DR, C0, 20, 3, quant=False)['pes']
    p1 = M.simulate_pes(DR, C1, 20, 3, quant=False)['pes']
    pr = M.predict_pes(p0, TS, S(C0), S(C1))
    assert abs(np.std(pr) / np.std(p1) - 1) < 0.02
