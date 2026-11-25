# Author: Xu Chen, University of Washington
"""Write the simulated demo logs of drive D1 (seeded and deterministic): data/07_hdd_frf (swept-sine
frequency response of the VCM plant, as exported by a servo analyzer) and data/07_hdd_pes (PES log during
track following). (The homework logs are generated separately and their sidecars carry no model parameters.)
Both tests run with the drive's firmware track-following controller C0 in the loop."""
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
from ctrlprac import data
import _model as M

# firmware track-following controller C0 of each drive: lead x integral x notch (see M.controller)
FIRMWARE = {'D1': dict(K=0.001426, fz=210.0, fp=3049.0, fi=80.0, notch=(6200.0, 20.0, 0.2))}


def firmware(dr):
    return M.controller(M.sample_time(dr), **FIRMWARE[dr['name']])[0]


def noise_psd(dr, f):
    """One-sided PSD (track^2/Hz) at f of the output-equivalent noise seen by the analyzer: PES noise, disk
    flutter and windage torque through the plant (RRO is synchronous and subtracted before the correlation)."""
    Ts = M.sample_time(dr)
    fs = 1 / Ts
    psd = np.full(len(f), 2 * dr['noise_track'] ** 2 / fs)
    fg = np.linspace(1, fs / 2, 200000)
    for f0, z, rms in dr['flutter']:
        H = lambda x: np.abs(2 * z * (1j * x / f0) / (1 - (x / f0) ** 2 + 2j * z * x / f0)) ** 2
        psd += rms ** 2 * H(f) / np.trapezoid(H(fg), fg)
    Hw = lambda x: 1 / (1 + (x / dr['windage_fc']) ** 2)
    G = np.abs(M.frf(M.plant_dt(dr), f)) ** 2
    psd += dr['windage_A'] ** 2 * Hw(f) / np.trapezoid(Hw(fg), fg) * G
    return psd


def analyzer_frf(dr, f, amp, n_avg, t_min, cycles, seed):
    """Swept sine with a servo analyzer, excitation r added to the controller output (closed loop):
    y = G u + D, u = r - C y. With the excitation as reference the estimate over n_avg blocks is
    G_hat = (G a + mean D)/(a - C mean D); each block's noise amplitude has E|D|^2 = 2 PSD/T (T = block
    length). The coherence is |sum Y|^2 / (n sum |Y|^2) between excitation and response."""
    rng = np.random.default_rng(seed)
    G = M.frf(M.plant_dt(dr), f)
    C = M.frf(firmware(dr), f)
    T = np.maximum(t_min, cycles / f)
    sd = np.sqrt(2 * noise_psd(dr, f) / T / 2)            # per real and imaginary part
    D = sd[:, None] * (rng.standard_normal((len(f), n_avg)) + 1j * rng.standard_normal((len(f), n_avg)))
    a = amp(f)[:, None]
    Y = G[:, None] * a + D                                 # common factor 1/(1 + G C) cancels in both ratios
    U = a - C[:, None] * D
    Gh = Y.sum(1) / U.sum(1)
    coh = np.abs(Y.sum(1)) ** 2 / (n_avg * (np.abs(Y) ** 2).sum(1))
    return Gh, coh, T


def frf_frame(f, Gh, coh):
    return pd.DataFrame({'f_Hz': np.round(f, 2), 'mag_dB': np.round(20 * np.log10(np.abs(Gh)), 3),
                         'phase_deg': np.round(np.degrees(np.angle(Gh)), 2), 'coherence': np.round(coh, 4)})


def pes_frame(dr, res):
    n = len(res['pes'])
    k = np.arange(n)
    return pd.DataFrame({'time_s': np.round(res['t'], 7), 'rev': k // dr['sectors'], 'sector': k % dr['sectors'],
                         'pes_pct': np.round(100 * res['pes'], 4), 'u_mA': np.round(1e3 * res['u'], 4)})


def drive_meta(dr):
    """Drive constants shared by the sidecars."""
    return dict(drive=dr['name'], spindle_rpm=dr['rpm'], servo_sectors=dr['sectors'],
                sample_time_s=M.sample_time(dr), track_pitch_nm=dr['track_pitch_nm'],
                pes_resolution_pct=round(100 * dr['pes_res_track'], 4),
                firmware_controller=dict(FIRMWARE[dr['name']], form='C(z) = lead K (s/wz + 1)/(s/wp + 1) x integral '
                                         '(s + wi)/s x notch (f0, depth_dB, zeta_p); corner frequencies in Hz; Tustin, '
                                         'prewarped (labs/07_hdd/_model.py controller())', units='K in A/track'))


def truth_meta(dr):
    p = {k: dr[k] for k in ('Kt', 'l_arm', 'J', 'tau_s', 'modes')}
    return dict(true_model='y/u = Kp [1/s^2 + sum r_i/(s^2 + 2 zeta_i w_i s + w_i^2)] (tracks/A), sampled with '
                           'a zero-order hold and a computation delay tau; see labs/07_hdd/_model.py',
                parameters=p, units=dict(Kt='N m/A', l_arm='m', J='kg m^2', tau_s='s', f='Hz'),
                true_derived=dict(Kp_track_per_As2=round(M.rigid_gain(dr), 1),
                                  Td_us=round(1e6 * (M.sample_time(dr) / 2 + dr['tau_s']), 3)),
                disturbances={k: dr[k] for k in ('rro_track', 'rro_wi_track', 'flutter', 'windage_A',
                                                  'windage_fc', 'noise_track', 'bias_A')},
                disturbance_units=dict(rro_track='track (1x, 2x, ...)', rro_wi_track='track at 6x, ~1/sqrt(k)',
                                       flutter='(f Hz, zeta, rms track)', windage_A='A rms (current-equivalent)',
                                       windage_fc='Hz', noise_track='track rms', bias_A='A'))


FRF_TEST = dict(points=400, f_lo=100.0, f_hi=17000.0, n_avg=8, t_min=0.004, cycles=10, amp_mA_1k=0.6, amp_max_mA=8.0)


def make_frf(dr, seed, t=None):
    t = t or FRF_TEST
    f = np.logspace(np.log10(t['f_lo']), np.log10(t['f_hi']), t['points'])
    amp = lambda x: 1e-3 * np.minimum(t['amp_mA_1k'] * np.maximum(1.0, (x / 1000) ** 2), t['amp_max_mA'])
    Gh, coh, T = analyzer_frf(dr, f, amp, t['n_avg'], t['t_min'], t['cycles'], seed)
    test = dict(excitation='stepped swept sine added to the controller output (track following, firmware '
                           'controller in the loop)', points=t['points'], f_range_Hz=[t['f_lo'], t['f_hi']],
                amplitude='%g mA up to 1 kHz, then rising as f^2, at most %g mA' % (t['amp_mA_1k'], t['amp_max_mA']),
                averages=t['n_avg'], integration='max(%g ms, %d periods) per average' % (1e3 * t['t_min'], t['cycles']),
                estimate='excitation as reference: G = H(exc -> position)/H(exc -> current); synchronous '
                         'runout removed before the correlation; coherence between excitation and position')
    return frf_frame(f, Gh, coh), test


PES_TEST = dict(revolutions=50)


def make_pes(dr, seed):
    res = M.simulate_pes(dr, firmware(dr), PES_TEST['revolutions'], seed, settle=False)   # as the log was written
    test = dict(revolutions=PES_TEST['revolutions'], note='track following on one track with the firmware '
                'controller; PES quantized to the PES resolution; u is the VCM current command')
    return pes_frame(dr, res), test, res


def make_all():
    dr, seeds = M.DRIVE_A, (701, 702)
    df, test = make_frf(dr, seeds[0])
    meta = dict(source='simulated', description='Swept-sine frequency response of the VCM actuator of drive D1, '
                'VCM current to head position, exported by the servo analyzer (demo).', **drive_meta(dr),
                **truth_meta(dr), test=test, seed=seeds[0],
                columns=dict(f_Hz='Hz', mag_dB='dB re 1 track/A', phase_deg='deg, wrapped to +/-180',
                             coherence='-'))
    data.save(df, '07_hdd_frf', meta)
    dp, test, res = make_pes(dr, seeds[1])
    meta = dict(source='simulated', description='PES log of drive D1 during track following, 50 revolutions, '
                'firmware controller (demo).', **drive_meta(dr), **truth_meta(dr), test=test, seed=seeds[1],
                columns=dict(time_s='s', rev='-', sector='-', pes_pct='% of track pitch', u_mA='mA'))
    data.save(dp, '07_hdd_pes', meta)
    return df, dp


if __name__ == '__main__':
    df, dp = make_all()
    print('07_hdd_frf: %d points, coherence min %.3f; 07_hdd_pes: %d sectors, PES 3 sigma %.1f %%'
          % (len(df), df['coherence'].min(), len(dp), 3 * dp['pes_pct'].std()))
    print('wrote data/07_hdd_*')
