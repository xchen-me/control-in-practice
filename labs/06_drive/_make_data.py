# Author: Xu Chen, University of Washington
"""Write the simulated belt-stage demo identification logs (stage S1, 4 kHz drive logs, seeded and
deterministic): data/06_drive_prbs (PRBS), data/06_drive_steppedsine. (The homework logs are generated
separately and their sidecars carry no model parameters.)
All tests run inside a low-bandwidth 'parking' PD loop (a double integrator can not be tested open loop);
the excitation current i_exc is added to the controller output in the drive."""
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
from ctrlprac import data, sim
import _model as M

PARK = {'S1': (0.02, 0.003, 0.004)}   # Kp (A/mm), Kd (A s/mm), Tf (s)


class Excited:
    """Parking PD plus a logged excitation sequence."""

    def __init__(self, pid, exc):
        self.pid, self.exc, self.k = pid, exc, 0

    def __call__(self, r, y):
        u = self.pid(r, y) + self.exc[min(self.k, len(self.exc) - 1)]
        self.k += 1
        return u


def run(st, exc, seed):
    Kp, Kd, Tf = PARK[st['name']]
    rng = np.random.default_rng(seed + 7)
    noise = M.I_NOISE_A * rng.standard_normal(len(exc) + 1)
    ctrl = Excited(sim.PID(Kp, Kd=Kd, Ts=M.TS, Tf=Tf, u_lim=M.I_LIM), exc)
    res = sim.simulate_loop(M.plant_ss(st), ctrl, M.TS, (len(exc) - 1) * M.TS, 0.0,
                            d=lambda t: noise[int(round(t / M.TS))], u_lim=M.I_LIM,
                            noise_std=M.ENC_NOISE_MM, quant=M.ENC_RES_MM,
                            delay=(M.N_DELAY - 1) * M.TS, seed=seed)
    return res


def meta(st, desc, test, seed):
    d = M.derived(st)
    Kp, Kd, Tf = PARK[st['name']]
    return dict(source='simulated', description=desc,
                true_model='x/i = current loop wc/(s+wc) * two-mass belt drive (non-collocated encoder), '
                           'see labs/06_drive/_model.py; plus %d samples of delay and ZOH' % M.N_DELAY,
                parameters={k: v for k, v in st.items() if k != 'name'},
                units=dict(J_m='kg m^2', r='m', m_L='kg', k_belt='N/m', b_L='N s/m', Kt='N m/A', f_cl='Hz'),
                true_derived=dict(f_r_Hz=round(d['f_r'], 3), zeta_r=st['zeta_belt'],
                                  J_total_kgm2=d['J'], K_rigid_mm_per_As2=round(d['K_rigid'], 2),
                                  a_rigid_rad_s=round(d['a_rigid'], 4)),
                stage=st['name'], sample_time_s=M.TS, delay_samples=M.N_DELAY,
                parking_loop=dict(Kp_A_per_mm=Kp, Kd_As_per_mm=Kd, Tf_s=Tf),
                encoder_resolution_mm=M.ENC_RES_MM, encoder_noise_mm=M.ENC_NOISE_MM,
                current_noise_A=M.I_NOISE_A, current_limit_A=list(M.I_LIM), test=test, seed=seed)


def frame(res, exc, extra=None):
    df = pd.DataFrame({'time_s': np.round(res['t'], 6), 'i_exc_A': np.round(exc, 4),
                       'i_cmd_A': np.round(res['u'], 4), 'x_mm': np.round(res['y_meas'], 4)})
    if extra is not None:
        for k, v in extra.items():
            df[k] = v
    return df


def prbs(n, amp, hold, seed):
    bits = np.sign(np.random.default_rng(seed).standard_normal(n // hold + 1))
    return amp * np.repeat(bits, hold)[:n]


def stepped(freqs, amp_of_f, settle_of_f, min_meas=0.4):
    """Each segment: amplitude faded in over 0.3 s (half cosine), held until the settle time, then an
    integer number of periods (at least 3 and at least min_meas s) flagged measure = 1."""
    sig, fcol, seg, meas = [], [], [], []
    for i, f in enumerate(freqs):
        settle = settle_of_f(f)
        n_set = int(round(settle * M.FS))
        n_per = max(3, int(np.ceil(min_meas * f)))
        n = n_set + int(round(n_per / f * M.FS))
        t = np.arange(n) * M.TS
        env = np.where(t < 0.3, 0.5 * (1 - np.cos(np.pi * t / 0.3)), 1.0)
        sig.append(amp_of_f(f) * env * np.sin(2 * np.pi * f * t))
        fcol.append(np.full(n, f))
        seg.append(np.full(n, i))
        meas.append((np.arange(n) >= n_set).astype(int))
    return np.concatenate(sig), np.concatenate(fcol), np.concatenate(seg), np.concatenate(meas)


def make_all():
    out = {}
    st, prefix, tag, seeds = M.STAGE_A, '06_drive', 'demo', (601, 602)
    exc = prbs(int(8 * M.FS) + 1, 0.4, 4, seeds[0] + 100)
    test = dict(excitation='PRBS', amplitude_A=0.4, bit_s=4 * M.TS, duration_s=8)
    desc = 'PRBS current excitation (+/-0.4 A, 1 ms bits) on stage %s inside the parking loop.'
    res = run(st, exc, seeds[0])
    data.save(frame(res, exc), prefix + '_prbs', meta(st, desc % st['name'] + ' (%s)' % tag, test, seeds[0]))
    out[prefix + '_prbs'] = res
    freqs = [2, 3, 5, 8, 12, 18, 26, 36, 44, 52, 56, 58, 60, 62, 65, 70, 80, 100, 130, 170]
    f_pk = 60.0                                        # resonance seen in the broadband test
    # amplitude grows with frequency (rigid-body roll-off) and is cut to 25 % within 5 Hz of the peak
    amp = lambda f: float(np.clip(0.15 * (f / 5) ** 0.8, 0.15, 1.5)) * (0.25 if abs(f - f_pk) < 5 else 1.0)
    settle = lambda f: 1.0 if abs(f - f_pk) < 5 else 0.3   # light damping rings for 1/(zeta w_r) s
    exc, fcol, seg, meas = stepped(freqs, amp, settle)
    res = run(st, exc, seeds[1])
    test = dict(excitation='stepped sine', frequencies_Hz=freqs,
                note='each segment: amplitude faded in over 0.3 s, settling time 0.3 s (1.0 s within 5 Hz '
                     'of the resonance), then the analysis window (measure = 1) of an integer number of periods '
                     '(at least 3, at least 0.4 s) at full amplitude; amplitude grows with frequency and is '
                     'reduced to 25 % within 5 Hz of the resonance')
    data.save(frame(res, exc, {'f_exc_Hz': fcol, 'segment': seg, 'measure': meas}), prefix + '_steppedsine',
              meta(st, 'Stepped-sine current excitation on stage %s inside the parking loop (%s).'
                   % (st['name'], tag), test, seeds[1]))
    out[prefix + '_steppedsine'] = res
    return out


if __name__ == '__main__':
    out = make_all()
    for k, r in out.items():
        print('%-22s x range %.2f..%.2f mm, |i| max %.2f A' % (k, r['y'].min(), r['y'].max(), np.abs(r['u']).max()))
    print('wrote data/06_drive_*')
