# Author: Xu Chen, University of Washington
"""Write the simulated balancing-robot demo logs (200 Hz firmware logs, seeded and deterministic):
  robot R1: data/05_balance_drop, data/05_balance_push, data/05_balance_tipover
(The homework logs are generated separately and their sidecars carry no model parameters.)"""
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
from ctrlprac import data, sim
import _model as M


def gains_A():
    """R1 firmware gains: PD placed at zeta = 0.7, omega_n = 2p (designed in demo_balance_2)."""
    p, b0 = M.pitch_params(M.ROBOT_A)
    wn = 2 * p
    return (wn ** 2 + p ** 2) / b0, 2 * 0.7 * wn / b0



def to_df(o, dither=None):
    df = pd.DataFrame({'time_s': np.round(o['t'], 4),
                         'theta_deg': np.round(np.degrees(o['theta_hat']) / M.LOG_RES_DEG) * M.LOG_RES_DEG,
                         'gyro_dps': np.round(np.degrees(o['gyro']), 2),
                         'acc_tilt_deg': np.round(np.degrees(o['acc_tilt']), 2),
                         'tau_cmd_Nm': np.round(o['tau_cmd'], 4)})
    if dither is not None:
        df['dither_Nm'] = np.round(dither[:len(df)], 4)
    return df


class Dithered:
    """Firmware PD plus a PRBS identification dither (amplitude amp N m, one bit = hold samples)."""

    def __init__(self, pid, amp, hold, n, seed):
        self.pid, self.k = pid, 0
        bits = np.sign(np.random.default_rng(seed).standard_normal(n // hold + 1))
        self.d = amp * np.repeat(bits, hold)[:n]

    def __call__(self, r, y):
        u = self.pid(r, y) + self.d[self.k]
        self.k += 1
        return u


def meta(rb, desc, test, seed):
    p, b0 = M.pitch_params(rb)
    return dict(source='simulated', description=desc,
                true_model='theta/tau = b0/(s^2 - p^2) (linearized); simulated with the nonlinear '
                           'wheel + body equations (sin theta), RK4, zero-order-hold torque',
                true_p_rad_s=round(float(p), 4), true_b0_rad_s2_per_Nm=round(float(b0), 3),
                parameters={k: v for k, v in rb.items() if k != 'name'},
                units=dict(m_b='kg', l='m (axle to COM)', I_b='kg m^2 (about COM)', m_w='kg (both wheels)',
                           r='m (wheel radius)', tau_max='N m (both motors)', h_push='m',
                           n_delay='samples of sensor + computation delay'),
                robot=rb['name'], sample_time_s=M.TS, imu=M.IMU,
                sign_convention='theta > 0 when the body leans backward; tau > 0 drives the wheels forward',
                test=test, seed=seed)


def push(F, t0, dur=0.1, F2=None, t2=None):
    return lambda t: (F if t0 <= t < t0 + dur else 0.0) + (F2 if F2 and t2 <= t < t2 + dur else 0.0)


def make(rb, prefix, gains, n_bad, seeds, tag):
    Kp, Kd = gains
    Tf = M.TS
    # 1) drop test: motors in coast mode, released by hand near upright, caught by a stand at 30 deg
    o = M.simulate(rb, None, 3.0, theta0_deg=0.8, motors_on=False, stop_deg=30.0, seed=seeds[0])
    data.save(to_df(o), prefix + '_drop',
              meta(rb, 'Drop test on robot %s (%s): motors off (coast), body released near upright, logged '
                   'until it reaches the 30 deg stand.' % (rb['name'], tag),
                   dict(motors='off', release_deg=0.8, stop_deg=30), seeds[0]))
    # 2) balancing with two pushes, nominal firmware (one sample of delay)
    pf = push(4.0 * rb['m_b'], 3.0, F2=-3.0 * rb['m_b'], t2=6.5)
    amp = 0.1 * rb['tau_max']
    ctrl = Dithered(sim.PID(Kp, Kd=Kd, Ts=M.TS, Tf=Tf), amp, 4, 2001, seeds[1] + 1000)
    o = M.simulate(rb, ctrl, 10.0, push=pf, theta0_deg=1.0, seed=seeds[1])
    data.save(to_df(o, ctrl.d), prefix + '_push',
              meta(rb, 'Robot %s (%s) balancing under PD control; pushed forward at the top at t = 3 s and '
                   'backward at t = 6.5 s (push force not logged).' % (rb['name'], tag),
                   dict(controller='PD, tau = -(Kp theta + Kd dtheta/dt filtered), derivative filter Tf',
                        Kp_Nm_per_rad=round(Kp, 4), Kd_Nms_per_rad=round(Kd, 5), Tf_s=Tf,
                        loop_delay_samples=rb['n_delay'], push_N=[4.0 * rb['m_b'], -3.0 * rb['m_b']],
                        push_t_s=[3.0, 6.5], push_duration_s=0.1, dither_Nm=amp, dither_bit_s=4 * M.TS),
                   seeds[1]))
    # 3) same gains, firmware v2 adds a blocking telemetry call: n_bad samples of delay; it tips over
    o = M.simulate(rb, sim.PID(Kp, Kd=Kd, Ts=M.TS, Tf=Tf), 10.0,
                   theta0_deg=0.1, n_delay=n_bad, seed=seeds[2])
    fell = abs(np.degrees(o['theta'][-1])) > 44
    data.save(to_df(o), prefix + '_tipover',
              meta(rb, 'Robot %s (%s), same PD gains, firmware v2 with blocking telemetry: loop delay %d '
                   'samples, started from upright with no push. The robot tips over; the log ends when the 45 deg tilt cutoff trips.'
                   % (rb['name'], tag, n_bad),
                   dict(Kp_Nm_per_rad=round(Kp, 4), Kd_Nms_per_rad=round(Kd, 5), Tf_s=Tf,
                        loop_delay_samples=n_bad, loop_delay_ms=n_bad * M.TS * 1e3, fell=bool(fell),
                        t_fall_s=float(o['t'][-1]) if fell else None, release_deg=0.1), seeds[2]))
    return fell


if __name__ == '__main__':
    fa = make(M.ROBOT_A, '05_balance', gains_A(), 6, (501, 502, 503), 'demo')
    print('wrote data/05_balance_*; tip-over log fell:', fa)
