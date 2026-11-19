# Author: Xu Chen, University of Washington
"""Write the simulated vehicle demo log data/02_cruise_step (sedan). Seeded and deterministic.
(The homework logs are generated separately and their sidecars carry no model parameters.)"""
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ctrlprac import data
from _model import CAR_A, TS, U_LIM, WHEEL_Q, WHEEL_N, GPS_N, drive, v_steady, linearize


def meta(car, desc, test, seed):
    lin = linearize(car)
    return dict(source='simulated', description=desc,
                true_model='m_eff dv/dt = F - 0.5 rho Cd A v^2 - Crr m g cos(th) - m g sin(th) - F_gust; '
                           'tau_e dF/dt = Fmax u/100 - F (after a %.1f s delay)' % car['delay'],
                parameters={k: v for k, v in car.items() if k != 'name'}, vehicle=car['name'],
                linearized_about_v0={k: lin[k] for k in ('b', 'K', 'tau', 'u0')},
                units=dict(m='kg', m_eff='kg', A='m^2', Fmax='N', tau_e='s', delay='s', v0='m/s',
                           b='N s/m', K='(m/s)/%', tau='s', u0='%'),
                test=test, sample_time_s=TS, wheel_speed_resolution_mps=WHEEL_Q, wheel_speed_noise_std_mps=WHEEL_N,
                gps_speed_noise_std_mps=GPS_N, gust_force_std_N=15.0, throttle_limits_pct=list(U_LIM), seed=seed)


def frame(res, with_setpoint=False):
    cols = {'time_s': np.round(res['t'], 2)}
    if with_setpoint:
        cols['setpoint_mps'] = res['r']
    cols.update(throttle_pct=np.round(res['u'], 2), speed_wheel_mps=np.round(res['v_wheel'], 3),
                speed_gps_mps=np.round(res['v_gps'], 2))
    return pd.DataFrame(cols)


def step_test(car, u0, du, t_step, T, seed):
    n = int(round(T / TS)) + 1
    u = np.where(np.arange(n) * TS < t_step, u0, u0 + du)
    return drive(car, T, u_seq=u, v_init=v_steady(car, u0), seed=seed)


if __name__ == '__main__':
    # demo: sedan cruising on a flat road at 28 % throttle, throttle step +5 % at t = 20 s
    res = step_test(CAR_A, 28.0, 5.0, 20.0, 360.0, seed=301)
    data.save(frame(res), '02_cruise_step', meta(CAR_A, 'Throttle step on a flat road, sedan (demo).',
                                                  dict(throttle_pct=[28, 33], t_step_s=20, road='flat'), 301))
    print('wrote data/02_cruise_step')
