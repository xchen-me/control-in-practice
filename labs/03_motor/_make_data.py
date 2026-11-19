# Author: Xu Chen, University of Washington
"""Write the Lab 03 demo log: voltage step test of a gearmotor read through a quadrature encoder.

data/03_motor_step     Motor A (demos)     3 V -> 9 V
(The homework log is generated separately; its sidecar carries no hidden parameters.)
The motor runs 0.5 s before logging starts, so the log begins at a steady speed (operating point)."""
import numpy as np
import pandas as pd
from ctrlprac import data
from _model import MOTOR_A, simulate_log, counts_per_rev, full_tf, reduced_params


def make(name, p, v0, v1, t_step, T, seed, count0, desc):
    pre = 0.5
    t, V, c, th = simulate_log(p, lambda tk: v0 if tk < pre + t_step else v1, pre + T, seed, count0=count0)
    k = t >= pre - 1e-9
    df = pd.DataFrame({'t (s)': np.round(t[k] - pre, 6), 'V_cmd (V)': V[k], 'encoder (counts)': c[k]})
    K, tau = reduced_params(p)
    poles = np.sort(np.roots(full_tf(p).den[0][0]).real)
    data.save(df, name, dict(
        source='simulated', description=desc,
        model='i: L di/dt = V(t - Ts) - R i - kt w_m;  w_m: J dw_m/dt = kt i - b w_m - Tc tanh(w_m/0.5);'
              '  theta_out = theta_m / N;  J = Jm + JL/N^2; encoder counts = floor(theta_out cpr/(2 pi))',
        parameters=p, counts_per_rev=counts_per_rev(p), sample_time_s=p['Ts'],
        step_V=[v0, v1], step_time_s=t_step,
        true_reduced_model=dict(K_rad_per_s_per_V=K, tau_s=tau,
                                poles_of_full_model=[float(x) for x in poles]),
        nonidealities='encoder quantization, 1-sample PWM update delay, 0.5 %% PWM amplitude ripple, '
                      'Coulomb friction (equivalent %.2f V), electrical lag L/R' % (p['Tc'] * p['R'] / p['kt']),
        seed=seed))


if __name__ == '__main__':
    make('03_motor_step', MOTOR_A, 3.0, 9.0, 0.3, 1.0, 303, 51234,
         'Motor A (12 V, 30:1 gearmotor, 64 CPR encoder) voltage step 3 V -> 9 V, encoder logged at 1 kHz')
