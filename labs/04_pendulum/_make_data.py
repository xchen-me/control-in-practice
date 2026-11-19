# Author: Xu Chen, University of Washington
"""Write the Lab 04 SIMULATED phone-pendulum logs in phyphox 'Gyroscope (rotation rate)' CSV format.

data/04_pendulum_sample      demo sample (replace it with your own recording), labelled "source": "simulated".
(The homework backup log is generated separately and its sidecar carries no model parameters.)"""
import numpy as np
from ctrlprac import data
from _model import simulate_phyphox, g

SETUPS = {
    '04_pendulum_sample': dict(d=0.75, theta0_deg=35.0, c_quad=0.08, zeta_lin=0.006, axis='z', seed=404,
                               desc='SIMULATED phyphox gyroscope export: phone on a string, pivot to phone centre '
                                    '0.75 m, released from 35 deg, swinging in the plane of the screen. '
                                    'Replace with your own recording.'),
}

if __name__ == '__main__':
    for name, s in SETUPS.items():
        df, tr = simulate_phyphox(s['d'], np.radians(s['theta0_deg']), 180.0, s['seed'], c_quad=s['c_quad'],
                                  zeta_lin=s['zeta_lin'], axis=s['axis'])
        data.save(df, name, dict(
            source='simulated', description=s['desc'], format='phyphox CSV export (comma, decimal point)',
            model="theta'' = -(g/L_eq) sin(theta) - 2 zeta_lin wn theta' - c_quad |theta'| theta';  "
                  "L_eq = d + (h^2 + w^2)/(12 d) (phone as an h x w plate swinging in its own plane)",
            parameters=dict(g=g, d_m=s['d'], phone_h_m=0.15, phone_w_m=0.072, theta0_deg=s['theta0_deg'],
                            zeta_lin=s['zeta_lin'], c_quad_per_rad=s['c_quad'], swing_axis=s['axis'],
                            hold_before_release_s=3.0, sample_rate_Hz=100.0),
            true_values=dict(L_eq_m=tr['L_eq'], wn_small_angle_rad_s=float(tr['wn_small']),
                             T0_small_angle_s=float(tr['T0']), wn_simple_pendulum_d_rad_s=float(np.sqrt(g / s['d']))),
            nonidealities='gyro bias (-5, 1.5, -2.5) mrad/s, white noise 1.5 mrad/s rms, 16-bit quantization '
                          '(1.06 mrad/s), 4 deg axis misalignment, string twist mode 0.55 Hz, hand tremor '
                          'before release, time-stamp jitter 0.15 ms rms',
            seed=s['seed']))
