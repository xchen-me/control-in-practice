# Author: Xu Chen, University of Washington
# %% [markdown]
# # Homework starter, Lab 05: balancing robot R2
#
# ![scene](figures/scene_hw.png)
#
# Loads the three R2 logs and plots them. Fill in the TODO steps (see homework.md).

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import control as ct
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
from ctrlprac import data, sim
from _model import TS, pd_frf, winding_ccw     # exact firmware-PD frequency response; encirclement count

drop, mdrop = data.load('05_hw_drop')
push, mpush = data.load('05_hw_push')
tip, mtip = data.load('05_hw_tipover')
Kp, Kd, Tf = (mpush['test'][k] for k in ('Kp_Nm_per_rad', 'Kd_Nms_per_rad', 'Tf_s'))
print('team gains: Kp = %.3f N m/rad, Kd = %.3f N m s/rad, Tf = %.3f s; loop delay %d sample(s); '
      'tip-over log delay %d samples' % (Kp, Kd, Tf, mpush['test']['loop_delay_samples'],
                                         mtip['test']['loop_delay_samples']))

fig, ax = plt.subplots(3, 1, figsize=(9, 8))
ax[0].plot(drop['time_s'], drop['theta_deg'], '.', ms=3)
ax[0].set(xlabel='time (s)', ylabel='theta (deg)', title='05_hw_drop: motors off')
ax[1].plot(push['time_s'], push['theta_deg'], lw=0.8, label='theta (deg)')
ax[1].plot(push['time_s'], 100 * push['tau_cmd_Nm'], lw=0.5, alpha=0.7, label='100 x tau_cmd (N m)')
ax[1].set(xlabel='time (s)', ylabel='theta (deg) / scaled torque', title='05_hw_push: balancing, two pushes')
ax[2].plot(tip['time_s'], tip['theta_deg'], lw=0.8)
ax[2].set(xlabel='time (s)', ylabel='theta (deg)', title='05_hw_tipover: firmware v2')
for a in ax:
    a.grid(True, alpha=0.4)
ax[1].legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 1. Drop test
# TODO: least-squares fit theta = A e^{pt} + B e^{-pt} (+ offset) on |theta| < 10 deg, scanning p.
p_drop = None

# %% [markdown]
# ## 2. Closed-loop identification of p^2 and b0
# TODO: integrated equation-error fit on push-free windows; torque shifted by the loop delay.
p_id, b0_id = None, None

# %% [markdown]
# ## 3-4. Routh range and PD design on the root locus
# TODO: Routh conditions; Kp, Kd for zeta = 0.7, omega_n = 2p; root locus plot.

# %% [markdown]
# ## 5. Bode plot of G(s) = b0/(s^2 - p^2) and G_s(s) = b0/(s + p)^2
# TODO

# %% [markdown]
# ## 6. Nyquist, gain margins (lower and upper), phase margin, delay limit
# L(jw) = pd_frf(Kp, Kd, Tf, TS, w) * G(jw) * exp(-1j w (n + 0.5) TS)
# TODO: encirclements (winding_ccw), GMs, PM, max delay; check with sim.simulate_loop(..., delay=(n-1)*TS)

# %% [markdown]
# ## 7. Tip-over log
# TODO: oscillation frequency from the log vs. unstable closed-loop pole of your model.
