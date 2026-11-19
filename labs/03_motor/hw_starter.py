# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 03 homework starter: Motor B
#
# ![scene](figures/scene_hw.png)
#
# Loads the logged voltage step, plots the raw log, and leaves the steps of homework.md as TODOs.

# %%
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.modules.pop('_model', None)          # run_all.py runs several labs in one process
import numpy as np
import matplotlib.pyplot as plt
import control as ct
from ctrlprac import data, ident
from ctrlprac.sim import PID, LTIController, simulate_loop
from _model import MOTOR_B, full_tf, counts_per_rev, speed_from_counts, pd_by_angle, lead_by_angle, step_specs

p = MOTOR_B
df, meta = data.load('03_motor_hw_step')
t, V, c = df['t (s)'].values, df['V_cmd (V)'].values, df['encoder (counts)'].values
cpr = counts_per_rev(p)
print('counts per output revolution:', cpr, '-> one count =', 2 * np.pi / cpr, 'rad')

fig, ax = plt.subplots(2, 1, figsize=(8, 5), sharex=True)
ax[0].plot(t, V, 'k')
ax[0].set(ylabel='V_cmd (V)', title='Motor B: raw log')
ax[1].plot(t, c)
ax[1].set(xlabel='time (s)', ylabel='encoder (counts)')
for a in ax:
    a.grid(True)
fig.tight_layout()
plt.show()

# %% Part 1: model and dominant pole
G3 = full_tf(p)
print('G3(s) =', G3)
# TODO: poles of G3, electrical/mechanical ratio, K and tau from the parameter formulas

# %% Part 2: identification
w1 = speed_from_counts(t, c, cpr, window=1, causal=True)
# TODO: 10-sample speed, resolution of each, fit with ident.fit_fopdt on the increment, plot, errors,
#       friction voltage from the operating point

# %% Part 3: P control and the root locus
# TODO: spec region (zeta_min, sigma_min), ct.root_locus_map, largest Kp for the overshoot spec

# %% Part 4: PD design
# TODO: pd_by_angle(K, tau, sd) for zeta = 0.70, 0.80, 0.85 and sigma = 28; responses with and without the zero

# %% Part 5: lead and the encoder
q, Ts, ulim = 2 * np.pi / cpr, p['Ts'], (-p['V_max'], p['V_max'])
# TODO: lead_by_angle(K, tau, sd, 1/tau); simulate_loop(G3, LTIController(C, Ts), Ts, 0.5, 0.3,
#       u_lim=ulim, quant=q, delay=Ts) and the same with PID(Kp, 0, Kd, Ts, Tf=Ts/50)

# %% Part 6: load torque
# TODO: equivalent voltage of 0.8 N m, lag and PI designs, simulate_loop(..., d=lambda t: ...)

# %% Part 7: engineering judgement (write it in your report)
