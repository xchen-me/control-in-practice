# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 04 homework starter: your phone pendulum
#
# ![scene](figures/scene_hw.png)
#
# Loads YOUR recording `data/04_pendulum_my_recording.csv` (phyphox gyroscope CSV). If it is not there
# yet, it falls back to the simulated backup and says so loudly. Fill in the TODOs (see homework.md).
# This script is part of what you submit: it must run as submitted and print every number in your report.
# The web demo may be used to check your numbers, not to replace this code.

# %%
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.modules.pop('_model', None)          # run_all.py runs several labs in one process
import numpy as np
import matplotlib.pyplot as plt
from ctrlprac import ident
from _model import load_phyphox, find_release, angle_from_gyro, cycles, L_equivalent, period_large_angle, \
    data_path, g

# ---- your setup (m): edit these
d, d_unc = 0.80, 0.005                   # pivot to phone centre and its uncertainty
h, w_ph = 0.150, 0.072                   # phone height and width

path = data_path('04_pendulum_my_recording')
if not os.path.exists(path):
    path = data_path('04_pendulum_hw_backup')
    d, h, w_ph = 1.05, 0.150, 0.072
    print('*** data/04_pendulum_my_recording.csv not found: using the SIMULATED backup (d = 1.05 m).')
    print('*** The homework requires your own recording unless you have no phone (see homework.md).')
rec = load_phyphox(path)
t, w = rec['t'], rec['w']
print('%s: %.1f Hz, %.0f s, swing axis %s' % (os.path.basename(path), rec['fs'], t[-1], rec['axis']))

fig, ax = plt.subplots(figsize=(9, 3.5))
for a, c in rec['columns'].items():
    ax.plot(rec['raw'].iloc[:, 0], rec['raw'][c], lw=0.6, label=c)
ax.set(xlabel='time (s)', ylabel='rotation rate (rad/s)', title='Raw gyroscope log')
ax.grid(True)
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% Problem 1: release, bias, angle
k = find_release(t, w)
th, bias = angle_from_gyro(t, w, k)
print('release at %.2f s, bias %.4f rad/s' % (t[k], bias))
# TODO: time-stamp jitter (std of np.diff of the raw time column); your own still segment, the bias and
#       its standard error, the drift that error causes over the record; release angle (angle at release
#       minus the final equilibrium)

# %% Problem 2: logarithmic decrement
# TODO: 0.1 s moving average of w - bias and its gain at omega_d; ident.log_decrement on the first 20 peaks
#       and on cycles 40-60; ln(peak) vs peak number with np.polyfit(..., cov=True) for the slope error;
#       repeat with 10 and 30 peaks

# %% Problem 3: second-order fit of the angle
# TODO: m = t < t[k] + 40; ident.fit_second_order_step(t[m], th[m], 1.0, t_step=t[k]); overlay and residual;
#       two more windows of your choice

# %% Problem 4: length
# TODO: sqrt(g/d), sqrt(g/L_equivalent(d, h, w_ph)), propagated uncertainty, implied lengths

# %% Problems 5 and 6: poles and settling time
# TODO: poles with error bars, omega_d vs omega_n, 4/(zeta wn) from 2a, 2b and 3 with uncertainties,
#       measured settling from cycles(t, th, k), envelopes on a log scale

# %% Problem 7: period and damping vs amplitude
# TODO: tc, T, A = cycles(t, th, k); least-squares fit T = T0 (1 + c A^2) over your amplitude range with the
#       standard error of c; zeta(A) from 5-cycle log decrements
