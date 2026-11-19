# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 03, demo 1: DC gearmotor, from physics to K/(s(tau s+1)) and back from a logged step
#
# ![scene](figures/scene.png)
#
# **Overview.** A brushed DC gearmotor turns an arm; an encoder counts its angle at 1 kHz. Physics gives a
# third-order transfer function from voltage to angle (L3). The electrical pole is about 40 times faster
# than the mechanical one, so the dominant-pole model K/(s(tau s+1)) is enough for design (L8).
# We then recover K and tau from a logged voltage step (L7), working only with what the hardware gives:
# the voltage command and integer encoder counts.

# %%
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.modules.pop('_model', None)          # run_all.py runs several labs in one process
import numpy as np
import matplotlib.pyplot as plt
import control as ct
from _model import MOTOR_A, full_tf, reduced_params, reduced_tf, J_total, counts_per_rev, \
    identify_from_log, speed_from_counts

p = MOTOR_A
print('datasheet-class values (representative):', {k: p[k] for k in ('R', 'L', 'kt', 'Jm', 'JL', 'N', 'b')})

# %% [markdown]
# ## Mechanism and model (L3)
# Armature: $L\,\dot i = V - R i - k_e\omega_m$. Rotor + reflected load: $J\dot\omega_m = k_t i - b\omega_m$,
# with $J = J_m + J_L/N^2$. Output angle $\theta = \theta_m/N$. With $k_e = k_t$ (SI units):
#
# $$G(s)=\frac{\Theta(s)}{V(s)}=\frac{k_t/N}{s\,[(Ls+R)(Js+b)+k_t^2]}.$$
#
# Setting $L\to 0$: $G(s)\approx \dfrac{K}{s(\tau s+1)}$, $K=\dfrac{k_t}{N(Rb+k_t^2)}$,
# $\tau=\dfrac{RJ}{Rb+k_t^2}$.

# %%
G3 = full_tf(p)
K0, tau0 = reduced_params(p)
poles = np.sort(np.roots(G3.den[0][0]).real)
print('J = Jm + JL/N^2 = %.3g kg m^2 (load share %.0f %%)' % (J_total(p), 100 * p['JL'] / p['N'] ** 2 / J_total(p)))
print('poles of the full model:', np.round(poles, 1), '1/s')
print('electrical pole / mechanical pole = %.0f' % (poles[0] / poles[1]))
print('reduced model (L -> 0): K = %.4f rad/s/V, tau = %.4f s;  1/|mechanical pole| = %.4f s'
      % (K0, tau0, -1 / poles[1]))

# %% Dominant poles: speed step responses of the full and reduced models (L8)
t = np.linspace(0, 0.2, 2001)
s = ct.tf('s')
_, w3 = ct.step_response(G3 * s, t)                       # speed = s * angle
_, w2 = ct.step_response(reduced_tf(K0, tau0) * s, t)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
ax[0].plot(t * 1e3, w3, label='third order (with L)')
ax[0].plot(t * 1e3, w2, '--', label='K/(tau s+1) (L neglected)')
ax[0].set(xlabel='time (ms)', ylabel='output speed per volt (rad/s/V)', title='Unit voltage step')
ax[1].plot(t * 1e3, w3, t * 1e3, w2, '--')
ax[1].set(xlim=(0, 10), ylim=(0, 0.4), xlabel='time (ms)', ylabel='(rad/s/V)', title='First 10 ms')
for a in ax:
    a.grid(True)
ax[0].legend()
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Identification from the logged step (L7)
# The log holds the voltage command and the encoder count. Speed is a difference of counts; one count is
# $2\pi/1920$ rad, so a one-sample difference has a resolution of 3.3 rad/s (useless), while a centered
# 10-sample difference has 0.33 rad/s and no lag (fine offline). `fit_fopdt` then gives K, tau and a small
# dead time (PWM update delay plus the neglected electrical lag).

# %%
fit, (t, V, c, w10, k_step, meta) = identify_from_log('03_motor_step', window=10)
cpr = meta['counts_per_rev']
w1 = speed_from_counts(t, c, cpr, window=1, causal=True)
true = meta['true_reduced_model']
tau_dom = -1 / meta['true_reduced_model']['poles_of_full_model'][1]
print('identified: K = %.4f rad/s/V (true %.4f, error %+.2f %%)' % (fit['K'], true['K_rad_per_s_per_V'],
      100 * (fit['K'] / true['K_rad_per_s_per_V'] - 1)))
print('identified: tau = %.4f s (dominant pole %.4f s, error %+.1f %%; L-neglected formula %.4f s)'
      % (fit['tau'], tau_dom, 100 * (fit['tau'] / tau_dom - 1), true['tau_s']))
print('dead time theta = %.1f ms, fit rmse = %.3f rad/s' % (fit['theta'] * 1e3, fit['rmse']))
print('operating point before the step: %.2f rad/s at %.1f V (K*V = %.2f: Coulomb friction offset)'
      % (fit['y0'], V[0], fit['K'] * V[0]))

tt = t - t[k_step]
yfit = fit['y0'] + fit['K'] * (V[k_step] - V[0]) * (1 - np.exp(-np.clip(tt - fit['theta'], 0, None) / fit['tau']))
fig, ax = plt.subplots(2, 1, figsize=(8, 5.5), sharex=True)
ax[0].plot(t, V, 'k')
ax[0].set(ylabel='V_cmd (V)', title='Logged voltage step, Motor A')
ax[1].plot(t, w1, '.', ms=2, color='0.7', label='1-sample count difference')
ax[1].plot(t, w10, label='centered 10-sample difference')
ax[1].plot(t, yfit, 'r--', label='fit: K=%.3f, tau=%.1f ms' % (fit['K'], fit['tau'] * 1e3))
ax[1].set(xlabel='time (s)', ylabel='output speed (rad/s)', ylim=(0, 20))
for a in ax:
    a.grid(True)
ax[1].legend(loc='lower right')
fig.tight_layout()
plt.show()

# %% [markdown]
# **Takeaways.** The identified K and tau match the dominant pole of the physics model to well under 1 %.
# Note the offset: the speed at 3 V is below K x 3 V because Coulomb friction eats about 0.25 V; this is
# why the test steps between two nonzero voltages and fits the *increment*. The design model for the
# position loop is $G(s)=K/(s(\tau s+1))$ with these numbers.
