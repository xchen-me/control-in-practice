# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 04, demo 2: where the linear second-order model stops: large angles and air drag
#
# ![scene](figures/scene.png)
#
# **Overview.** The model $\ddot\theta + 2\zeta\omega_n\dot\theta + \omega_n^2\theta = 0$ assumes
# $\sin\theta\approx\theta$ and damping proportional to speed. A real swing breaks both: at large amplitude
# the period grows, $T \approx T_0(1+\theta_a^2/16)$, and air drag grows with speed squared, so the
# damping ratio depends on the amplitude. Cycle-by-cycle analysis of the same free decay shows both, and
# explains why demo 1's settling-time estimate missed (L6-L8, L16).

# %%
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.modules.pop('_model', None)          # run_all.py runs several labs in one process
import json
import numpy as np
import matplotlib.pyplot as plt
from _model import load_phyphox, find_release, angle_from_gyro, cycles, period_large_angle, L_equivalent, \
    data_path, g

name = '04_pendulum_sample'              # your file: '04_pendulum_my_recording'
par = dict(d_m=0.75, phone_h_m=0.150, phone_w_m=0.072)    # your measured setup (m)
rec = load_phyphox(data_path(name))
jpath = data_path(name).replace('.csv', '.json')
if os.path.exists(jpath):                # simulated sample: setup and true damping from the JSON
    par = json.load(open(jpath))['parameters']
t, w = rec['t'], rec['w']
k = find_release(t, w)
th, bias = angle_from_gyro(t, w, k)
tc, T, A = cycles(t, th, k)
if len(T) < 8:
    raise ValueError('only %d clean swing cycles after the release at %.1f s; check that the phone swings freely '
                     'and that the record covers the decay' % (len(T), t[k]))
A_rel = abs(th[k] - np.mean(th[-len(th) // 5:]))     # release amplitude, as in demo 1
print('%d cycles with amplitude from %.1f to %.1f deg' % (len(T), np.degrees(A[0]), np.degrees(A[-1])))

# %% [markdown]
# ## Period versus amplitude
# Fit $T = T_0(1 + c\,\theta_a^2)$ by least squares (linear in $T_0$ and $T_0c$). Theory: $c = 1/16$.
# The intercept $T_0$ is the small-angle period, which gives the length $L = g(T_0/2\pi)^2$.

# %%
M = np.column_stack([np.ones_like(A), A ** 2])
b0, b1 = np.linalg.lstsq(M, T, rcond=None)[0]
T0, c = b0, b1 / b0
Leq = L_equivalent(par['d_m'], par['phone_h_m'], par['phone_w_m'])
T0_th = 2 * np.pi * np.sqrt(Leq / g)
print('fit: T0 = %.4f s (theory %.4f s, %+.2f %%), c = %.4f (theory 1/16 = 0.0625)'
      % (T0, T0_th, 100 * (T0 / T0_th - 1), c))
print('length from T0: %.3f m (L_eq %.3f m); first cycle at %.0f deg is %.1f %% slower than T0'
      % (g * (T0 / 2 / np.pi) ** 2, Leq, np.degrees(A[0]), 100 * (T[0] / T0 - 1)))

# %% [markdown]
# ## Damping versus amplitude
# Per cycle, $\zeta_k \approx \ln(A_k/A_{k+1})/(2\pi)$ (averaged over 5 cycles to tame noise). For damping
# $2\zeta_l\omega_n\dot\theta + c_q|\dot\theta|\dot\theta$, energy balance over one cycle gives the
# equivalent $\zeta(A) = \zeta_l + \tfrac{4}{3\pi}c_qA$: a straight line in the amplitude.

# %%
n = 5
zk = np.log(A[:-n] / A[n:]) / (2 * np.pi * n)
Am = np.sqrt(A[:-n] * A[n:])
p1 = np.polyfit(Am, zk, 1)
zeta_l, c_q = p1[1], p1[0] * 3 * np.pi / 4
print('zeta(A) = %.4f + %.4f A  ->  zeta_lin = %.4f, c_quad = %.3f 1/rad' % (p1[1], p1[0], zeta_l, c_q))
if 'zeta_lin' in par:
    print('   simulation truth: zeta_lin = %.4f, c_quad = %.3f 1/rad' % (par['zeta_lin'], par['c_quad_per_rad']))

# envelope prediction: dA/dt = -zeta(A) wn A, from the release amplitude to 2 % of it (the threshold of demo 1)
wn = 2 * np.pi / T0
T_env = tc[-1] - t[k] + 60
if zeta_l > 0:                           # long enough for the linear damping alone to reach 2 %
    T_env = max(T_env, 1.2 * np.log(50) / (zeta_l * wn))
te = np.linspace(0, T_env, 4000)
Ae = np.zeros_like(te)
Ae[0] = A_rel
for i in range(1, len(te)):
    Ae[i] = max(Ae[i - 1] - (te[i] - te[i - 1]) * (zeta_l + p1[0] * Ae[i - 1]) * wn * Ae[i - 1], 0.0)
hit = np.nonzero(Ae < 0.02 * A_rel)[0]
Ts_nl = te[hit[0]] if len(hit) else np.nan
z_early = np.mean(zk[:5])
print('2 %% settling: linear model with early zeta %.4f -> %.0f s; amplitude-dependent model -> %s'
      % (z_early, 4 / (z_early * wn), '%.0f s' % Ts_nl if np.isfinite(Ts_nl) else 'beyond %.0f s' % te[-1]))

# %%
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
aa = np.linspace(0, A[0] * 1.05, 100)
ax[0].plot(np.degrees(A), T, '.', ms=4, label='measured cycles')
ax[0].plot(np.degrees(aa), period_large_angle(T0_th, aa), 'k', lw=1, label=r'$T_0(1+\theta_a^2/16+\dots)$, $L_{eq}$')
ax[0].plot(np.degrees(aa), T0 * (1 + c * aa ** 2), 'r--', lw=1, label='fit')
ax[0].set(xlabel='amplitude (deg)', ylabel='period (s)', title='Period grows with amplitude')
ax[1].plot(np.degrees(Am), zk, '.', ms=4, label='5-cycle log decrement')
ax[1].plot(np.degrees(aa), np.polyval(p1, aa), 'r--', lw=1, label=r'$\zeta_l + \frac{4}{3\pi}c_q A$')
ax[1].set(xlabel='amplitude (deg)', ylabel='damping ratio', title='Damping grows with amplitude')
ax[2].semilogy(tc - t[k], np.degrees(A), '.', ms=4, label='measured amplitude')
ax[2].semilogy(te, np.degrees(A_rel * np.exp(-z_early * wn * te)), 'k:', label='linear, early zeta')
ax[2].semilogy(te[Ae > 0], np.degrees(Ae[Ae > 0]), 'r--', label='amplitude-dependent zeta')
ax[2].axhline(np.degrees(0.02 * A_rel), color='0.5', lw=0.8, label='2 % of release')
ax[2].set(ylim=(0.05, max(60, 1.5 * np.degrees(A_rel))), xlabel='time after release (s)', ylabel='amplitude (deg)',
          title='Envelope')
for a in ax:
    a.grid(True)
    a.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# **Takeaways.** The period-amplitude law is visible with a phone; extrapolating to zero amplitude recovers
# the small-angle period, which is what $\sqrt{g/L}$ predicts. The damping is not one number: a single
# $\zeta$ from the first cycles predicts settling far too early, while $\zeta(A)$ reproduces the slow tail.
# The linear model is still the right local description: around any amplitude it gives the correct
# frequency and decay rate, which is how linearization is used in control design.
