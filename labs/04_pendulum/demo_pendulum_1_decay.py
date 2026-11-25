# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 04, demo 1: a phone pendulum, damping and natural frequency from a free decay
#
# ![scene](figures/scene.png)
#
# **Overview.** Hang a phone on a string, pull it aside, hold it still, release it, and record the
# gyroscope with phyphox. The gyro measures the swing rate $\dot\theta$ directly. Near the bottom the
# pendulum is a second-order system $\ddot\theta + 2\zeta\omega_n\dot\theta + \omega_n^2\theta = 0$ (L6),
# so the free decay gives $\zeta$ and $\omega_n$ two ways: the logarithmic decrement and a second-order
# fit (L7, L8). Physics predicts $\omega_n$ from the length; the poles $-\zeta\omega_n \pm j\omega_d$
# (L16) predict the settling time, which we check against the record.
#
# The sample file is **simulated** (see its JSON). Replace it with your own recording
# (README: how to record).

# %%
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.modules.pop('_model', None)          # run_all.py runs several labs in one process
import json
import numpy as np
import matplotlib.pyplot as plt
from ctrlprac import ident
from _model import load_phyphox, find_release, angle_from_gyro, L_equivalent, data_path, g

name = '04_pendulum_sample'              # your file: '04_pendulum_my_recording'
d, h, w_ph = 0.75, 0.150, 0.072          # your measured d, phone height and width (m)
rec = load_phyphox(data_path(name))
jpath = data_path(name).replace('.csv', '.json')
if os.path.exists(jpath):                # the simulated sample carries its setup in the JSON
    meta = json.load(open(jpath))
    print('%s: source = %s' % (name, meta['source']))
    d, h, w_ph = meta['parameters']['d_m'], meta['parameters']['phone_h_m'], meta['parameters']['phone_w_m']
print('measured setup: pivot to phone centre d = %.3f m, phone %.3f x %.3f m' % (d, h, w_ph))
print('sampling %.1f Hz, axis energies (rad/s)^2: %s -> swing axis "%s"'
      % (rec['fs'], {k: round(float(v), 4) for k, v in rec['energy'].items()}, rec['axis']))
t, w = rec['t'], rec['w']

# %% Raw log: all three axes
fig, ax = plt.subplots(figsize=(9, 3.5))
for a, c in rec['columns'].items():
    ax.plot(rec['raw'].iloc[:, 0], rec['raw'][c], lw=0.6, label=c)
ax.set(xlabel='time (s)', ylabel='rotation rate (rad/s)',
       title='phyphox gyroscope export (%s)' % ('simulated sample' if os.path.exists(jpath) else name))
ax.grid(True)
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Release, bias and angle
# The phone is held still for a few seconds before the release; the gyro mean there is its bias. Removing
# it and integrating gives the angle (zero at the start of the record). A residual bias would make the
# angle drift; the pendulum's equilibrium cannot move, so a fitted straight-line drift is removed too.

# %%
k = find_release(t, w)
th, bias = angle_from_gyro(t, w, k)
w0 = w - bias
print('release at t = %.2f s, gyro bias %.4f rad/s, release angle %.1f deg'
      % (t[k], bias, np.degrees(abs(th[k] - np.mean(th[-len(th) // 5:])))))

# %% [markdown]
# ## Logarithmic decrement (L7)
# For successive peaks $x_i$: $\delta = \ln(x_i/x_{i+1}) = 2\pi\zeta/\sqrt{1-\zeta^2}$,
# $\omega_d = 2\pi/T_d$, $\omega_n = \omega_d/\sqrt{1-\zeta^2}$. Use a fixed number of peaks after the
# release (the amplitude range matters, see demo 2).

# %%
ld = ident.log_decrement(t[k:], w0[k:], n_peaks=30)
print('log decrement (first 30 peaks): zeta = %.4f, omega_n = %.4f rad/s, T_d = %.4f s'
      % (ld['zeta'], ld['wn'], 2 * np.pi / ld['wd']))

# %% [markdown]
# ## Second-order fit (L6, L8)
# Released from rest at $\theta_0$, the angle is $\theta_0$ minus $\theta_0$ times the unit step response
# of $\omega_n^2/(s^2+2\zeta\omega_n s+\omega_n^2)$. So `fit_second_order_step` on the angle, with the
# release as the "step", returns $\zeta$, $\omega_n$ and $K=-\theta_0$. A 40 s window keeps many cycles.

# %%
m = t < t[k] + 40
fit = ident.fit_second_order_step(t[m], th[m], 1.0, t_step=t[k])
print('second-order fit (40 s): zeta = %.4f, omega_n = %.4f rad/s, K = %.3f rad (theta0 = %.1f deg), rmse = %.4f rad'
      % (fit['zeta'], fit['wn'], fit['K'], np.degrees(-fit['K']), fit['rmse']))

fig, ax = plt.subplots(2, 1, figsize=(9, 5.5))
ax[0].plot(t, w0, lw=0.6, label='gyro, bias removed')
ax[0].plot(ld['peaks_t'], ld['peaks_y'], 'r.', label='peaks for the log decrement')
ax[0].set(xlabel='time (s)', ylabel='swing rate (rad/s)', title='Free decay')
tt = np.linspace(0, 40, 4001)
from control import step_response
_, sr = step_response(fit['model'], tt)
ax[1].plot(t[m] - t[k], th[m], lw=1.2, label='angle (integrated gyro)')
ax[1].plot(tt, sr, 'r--', lw=0.9, label='second-order fit')
ax[1].set(xlim=(-1, 12), xlabel='time after release (s)', ylabel='angle (rad)')
for a in ax:
    a.grid(True)
    a.legend(fontsize=8, loc='upper right')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Compare with the length
# Simple pendulum: $\omega_n = \sqrt{g/d}$. The phone is not a point mass: as an $h\times w$ plate swinging in
# its own plane, $I_{pivot} = m d^2 + m(h^2+w^2)/12$ and $\omega_n=\sqrt{g/L_{eq}}$, $L_{eq}=d+(h^2+w^2)/(12d)$.
# Conversely, the measured $\omega_n$ gives a length $L = g/\omega_n^2$.

# %%
Leq = L_equivalent(d, h, w_ph)
wn_simple, wn_phys = np.sqrt(g / d), np.sqrt(g / Leq)
print('omega_n: simple pendulum %.4f, physical pendulum %.4f rad/s (L_eq = %.4f m)' % (wn_simple, wn_phys, Leq))
for lab, wn in (('log decrement', ld['wn']), ('2nd-order fit', fit['wn'])):
    print('  %-14s %.4f rad/s: %+.2f %% vs physical; implied length g/wn^2 = %.3f m'
          % (lab, wn, 100 * (wn / wn_phys - 1), g / wn ** 2))
print('a 1 cm error in d changes omega_n by %.2f %%' % (50 * 0.01 / d))

# %% [markdown]
# ## Poles and settling time (L16)
# $s = -\zeta\omega_n \pm j\omega_n\sqrt{1-\zeta^2}$; $T_s \approx 4/(\zeta\omega_n)$ for the 2 % band.
# Measured: first time the cycle amplitude stays below 2 % of the release amplitude.

# %%
from _model import cycles
tc, T, A = cycles(t, th, k)
A0 = abs(th[k] - np.mean(th[-len(th) // 5:]))     # release amplitude
below = np.nonzero(A > 0.02 * A0)[0]
Ts_meas = tc[below[-1] + 1] - t[k] if len(below) and below[-1] + 1 < len(tc) else np.nan
fig, ax = plt.subplots(figsize=(5, 4.2))
for lab, z, wn, mk in (('log decrement', ld['zeta'], ld['wn'], 'o'), ('2nd-order fit', fit['zeta'], fit['wn'], 's')):
    p = complex(-z * wn, wn * np.sqrt(1 - z ** 2))
    print('%-14s poles %.4f +/- %.4fj 1/s -> Ts = 4/(zeta wn) = %.0f s' % (lab, p.real, p.imag, 4 / (z * wn)))
    ax.plot([p.real] * 2, [p.imag, -p.imag], mk, mfc='none', ms=9, label=lab)
print('measured 2 %% settling time: %s'
      % ('%.0f s' % Ts_meas if np.isfinite(Ts_meas) else 'not reached in the record'))
ax.axvline(0, color='k', lw=0.8)
ax.set(xlim=(-0.3, 0.05), ylim=(-4.5, 4.5), xlabel='Re(s) (1/s)', ylabel='Im(s) (rad/s)', title='Pendulum poles')
ax.grid(True)
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# **Takeaways.** Both methods give $\omega_n$ within about 1 % of $\sqrt{g/L_{eq}}$, slightly low because
# the release amplitude was large (demo 2). The damping is tiny and the poles hug the imaginary axis. The
# settling-time estimate from early-cycle damping is far shorter than the measured one: the damping is
# mostly air drag, which shrinks with amplitude, so one $\zeta$ does not describe the whole decay.
