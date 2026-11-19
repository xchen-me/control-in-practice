# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 01, demo 1: a first-order-plus-dead-time model from a heater step test (L6-L7)
#
# ![scene](figures/scene.png)
#
# **Overview.** We switch the heater of board A from 0 % to 50 % and log the temperature every second.
# The response looks like a first-order rise that starts a little late. We fit
# $G(s) = K e^{-\theta s}/(\tau s + 1)$ to the log and check it against the "true" two-lump plant
# that generated the (simulated) data.
#
# **Mechanism.** The transistor heats an aluminium fin (large thermal mass, slow, $\tau_1$); the
# thermistor sits on the fin and lags it a little ($\tau_2$); there is a few-second transport delay.
# Two slow lumps in series look like one lag plus an *apparent* dead time.
#
# **Math.** Step of size $\Delta u$ at $t_s$:
# $T(t) = T_0 + K \Delta u \,(1 - e^{-(t - t_s - \theta)/\tau})$ for $t > t_s + \theta$.
# $K$ is the steady rise per % heater, $\tau$ is where 63 % of the rise is reached (after $\theta$).

# %%
import os, sys
import numpy as np
import control as ct
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ctrlprac import data, ident
from _model import true_plant

df, meta = data.load('01_thermal_step')
t, u, T = df['time_s'].values, df['heater_pct'].values, df['T_C'].values
t_step = t[np.argmax(np.diff(u) != 0) + 1]
du = u[-1] - u[0]
print('step: %.0f %% -> %.0f %% at t = %.0f s' % (u[0], u[-1], t_step))

# %% [markdown]
# ## Fit the FOPDT model (least squares on the whole log)

# %%
fit = ident.fit_fopdt(t, T, du, t_step)
print('identified: K = %.3f degC/%%, tau = %.1f s, theta = %.1f s, T0 = %.2f degC, rmse = %.3f degC'
      % (fit['K'], fit['tau'], fit['theta'], fit['y0'], fit['rmse']))

# graphical check: 63 % point
dT = T[-30:].mean() - fit['y0']
t63 = t[np.argmax(T >= fit['y0'] + 0.632 * dT)] - t_step
print('graphical: dT = %.2f degC -> K = %.3f; t63 = %.0f s (= theta + tau = %.0f s)'
      % (dT, dT / du, t63, fit['theta'] + fit['tau']))

# %% [markdown]
# ## Compare with the true plant (known only because the data are simulated)

# %%
p = meta['parameters']
G_true = true_plant(p)
print('true: K = %.3f, tau1 = %.0f s, tau2 = %.0f s, transport delay = %.0f s'
      % (p['K'], p['tau1'], p['tau2'], p['delay']))
print('K error = %.1f %%; fit tau = %.0f s vs tau1 = %.0f s (%+.1f %%)'
      % (100 * (fit['K'] / p['K'] - 1), fit['tau'], p['tau1'], 100 * (fit['tau'] / p['tau1'] - 1)))
print('apparent delay %.1f s = %.0f s transport + %.1f s that stands in for the %.0f s sensor lag'
      % (fit['theta'], p['delay'], fit['theta'] - p['delay'], p['tau2']))

tt = t - t_step
y_fit = fit['y0'] + fit['K'] * du * (1 - np.exp(-np.clip(tt - fit['theta'], 0, None) / fit['tau']))
tg = np.arange(len(t)) * 1.0
_, y_true = ct.step_response(G_true, tg)
y_true = p['amb'] + du * np.interp(tt - p['delay'], tg, y_true, left=0)

fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
ax[0].plot(t, T, '.', ms=2, color='0.5', label='logged T (0.1 degC steps)')
ax[0].plot(t, y_fit, 'C0', lw=2, label='FOPDT fit: K=%.3f, tau=%.0f s, theta=%.0f s'
           % (fit['K'], fit['tau'], fit['theta']))
ax[0].plot(t, y_true, 'C3--', lw=1.2, label='true two-lump plant')
ax[0].set_ylabel('temperature [degC]'); ax[0].grid(alpha=.3); ax[0].legend(loc='lower right')
ax[0].set_title('Board A: heater 0 -> 50 % step')
ax[1].plot(t, T - y_fit, 'C0', lw=.8)
ax[1].set_ylabel('residual [degC]'); ax[1].set_xlabel('time [s]'); ax[1].grid(alpha=.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# Zoom on the start: the two-lump plant starts with zero slope (an S-shape), which the FOPDT
# fit explains with a longer dead time than the physical transport delay.

# %%
fig, ax = plt.subplots(figsize=(7, 4))
m = (t > t_step - 20) & (t < t_step + 120)
ax.plot(t[m] - t_step, T[m], '.', color='0.5', label='logged')
ax.plot(t[m] - t_step, y_fit[m], 'C0', label='FOPDT fit')
ax.plot(t[m] - t_step, y_true[m], 'C3--', label='true')
ax.axvline(p['delay'], color='k', lw=.8, ls=':'); ax.text(p['delay'] + 1, fit['y0'] + 3, 'transport delay')
ax.axvline(fit['theta'], color='C0', lw=.8, ls=':'); ax.text(fit['theta'] + 1, fit['y0'] + 5, 'fitted theta')
ax.set_xlabel('time after step [s]'); ax.set_ylabel('temperature [degC]'); ax.grid(alpha=.3); ax.legend()
fig.tight_layout()
plt.show()
