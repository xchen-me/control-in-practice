# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 02, demo 1: a cruise-control model from physics and from a logged throttle step (L3, L6-L7)
#
# ![scene](figures/scene.png)
#
# **Overview.** A car at highway speed behaves like a mass pushed by the engine against drag. Linearized
# about the cruise speed it is a first-order lag (time constant over a minute) in series with a short
# engine lag. We first predict the model from the physics, then identify it from a logged throttle step on a
# flat road, and explain the difference.
#
# **Mechanism.** Throttle $u$ (%) sets the engine force $F = F_{max} u/100$ after a lag $\tau_e$. Aerodynamic
# drag $\tfrac12\rho C_d A v^2$ grows with speed, so the car settles where force and road load balance.
#
# **Math.** $m_{eff}\dot v = F - \tfrac12\rho C_d A v^2 - C_{rr} m g - m g\sin\vartheta$. About $v_0$ on a flat
# road, $\delta(\text{drag}) \approx b\,\delta v$ with $b = \rho C_d A v_0$, so
# $G(s) = \dfrac{\delta V}{\delta U} = \dfrac{K}{(\tau s + 1)(\tau_e s + 1)}$, $K = \dfrac{F_{max}/100}{b}$, $\tau = \dfrac{m_{eff}}{b}$.

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ctrlprac import data, ident
from _model import TS, linearize

df, meta = data.load('02_cruise_step')
car = meta['parameters']
lin = linearize(car)
print('%s: m = %.0f kg (m_eff %.0f kg), CdA = %.3f m^2, Crr = %.3f, Fmax = %.0f N, tau_e = %.1f s'
      % (meta['vehicle'], car['m'], car['m_eff'], car['Cd'] * car['A'], car['Crr'], car['Fmax'], car['tau_e']))
print('physics at v0 = %.0f m/s: b = %.1f N s/m, K = %.3f (m/s)/%%, tau = %.1f s, trim throttle = %.1f %%'
      % (car['v0'], lin['b'], lin['K'], lin['tau'], lin['u0']))

# %% [markdown]
# ## Identify from the log (wheel speed and GPS speed)

# %%
t, u = df['time_s'].values, df['throttle_pct'].values
t_step = t[np.argmax(np.diff(u) != 0) + 1]
du = u[-1] - u[0]
fits = {}
for col in ('speed_wheel_mps', 'speed_gps_mps'):
    fits[col] = f = ident.fit_fopdt(t, df[col].values, du, t_step)
    print('%-16s: K = %.3f (m/s)/%%, tau = %.1f s, theta = %.2f s, v0 = %.2f m/s, rmse = %.3f m/s'
          % (col, f['K'], f['tau'], f['theta'], f['y0'], f['rmse']))
fw = fits['speed_wheel_mps']
v1 = fw['y0'] + fw['K'] * du
b_data, m_data = car['Fmax'] / 100 / fw['K'], fw['tau'] * car['Fmax'] / 100 / fw['K']
print('from data: b = Fmax/(100 K) = %.1f N s/m, m_eff = tau b = %.0f kg (spec %.0f kg)'
      % (b_data, m_data, car['m_eff']))

# %% [markdown]
# ## Physics vs data
# The step moved the car from about 27.2 to 30.5 m/s. Drag is quadratic, so its slope over the step is the
# chord $\tfrac12\rho C_dA(v_0 + v_1)$, larger than the tangent at 27 m/s: the car is "stiffer" (smaller $K$,
# shorter $\tau$) than the tangent model predicts.

# %%
v_mid = 0.5 * (fw['y0'] + v1)
lin_mid = linearize(car, v_mid)      # for a quadratic, the chord slope = tangent slope at the mid-speed
print('chord over the step (tangent at %.2f m/s): b = %.1f N s/m -> K = %.3f, tau = %.1f s'
      % (v_mid, lin_mid['b'], lin_mid['K'], lin_mid['tau']))
print('identified vs tangent at v0: K %+.1f %%, tau %+.1f %%; engine lag + delay + one loop sample %.1f s vs fitted'
      ' theta %.2f s' % (100 * (fw['K'] / lin['K'] - 1), 100 * (fw['tau'] / lin['tau'] - 1),
                          car['tau_e'] + car['delay'] + TS, fw['theta']))

tt = t - t_step
model = lambda f, K, tau, th: f['y0'] + K * du * (1 - np.exp(-np.clip(tt - th, 0, None) / tau))
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
ax[0].plot(t, df['speed_gps_mps'], '.', ms=1.5, color='0.75', label='GPS speed')
ax[0].plot(t, df['speed_wheel_mps'], '.', ms=1.5, color='0.4', label='wheel speed')
ax[0].plot(t, model(fw, fw['K'], fw['tau'], fw['theta']), 'C0', lw=2,
           label='FOPDT fit (wheel): K = %.3f, tau = %.0f s' % (fw['K'], fw['tau']))
ax[0].plot(t, model(fw, lin['K'], lin['tau'], car['tau_e']), 'C3--',
           label='physics, tangent at 27 m/s: K = %.3f, tau = %.0f s' % (lin['K'], lin['tau']))
ax[0].set_ylabel('speed [m/s]'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8, loc='lower right')
ax[0].set_title('Throttle step %.0f -> %.0f %% on a flat road' % (u[0], u[-1]))
ax[1].plot(t, u, 'k'); ax[1].set_ylabel('throttle [%]'); ax[1].set_xlabel('time [s]'); ax[1].grid(alpha=.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# The chord model explains most of the gap; the remaining few percent come from wind gusts during the test
# (a slow random force the car integrates). The GPS and wheel-speed fits agree on $K$ and $\tau$. The car applies
# each throttle command one 0.1 s loop sample after it is logged and then waits its 0.1 s delay, so it has 0.8 s of
# engine lag and delay in all. The fitted dead time of 2.4 s is much longer and poorly determined by a test whose
# time constant is over a minute, so a design based on it is conservative.
