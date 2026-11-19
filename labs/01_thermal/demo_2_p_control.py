# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 01, demo 2: proportional control leaves a steady-state error (L12)
#
# ![scene](figures/scene.png)
#
# **Overview.** With $u = K_p (r - T)$ the heater power is zero when the error is zero, but holding the
# board above ambient needs power. So the loop settles where the error is just large enough to
# supply that power. Raising $K_p$ shrinks the error but the delay makes the loop oscillate.
#
# **Math.** Around ambient, $T - T_a = G(s) u$ with $G(0) = K$. The loop $L = K_p G$ is type 0, so for
# a set-point step of size $r - T_a$:
# $e_{ss} = \dfrac{r - T_a}{1 + K K_p}$ (position error constant $K_{pos} = K K_p$).

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ctrlprac import data, ident, sim
from _model import TS, QUANT, NOISE, U_LIM, true_plant

df, meta = data.load('01_thermal_step')
t, u, T = df['time_s'].values, df['heater_pct'].values, df['T_C'].values
fit = ident.fit_fopdt(t, T, u[-1] - u[0], 60.0)
p = meta['parameters']
amb, r = p['amb'], 50.0
print('identified K = %.3f degC/%%; set point %.0f degC, ambient %.0f degC' % (fit['K'], r, amb))

# %% [markdown]
# ## Predicted vs simulated steady-state error (simulation uses the true plant with delay)

# %%
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for i, Kp in enumerate([2.0, 5.0, 10.0, 20.0]):
    ess = (r - amb) / (1 + fit['K'] * Kp)
    res = sim.simulate_loop(true_plant(p), sim.PID(Kp, Ts=TS, u_lim=U_LIM), TS, 1500, r, u_lim=U_LIM,
                            noise_std=NOISE, quant=QUANT, delay=p['delay'], y0=amb, seed=10 + i)
    ess_sim = r - res['y'][-200:].mean()
    print('Kp = %4.1f %%/degC: K*Kp = %5.2f, predicted e_ss = %.2f degC, simulated e_ss = %.2f degC'
          % (Kp, fit['K'] * Kp, ess, ess_sim))
    ax[0].plot(res['t'], res['y_meas'], 'C%d' % i, lw=1, label='Kp = %g (e_ss pred. %.1f degC)' % (Kp, ess))
    ax[0].axhline(r - ess, color='C%d' % i, lw=.6, ls=':')
    ax[1].plot(res['t'], res['u'], 'C%d' % i, lw=1)
ax[0].axhline(r, color='k', lw=1, label='set point')
ax[0].set_ylabel('T [degC]'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8, loc='lower right')
ax[0].set_title('P control of board A: higher Kp, smaller error, more oscillation')
ax[1].set_ylabel('heater [%]'); ax[1].set_xlabel('time [s]'); ax[1].grid(alpha=.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# The error never reaches zero, and $K_p = 20$ already rings because the phase lag of the delay
# grows with frequency (demo 4). Removing the error needs an integrator: PI control (demo 3).
