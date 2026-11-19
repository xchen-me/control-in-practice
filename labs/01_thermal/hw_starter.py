# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 01 homework starter: heater board B
#
# ![scene](figures/scene_hw.png)
#
# Loads both logs, plots them, and leaves the analysis to you (see homework.md). Runs as shipped.

# %%
import numpy as np
import control as ct
import matplotlib.pyplot as plt
from ctrlprac import data, ident, sim

TS, U_LIM = 1.0, (0.0, 100.0)
step, _ = data.load('01_hw_step')
pctl, _ = data.load('01_hw_pcontrol')
t, u, T = step['time_s'].values, step['heater_pct'].values, step['T_C'].values

fig, ax = plt.subplots(2, 2, figsize=(10, 6), sharex='col')
ax[0, 0].plot(t, T, '.', ms=2); ax[0, 0].set_ylabel('T [degC]'); ax[0, 0].set_title('01_hw_step')
ax[1, 0].plot(t, u); ax[1, 0].set_ylabel('heater [%]'); ax[1, 0].set_xlabel('time [s]')
ax[0, 1].plot(pctl['time_s'], pctl['T_C'], '.', ms=2, label='T')
ax[0, 1].plot(pctl['time_s'], pctl['setpoint_C'], 'k', lw=.8, label='set point')
ax[0, 1].set_title('01_hw_pcontrol (Kp = 4 %/degC)'); ax[0, 1].legend()
ax[1, 1].plot(pctl['time_s'], pctl['heater_pct']); ax[1, 1].set_xlabel('time [s]')
for a in ax.ravel():
    a.grid(alpha=.3)
fig.tight_layout()
plt.show()

# %% Problem 1: identification
# TODO: find the step time and size, then  fit = ident.fit_fopdt(t, T, du, t_step)
# TODO: plot the fit and the residual; graphical check of K and theta + tau

# %% Problem 2: P control
# TODO: measured vs predicted e_ss; Kp for e_ss <= 1 degC; ultimate gain with the exact delay

# %% Problem 3: PI by root locus
# TODO: Gp = fit['model'] * ct.tf(*ct.pade(theta, 1)); locus of (s + 1/tau)/s * Gp; Kc for zeta = 0.7

# %% Problem 4: implementation
# TODO: sim.simulate_loop(fit['model'], sim.PID(Kp, Ki, Ts=TS, u_lim=U_LIM, antiwindup=...), TS, 1500, 55.0,
#       u_lim=U_LIM, noise_std=0.04, quant=0.1, delay=theta, y0=T_ambient)

# %% Problem 5: margins with the exact delay
# TODO: evaluate L(jw) = C(jw) G(jw) exp(-jw theta) on a frequency grid; find wc, PM, w180, GM

# %% Problem 6: engineering judgement
# TODO: 45 -> 55 degC simulations, time to stay within +-0.5 degC, and your recommendation
