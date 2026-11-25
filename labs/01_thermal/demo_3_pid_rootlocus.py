# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 01, demo 3: PI/PID by root locus, saturation and anti-windup (L13-L19)
#
# ![scene](figures/scene.png)
#
# **Overview.** The integrator of a PI controller removes the steady-state error of demo 2. We place the PI
# zero, pick the gain on the root locus for a damping target, and then run the controller the way the
# board runs it: sampled at 1 s, heater clipped to 0-100 %. A large set-point step saturates the heater;
# without anti-windup the integrator keeps charging and the temperature overshoots.
#
# **Mechanism.** The root locus needs a rational $G(s)$, so the dead time is replaced by a first-order
# Pade approximation $e^{-\theta s} \approx (1 - \theta s/2)/(1 + \theta s/2)$: a stable pole at
# $-2/\theta$ and a right-half-plane zero at $+2/\theta$ that bends the locus to the right.
#
# **Math.** $C(s) = K_c (s + z)/s = K_p + K_i/s$ with $K_p = K_c$, $K_i = K_c z$. Choosing $z = 1/\tau$
# cancels the slow plant pole, leaving $L(s) = \frac{K_c K}{\tau}\,\frac{1}{s}\,\frac{1 - \theta s/2}{1 + \theta s/2}$.

# %%
import os, sys
import numpy as np
import control as ct
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ctrlprac import data, ident, sim
from _model import TS, QUANT, NOISE, U_LIM, true_plant, fopdt_pade, gain_for_zeta, step_metrics

df, meta = data.load('01_thermal_step')
fit = ident.fit_fopdt(df['time_s'], df['T_C'], 50.0, 60.0)
p = meta['parameters']
amb = p['amb']
Gp = fopdt_pade(fit)
print('model: K = %.3f degC/%%, tau = %.1f s, theta = %.1f s; Pade(1) pole/zero at -/+ %.4f rad/s'
      % (fit['K'], fit['tau'], fit['theta'], 2 / fit['theta']))

# %% [markdown]
# ## PI design: zero at $1/\tau$, gain for $\zeta = 0.7$ on the root locus

# %%
ZETA = 0.7
z = 1 / fit['tau']
L1 = ct.tf([1, z], [1, 0]) * Gp            # PI with Kc = 1
Kc, cl = gain_for_zeta(L1, ZETA)
Kp, Ki = Kc, Kc * z
pd = cl[np.argmax(np.where(np.abs(cl.imag) > 1e-9, cl.real, -np.inf))]
wn = abs(pd)
print('PI: z = %.5f rad/s, Kc = %.3f -> Kp = %.3f %%/degC, Ki = %.5f %%/(degC s)' % (z, Kc, Kp, Ki))
print('closed-loop poles (Pade model):', np.round(cl, 4))
print('dominant pair: zeta = %.2f, wn = %.4f rad/s -> 2%% settling ~ 4/(zeta wn) = %.0f s'
      % (ZETA, wn, 4 / (ZETA * wn)))

num, den = np.ravel(L1.num[0][0]), np.ravel(L1.den[0][0])
ks = np.concatenate([[0], np.logspace(-3, 2, 1500)])
branches = np.array([np.roots(np.polyadd(den, k * num)) for k in ks])
fig, ax = plt.subplots(figsize=(7, 5.5))
ax.plot(branches.real.ravel(), branches.imag.ravel(), '.', color='C0', ms=1.5, label='root locus, Kc > 0')
ax.plot(np.roots(den).real, np.roots(den).imag, 'kx', ms=9, label='open-loop poles')
ax.plot(np.roots(num).real, np.roots(num).imag, 'ko', mfc='none', ms=9, label='open-loop zeros')
ax.plot(cl.real, cl.imag, 'C3s', ms=7, label='design, Kc = %.2f' % Kc)
r_ = np.linspace(0, 0.2, 2)
ax.plot(-r_ * ZETA, r_ * np.sqrt(1 - ZETA ** 2), 'k:', lw=.8)
ax.plot(-r_ * ZETA, -r_ * np.sqrt(1 - ZETA ** 2), 'k:', lw=.8, label='zeta = %.1f' % ZETA)
ax.set_xlim(-0.16, 0.14); ax.set_ylim(-0.13, 0.13)
ax.axhline(0, color='k', lw=.5); ax.axvline(0, color='k', lw=.5)
ax.set_xlabel('Re(s) [rad/s]'); ax.set_ylabel('Im(s) [rad/s]'); ax.grid(alpha=.3); ax.legend(fontsize=8)
ax.set_title('Root locus of PI x FOPDT (Pade 1): zero at -1/tau cancels the plant pole')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Implementation on the board: 21 -> 50 degC, heater clipped to 0-100 %, with and without anti-windup
# The simulation uses the **true** plant (two lumps + 4 s delay), not the model the controller was designed on.

# %%
r, T_END = 50.0, 1500.0
runs = {}
for aw in (True, False):
    ctrl = sim.PID(Kp, Ki, Ts=TS, u_lim=U_LIM, antiwindup=aw)
    runs[aw] = sim.simulate_loop(true_plant(p), ctrl, TS, T_END, r, u_lim=U_LIM, noise_std=NOISE,
                                 quant=QUANT, delay=p['delay'], y0=amb, seed=31)
    OS, ts = step_metrics(runs[aw]['t'], runs[aw]['y'], r, amb)
    print('anti-windup %-5s: overshoot = %5.1f %%, 2%% settling = %4.0f s, heater at 100 %% for %3.0f s'
          % (aw, OS, ts, TS * np.sum(runs[aw]['u'][:-1] >= 99.99)))
# linear loop without saturation, for reference
Lin = ct.feedback(ct.tf([Kp, Ki], [1, 0]) * true_plant(p) * ct.tf(*ct.pade(p['delay'], 3)), 1)
_, y_lin = ct.step_response(Lin, runs[True]['t'])
OS, ts = step_metrics(runs[True]['t'], amb + (r - amb) * y_lin, r, amb)
print('linear loop (no saturation): overshoot = %.1f %%, 2%% settling = %.0f s; the clipped heater'
      ' (and the imperfect pole-zero cancellation) make the real loop slower' % (OS, ts))

fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for aw, c in ((True, 'C0'), (False, 'C3')):
    lab = 'clamping anti-windup' if aw else 'no anti-windup (windup)'
    ax[0].plot(runs[aw]['t'], runs[aw]['y_meas'], c, lw=1, label=lab)
    ax[1].plot(runs[aw]['t'], runs[aw]['u'], c, lw=1)
ax[0].plot(runs[True]['t'], amb + (r - amb) * y_lin, 'k--', lw=.8, label='linear, no saturation')
ax[0].axhline(r, color='k', lw=.6)
ax[0].set_ylabel('T [degC]'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8, loc='lower right')
ax[0].set_title('PI (Kp = %.2f, Ki = %.4f) on the true board, heater saturates at 100 %%' % (Kp, Ki))
ax[1].set_ylabel('heater [%]'); ax[1].set_xlabel('time [s]'); ax[1].grid(alpha=.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Model-mismatch check: identified model vs true plant, with a cooling-fan disturbance
# At t = 900 s a fan starts; its extra heat loss is equivalent to -12 % heater power.

# %%
fan = lambda t: -12.0 if t >= 900 else 0.0
cases = [('identified FOPDT model', fit['model'], fit['theta']), ('true two-lump plant', true_plant(p), p['delay'])]
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for i, (name, plant, dly) in enumerate(cases):
    res = sim.simulate_loop(plant, sim.PID(Kp, Ki, Ts=TS, u_lim=U_LIM), TS, 1800, r, d=fan, u_lim=U_LIM,
                            noise_std=NOISE, quant=QUANT, delay=dly, y0=amb, seed=32)
    OS, ts = step_metrics(res['t'][res['t'] < 900], res['y'][res['t'] < 900], r, amb)
    m = res['t'] >= 900
    dip = r - res['y'][m].min()
    print('%-24s: overshoot %.1f %%, settling %.0f s; fan: max dip %.2f degC, e_ss %.3f degC'
          % (name, OS, ts, dip, r - res['y'][-100:].mean()))
    ax[0].plot(res['t'], res['y_meas'], 'C%d' % i, lw=1, label=name)
    ax[1].plot(res['t'], res['u'], 'C%d' % i, lw=1)
ax[0].axvline(900, color='k', lw=.6, ls=':'); ax[0].text(910, amb + 3, 'fan on')
ax[0].set_ylabel('T [degC]'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8, loc='lower right')
ax[0].set_title('Same PI controller on the model it was designed for and on the true plant')
ax[1].set_ylabel('heater [%]'); ax[1].set_xlabel('time [s]'); ax[1].grid(alpha=.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## PID: second zero at the Pade pole (IMC-style), derivative filter $T_f = \theta/10$
# $C(s) = K_c \frac{(s + z_1)(s + z_2)}{s (T_f s + 1)}$ with $z_1 = 1/\tau$, $z_2 = 2/\theta$. `sim.PID` realizes
# $K_p + K_i/s + K_d s/(T_f s + 1)$, whose numerator over $s(T_f s + 1)$ is
# $(K_p T_f + K_d)s^2 + (K_p + K_i T_f)s + K_i$, so matching coefficients gives $K_i = K_c z_1 z_2$,
# $K_p = K_c(z_1 + z_2) - K_i T_f$ and $K_d = K_c - K_p T_f$. Faster on paper, but the derivative turns the
# 0.1 degC sensor steps into heater chatter.

# %%
Tf = fit['theta'] / 10
z2 = 2 / fit['theta']
L2 = ct.tf(np.polymul([1, z], [1, z2]), np.polymul([1, 0], [Tf, 1])) * Gp
Kc2, cl2 = gain_for_zeta(L2, ZETA)
Ki2 = Kc2 * z * z2
Kp2 = Kc2 * (z + z2) - Ki2 * Tf
Kd2 = Kc2 - Kp2 * Tf
print('PID: Kp = %.3f, Ki = %.5f, Kd = %.2f %%s/degC, Tf = %.2f s; slowest complex pair %s'
      % (Kp2, Ki2, Kd2, Tf, np.round(cl2[np.argmax(np.where(np.abs(cl2.imag) > 1e-9, cl2.real, -9))], 4)))
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for name, ctrl, c in (('PI', sim.PID(Kp, Ki, Ts=TS, u_lim=U_LIM), 'C0'),
                      ('PID', sim.PID(Kp2, Ki2, Kd2, Ts=TS, Tf=Tf, u_lim=U_LIM), 'C2')):
    res = sim.simulate_loop(true_plant(p), ctrl, TS, 1200, lambda t: 35.0 if t < 600 else 40.0, u_lim=U_LIM,
                            noise_std=NOISE, quant=QUANT, delay=p['delay'], y0=amb, seed=33)
    m = res['t'] > 1000
    print('%-3s: heater std in steady state = %.2f %%' % (name, res['u'][m].std()))
    ax[0].plot(res['t'], res['y_meas'], c, lw=1, label=name)
    ax[1].plot(res['t'], res['u'], c, lw=.8, label=name)
ax[0].set_ylabel('T [degC]'); ax[0].grid(alpha=.3); ax[0].legend()
ax[0].set_title('PI vs PID on the true board: 35 -> 40 degC set-point change')
ax[1].set_ylabel('heater [%]'); ax[1].set_xlabel('time [s]'); ax[1].grid(alpha=.3)
fig.tight_layout()
plt.show()
