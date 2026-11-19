# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 02, demo 2: a hill is a disturbance, and P control cannot fully reject it (L2, L12, L13)
#
# ![scene](figures/scene.png)
#
# **Overview.** When the road tilts up, gravity pulls back with $m g \sin\vartheta$. In the block diagram this
# force enters at the plant input, like a throttle reduction. A proportional cruise controller answers with
# more throttle, but only in proportion to the speed error, so the car settles below the set speed.
#
# **Math.** With trim throttle $u_0$ and $u = u_0 + K_p(r - v)$, the grade adds
# $d = -100\, m g \sin\vartheta / F_{max}$ (% throttle) at the plant input:
# $V = \dfrac{K_p G}{1 + K_p G} R + \dfrac{G}{1 + K_p G} D$, so the steady speed drop is
# $\Delta v_{ss} = \dfrac{K |d|}{1 + K K_p}$. Higher $K_p$ helps, but the root locus shows how far we can go.

# %%
import os, sys
import numpy as np
import control as ct
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ctrlprac import data, ident, sim
from _model import TS, U_LIM, CAR_A, linearize, grade_to_pct, gain_for_zeta, margins, drive

df, meta = data.load('02_cruise_step')
fit = ident.fit_fopdt(df['time_s'], df['speed_wheel_mps'], 5.0, 20.0)
car, lin = CAR_A, linearize(CAR_A)
u0, vset, GRADE = round(lin['u0'], 1), car['v0'], 0.03
d = grade_to_pct(car, GRADE)
print('3 %% grade: m g sin = %.0f N = %.1f %% throttle; identified K = %.3f (m/s)/%%' % (-d * car['Fmax'] / 100, d, fit['K']))

# %% [markdown]
# ## How large can $K_p$ be? Root locus of $K_p G(s)$ with a Pade(1) model of the fitted dead time

# %%
Gp = fit['model'] * ct.tf(*ct.pade(fit['theta'], 1))
Kp7, cl = gain_for_zeta(Gp, 0.7)
mP = margins(Kp7 * fit['model'], fit['theta'])
print('Kp for zeta = 0.7: %.1f %%/(m/s), poles %s; PM = %.0f deg, GM = %.2f'
      % (Kp7, np.round(cl, 3), mP['PM'], mP['GM']))
num, den = np.ravel(Gp.num[0][0]), np.ravel(Gp.den[0][0])
br = np.array([np.roots(np.polyadd(den, k * num)) for k in np.logspace(-2, 2.3, 1500)])
fig, ax = plt.subplots(figsize=(6.5, 5))
ax.plot(br.real.ravel(), br.imag.ravel(), '.', ms=1.5, color='C0', label='locus, Kp > 0')
ax.plot(np.roots(den).real, np.roots(den).imag, 'kx', ms=9, label='poles')
ax.plot(np.roots(num).real, np.roots(num).imag, 'ko', mfc='none', ms=9, label='Pade zero (RHP)')
ax.plot(cl.real, cl.imag, 'C3s', label='Kp = %.1f (zeta 0.7)' % Kp7)
r_ = np.array([0, 1.2])
ax.plot(-0.7 * r_, 0.714 * r_, 'k:', lw=.8); ax.plot(-0.7 * r_, -0.714 * r_, 'k:', lw=.8)
ax.set_xlim(-1.0, 1.0); ax.set_ylim(-0.9, 0.9); ax.axhline(0, color='k', lw=.5); ax.axvline(0, color='k', lw=.5)
ax.set_xlabel('Re(s) [rad/s]'); ax.set_ylabel('Im(s) [rad/s]'); ax.grid(alpha=.3); ax.legend(fontsize=8)
ax.set_title('P cruise control: the lag/delay limits the gain')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Simulation on the nonlinear car: flat road, then a 3 % grade from t = 20 s

# %%
grade = lambda t: GRADE * np.clip((t - 20) / 3, 0, 1)
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for i, Kp in enumerate([5.0, 10.0, Kp7]):
    res = drive(car, 150, controller=sim.PID(Kp, Ts=TS, u_lim=U_LIM, u0=u0), r=vset, grade=grade, seed=20 + i)
    drop = vset - res['v'][-100:].mean()
    pred = fit['K'] * abs(d) / (1 + fit['K'] * Kp)
    print('Kp = %5.1f: speed drop predicted %.2f m/s, simulated %.2f m/s' % (Kp, pred, drop))
    ax[0].plot(res['t'], res['v_wheel'], 'C%d' % i, lw=.8, label='Kp = %.1f: drop %.2f m/s' % (Kp, drop))
    ax[1].plot(res['t'], res['u'], 'C%d' % i, lw=.8)
ax[0].axhline(vset, color='k', lw=.6); ax[0].axvline(20, color='k', lw=.6, ls=':')
ax[0].set_ylabel('speed [m/s]'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
ax[0].set_title('P cruise control on a 3 % grade: a steady speed deficit remains')
ax[1].set_ylabel('throttle [%]'); ax[1].set_xlabel('time [s]'); ax[1].grid(alpha=.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# The deficit falls like $1/(1 + K K_p)$ but never reaches zero, and $K_p$ cannot grow much beyond the
# $\zeta = 0.7$ value without ringing (and amplifying sensor noise into the throttle). Demo 3 adds an integrator.
