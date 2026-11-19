# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 02, demo 4: a steep hill saturates the throttle, and the integrator winds up (L19, L26)
#
# ![scene](figures/scene.png)
#
# **Overview.** On a 9 % grade the sedan needs more force than the engine gives in top gear, so the throttle
# sits at 100 % and the speed drops whatever the controller says. A plain PI integrator keeps accumulating
# the error during the climb; at the crest it is far too large and the car overshoots the set speed for a long
# time. Clamping anti-windup stops integrating while the throttle is saturated.
#
# **Math.** Needed throttle on a grade: $u = u_0 + 100\,m g\sin\vartheta/F_{max}$; above 100 % the loop is open
# and $\dot I = K_i e$ grows without bound. We finish with the Bode margins of the final loop.

# %%
import os, sys
import numpy as np
import control as ct
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ctrlprac import data, ident, sim
from _model import TS, U_LIM, CAR_A, linearize, grade_to_pct, gain_for_zeta, margins, drive, recovery

df, meta = data.load('02_cruise_step')
fit = ident.fit_fopdt(df['time_s'], df['speed_wheel_mps'], 5.0, 20.0)
car, lin = CAR_A, linearize(CAR_A)
u0, vset, HILL = round(lin['u0'], 1), car['v0'], 0.09
Kc, _ = gain_for_zeta(fit['model'] * ct.tf(*ct.pade(fit['theta'], 1)), 0.7)
z = 0.05                                   # PI zero of demo 3 (rounded)
print('PI: Kp = %.1f %%/(m/s), Ki = %.2f %%/(m s)' % (Kc, Kc * z))
print('9 %% grade needs %.0f %% throttle at %.0f m/s (limit 100 %%)' % (u0 - grade_to_pct(car, HILL), vset))

# %% [markdown]
# ## A 60 s climb at 9 % (ramps of 5 s), PI with and without anti-windup

# %%
hill = lambda t: HILL * np.clip(min((t - 20) / 5, (85 - t) / 5), 0, 1)
fig, ax = plt.subplots(3, 1, figsize=(8, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 2, 1]))
for aw, c in ((True, 'C0'), (False, 'C3')):
    res = drive(car, 260, controller=sim.PID(Kc, Kc * z, Ts=TS, u_lim=U_LIM, u0=u0, antiwindup=aw),
                r=vset, grade=hill, seed=50)
    m = res['t'] >= 85
    over = res['v'][m].max() - vset
    _, t_rec = recovery(res['t'], res['v'], vset, 85.0, band=0.3)
    print('anti-windup %-5s: min speed %.2f m/s, overshoot after crest %.2f m/s, back within 0.3 m/s %.0f s after crest,'
          ' throttle at 100 %% for %.0f s' % (aw, res['v'].min(), over, t_rec, TS * np.sum(res['u'] >= 99.99)))
    lab = 'clamping anti-windup' if aw else 'no anti-windup'
    ax[0].plot(res['t'], res['v'], c, label=lab)
    ax[1].plot(res['t'], res['u'], c, lw=1)
ax[2].plot(res['t'], 100 * res['grade'], 'k')
ax[0].axhline(vset, color='k', lw=.6); ax[0].set_ylabel('speed [m/s]'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
ax[0].set_title('Steep climb: the throttle saturates; windup causes a long overshoot after the crest')
ax[1].set_ylabel('throttle [%]'); ax[1].grid(alpha=.3)
ax[2].set_ylabel('grade [%]'); ax[2].set_xlabel('time [s]'); ax[2].grid(alpha=.3)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Margins of the final loop: identified model vs the linearized "true" car

# %%
C = Kc * ct.tf([1, z], [1, 0])
w = np.logspace(-3, 1, 2000)
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for (name, G, dly), c in zip((('identified FOPDT', fit['model'], fit['theta']),
                               ('true car linearized at 27 m/s', lin['G'], car['delay'])), ('C0', 'C3')):
    m = margins(C * G, dly)
    print('%-30s: wc = %.3f rad/s, PM = %.1f deg, GM = %.2f (%.1f dB), delay margin = %.2f s'
          % (name, m['wc'], m['PM'], m['GM'], 20 * np.log10(m['GM']), m['DM']))
    H = np.asarray((C * G)(1j * w)).ravel() * np.exp(-1j * w * dly)
    ax[0].semilogx(w, 20 * np.log10(abs(H)), c, label=name)
    ax[1].semilogx(w, np.degrees(np.unwrap(np.angle(H))), c, label='%s: PM = %.0f deg' % (name, m['PM']))
ax[0].axhline(0, color='k', lw=.6); ax[1].axhline(-180, color='k', lw=.6)
ax[0].set_ylabel('|L| [dB]'); ax[0].grid(alpha=.3, which='both'); ax[0].legend(fontsize=8)
ax[0].set_title('Loop gain of the cruise PI')
ax[1].set_ylabel('phase [deg]'); ax[1].set_xlabel('frequency [rad/s]'); ax[1].set_ylim(-300, -60)
ax[1].grid(alpha=.3, which='both'); ax[1].legend(fontsize=8)
fig.tight_layout()
plt.show()
