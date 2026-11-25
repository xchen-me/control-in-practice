# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 02, demo 3: where to put the PI zero, and lag vs PI (L13-L18)
#
# ![scene](figures/scene.png)
#
# **Overview.** Adding an integrator, $C(s) = K_c (s + z)/s$, removes the speed deficit on a hill. The zero
# location $z$ is the design choice. Cancelling the slow vehicle pole ($z = 1/\tau$) looks tidy, but the
# grade enters *before* the plant, so the cancelled pole still shows up in the disturbance response as a
# slow 1-minute recovery. Moving $z$ left speeds recovery and costs phase margin. A lag compensator
# $K_c (s + z)/(s + p)$ is the "almost PI": it raises the low-frequency gain by $z/p$ but keeps an error.
#
# **Math.** Disturbance-to-speed: $\dfrac{V}{D} = \dfrac{G}{1 + CG}$. With $C = K_c(s+z)/s$ its DC gain is 0
# (type 1). The closed-loop pole that sits near $-z$ sets the tail of the recovery, $\approx e^{-z t}$.
# Lag: $e_{ss} = K|d|/(1 + K K_c z/p)$.
#
# The root loci and phase margins below use the FOPDT model fitted to the logged step, which does not change
# if you edit `car`; only the simulations use `car` and `lin`.

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
u0, vset, GRADE = round(lin['u0'], 1), car['v0'], 0.03
d = grade_to_pct(car, GRADE)
Gp = fit['model'] * ct.tf(*ct.pade(fit['theta'], 1))
Kc, _ = gain_for_zeta(Gp, 0.7)            # the P design of demo 2; the PI keeps this gain
print('Kc = %.1f %%/(m/s) (P design, zeta = 0.7); slow plant pole 1/tau = %.4f rad/s' % (Kc, 1 / fit['tau']))

# %% [markdown]
# ## Sweep the PI zero at fixed $K_c$: phase margin vs recovery from a 3 % grade (linear model)

# %%
zs = np.logspace(np.log10(1 / fit['tau']), np.log10(0.2), 25)
rows = []
for z in zs:
    C = Kc * ct.tf([1, z], [1, 0])
    m = margins(C * fit['model'], fit['theta'])
    p = ct.poles(ct.feedback(C * Gp, 1))
    res = sim.simulate_loop(lin['G'], sim.PID(Kc, Kc * z, Ts=TS, u_lim=U_LIM, u0=u0), TS, 300, vset,
                            d=lambda t: -u0 + (d if t >= 20 else 0.0), u_lim=U_LIM, delay=car['delay'], y0=vset)
    rows.append((z, m['PM'], p.real.max(), *recovery(res['t'], res['y'], vset, 20.0, band=0.2)))
rows = np.array(rows)
ok = rows[:, 1] >= 50
z_star = rows[ok, 0].max()
print('   z [rad/s]  PM [deg]  slowest pole  max drop [m/s]  recovery to 0.2 m/s [s]')
for rw in rows[::4]:
    print('   %8.4f   %6.1f    %8.4f      %6.2f        %6.1f' % tuple(rw))
print('design: largest z with PM >= 50 deg -> z = %.4f rad/s, Ki = Kc z = %.3f %%/m' % (z_star, Kc * z_star))

fig, ax = plt.subplots(figsize=(7, 4))
ax.semilogx(rows[:, 0], rows[:, 4], 'C0o-', ms=3, label='recovery time [s]')
ax.set_xlabel('PI zero z [rad/s]'); ax.set_ylabel('recovery to +-0.2 m/s [s]', color='C0'); ax.grid(alpha=.3, which='both')
ax2 = ax.twinx(); ax2.semilogx(rows[:, 0], rows[:, 1], 'C3s-', ms=3); ax2.set_ylabel('phase margin [deg]', color='C3')
ax.axvline(1 / fit['tau'], color='k', lw=.6, ls=':'); ax.text(1 / fit['tau'] * 1.05, rows[:, 4].max() * 0.9, 'z = 1/tau')
ax.axvline(z_star, color='C2', lw=1); ax.text(z_star * 1.05, rows[:, 4].max() * 0.6, 'design', color='C2')
ax.set_title('PI zero placement: faster hill recovery costs phase margin')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Root loci for three zero locations (gain $K_c$ marked)

# %%
fig, axs = plt.subplots(1, 3, figsize=(12, 4.2), sharey=True)
for ax, z, lab in zip(axs, (1 / fit['tau'], z_star, 0.15), ('cancel 1/tau', 'design', 'too far left')):
    L = ct.tf([1, z], [1, 0]) * Gp
    num, den = np.ravel(L.num[0][0]), np.ravel(L.den[0][0])
    br = np.array([np.roots(np.polyadd(den, k * num)) for k in np.logspace(-2, 2.3, 1200)])
    cl = np.roots(np.polyadd(den, Kc * num))
    ax.plot(br.real.ravel(), br.imag.ravel(), '.', ms=1.2, color='C0')
    ax.plot(np.roots(den).real, np.roots(den).imag, 'kx', ms=8)
    ax.plot(np.roots(num).real, np.roots(num).imag, 'ko', mfc='none', ms=8)
    ax.plot(cl.real, cl.imag, 'C3s', ms=6)
    ax.set_xlim(-0.9, 0.15); ax.set_ylim(-0.6, 0.6); ax.grid(alpha=.3)
    ax.axhline(0, color='k', lw=.5); ax.axvline(0, color='k', lw=.5)
    ax.set_title('%s: z = %.3f rad/s\nslowest pole %.3f rad/s' % (lab, z, cl.real.max()), fontsize=9)
    ax.set_xlabel('Re(s) [rad/s]')
axs[0].set_ylabel('Im(s) [rad/s]')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## P vs lag vs PI on the nonlinear car, 3 % grade at t = 20 s

# %%
p_lag = z_star / 10
ctrls = {'P': sim.PID(Kc, Ts=TS, u_lim=U_LIM, u0=u0),
         'lag (z/p = 10)': sim.LTIController(Kc * ct.tf([1, z_star], [1, p_lag]), TS, u_lim=U_LIM, u0=u0),
         'PI, z = 1/tau': sim.PID(Kc, Kc / fit['tau'], Ts=TS, u_lim=U_LIM, u0=u0),
         'PI, design z': sim.PID(Kc, Kc * z_star, Ts=TS, u_lim=U_LIM, u0=u0)}
grade = lambda t: GRADE * np.clip((t - 20) / 3, 0, 1)
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
for i, (name, c) in enumerate(ctrls.items()):
    res = drive(car, 300, controller=c, r=vset, grade=grade, seed=40)
    drop, t_rec = recovery(res['t'], res['v'], vset, 20.0, band=0.2)
    print('%-15s: max drop %.2f m/s, error at 300 s %.2f m/s, recovery to 0.2 m/s %s'
          % (name, drop, vset - res['v'][-50:].mean(), 'never' if np.isnan(t_rec) else '%.0f s' % t_rec))
    ax[0].plot(res['t'], res['v'], 'C%d' % i, lw=1.2, label=name)
    ax[1].plot(res['t'], res['u'], 'C%d' % i, lw=.8)
print('lag predicted e_ss = %.2f m/s' % (fit['K'] * abs(d) / (1 + fit['K'] * Kc * z_star / p_lag)))
ax[0].axhline(vset, color='k', lw=.6)
ax[0].set_ylabel('true speed [m/s]'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
ax[0].set_title('3 %% grade: P keeps an error, lag shrinks it slowly, PI removes it (Kc = %.1f)' % Kc)
ax[1].set_ylabel('throttle [%]'); ax[1].set_xlabel('time [s]'); ax[1].grid(alpha=.3)
fig.tight_layout()
plt.show()
