# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 01, demo 4: the delay limits phase margin and bandwidth (L20-L26)
#
# ![scene](figures/scene.png)
#
# **Overview.** A dead time does not change $|L(j\omega)|$ but adds phase lag $-\omega\theta$ that grows
# without bound. Pushing the crossover frequency up therefore eats the phase margin, and a loop with dead
# time $\theta$ cannot have a crossover much above roughly $1/\theta$. We read this off the Bode plot,
# see it on the Nyquist plot, and compare the identified model with the true plant.
#
# **Math.** $L(j\omega) = C(j\omega)G(j\omega)e^{-j\omega\theta}$, $\angle e^{-j\omega\theta} = -\omega\theta$ (rad).
# $PM = 180^\circ + \angle L(j\omega_c)$ at $|L(j\omega_c)| = 1$; $GM = 1/|L(j\omega_{pc})|$ at
# $\angle L = -180^\circ$. Extra delay the loop tolerates (delay margin): $\Delta\theta = PM\,[\mathrm{rad}]/\omega_c$.

# %%
import os, sys
import numpy as np
import control as ct
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ctrlprac import data, ident
from _model import true_plant, fopdt_pade, gain_for_zeta, margins, freq

df, meta = data.load('01_thermal_step')
fit = ident.fit_fopdt(df['time_s'], df['T_C'], 50.0, 60.0)
p = meta['parameters']
z = 1 / fit['tau']
Kc, _ = gain_for_zeta(ct.tf([1, z], [1, 0]) * fopdt_pade(fit), 0.7)   # the PI of demo 3
C = Kc * ct.tf([1, z], [1, 0])
L_mod, L_true = C * fit['model'], C * true_plant(p)
w = np.logspace(-4, 0, 3000)

loops = [('identified FOPDT, exact delay %.1f s' % fit['theta'], L_mod, fit['theta']),
         ('true plant, delay %.0f s' % p['delay'], L_true, p['delay'])]
print('PI: Kp = %.3f %%/degC, Ki = %.5f %%/(degC s)' % (Kc, Kc * z))
for name, L, d in loops:
    m = margins(L, d)
    print('%-36s: wc = %.4f rad/s, PM = %.1f deg, w180 = %.4f rad/s, GM = %.2f (%.1f dB), delay margin = %.1f s'
          % (name, m['wc'], m['PM'], m['wpc'], m['GM'], 20 * np.log10(m['GM']), m['DM']))

# %% [markdown]
# ## Bode plot: same magnitude, the delay only adds phase

# %%
fig, ax = plt.subplots(2, 1, figsize=(8, 6.5), sharex=True)
for (name, L, d), c in zip(loops, ('C0', 'C3')):
    H = freq(L, w, d)
    ax[0].semilogx(w, 20 * np.log10(np.abs(H)), c, label=name)
    ax[1].semilogx(w, np.degrees(np.unwrap(np.angle(H))), c, label=name)
    m = margins(L, d)
    ax[1].plot(m['wc'], m['PM'] - 180, c + 'o')
    ax[1].annotate('PM = %.0f deg' % m['PM'], (m['wc'], m['PM'] - 180), (m['wc'] * 1.4, m['PM'] - 150),
                   color=c, fontsize=8, arrowprops=dict(arrowstyle='->', color=c))
H0 = freq(L_mod, w, 0.0)
ax[1].semilogx(w, np.degrees(np.unwrap(np.angle(H0))), 'k:', label='identified model without delay')
ax[0].axhline(0, color='k', lw=.6); ax[1].axhline(-180, color='k', lw=.6)
ax[0].set_ylabel('|L| [dB]'); ax[0].grid(alpha=.3, which='both'); ax[0].legend(fontsize=8)
ax[0].set_title('Loop gain L = C G e^{-s theta} with the PI of demo 3')
ax[1].set_ylabel('phase of L [deg]'); ax[1].set_xlabel('frequency [rad/s]'); ax[1].set_ylim(-400, -60)
ax[1].grid(alpha=.3, which='both'); ax[1].legend(fontsize=8, loc='lower left')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Achievable bandwidth: phase margin vs crossover when we raise the PI gain
# Scaling $K_c$ moves $\omega_c$; the delay costs $\omega_c\theta$ rad at crossover. We also ask what a
# shorter apparent delay (e.g. thermistor glued closer to the heater) would buy.

# %%
fig, ax = plt.subplots(figsize=(7, 4.5))
for th, c in ((fit['theta'], 'C0'), (fit['theta'] / 2, 'C2'), (p['delay'], 'C1')):
    wcs, pms = [], []
    for a in np.logspace(-1, 1.2, 120):
        m = margins(a * L_mod, th)
        wcs.append(m['wc']); pms.append(m['PM'])
    wcs, pms = np.array(wcs), np.array(pms)
    ok = pms >= 45
    wmax = wcs[ok].max()
    print('theta = %5.1f s: largest crossover with PM >= 45 deg is %.4f rad/s (= %.2f / theta)'
          % (th, wmax, wmax * th))
    ax.semilogx(wcs, pms, c, label='theta = %.1f s: wc,max = %.3f rad/s' % (th, wmax))
    ax.plot(wmax, 45, c + 'o')
mk = margins(L_mod, fit['theta'])
ax.plot(mk['wc'], mk['PM'], 'ks', label='design of demo 3')
ax.axhline(45, color='k', lw=.6, ls=':')
ax.set_xlabel('gain crossover frequency wc [rad/s]'); ax.set_ylabel('phase margin [deg]')
ax.set_ylim(-60, 100); ax.grid(alpha=.3, which='both'); ax.legend(fontsize=8)
ax.set_title('PI zero fixed at 1/tau, gain scaled: PM falls as wc grows')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Nyquist plot: the delay winds $L(j\omega)$ into a spiral around the origin
# The open loop is stable (integrator aside), so the closed loop is stable iff the plot does not encircle
# $-1$. Raising $K_c$ by more than GM swells the spiral until it encircles $-1$.

# %%
m = margins(L_mod, fit['theta'])
wn = np.logspace(-3.2, 0, 4000)
fig, ax = plt.subplots(figsize=(6.5, 6))
for a, c, lab in ((1.0, 'C0', 'design (GM = %.2f)' % m['GM']), (1.2 * m['GM'], 'C3', 'Kc x %.2f (unstable)' % (1.2 * m['GM']))):
    H = a * freq(L_mod, wn, fit['theta'])
    ax.plot(H.real, H.imag, c, lw=1.2, label=lab)
    ax.plot(H.real, -H.imag, c, lw=.6, ls='--')
th_ = np.linspace(0, 2 * np.pi, 200)
ax.plot(np.cos(th_), np.sin(th_), 'k:', lw=.8, label='unit circle')
ax.plot(-1, 0, 'k+', ms=12, mew=2)
Hc = freq(L_mod, np.array([m['wc']]), fit['theta'])[0]
ax.plot([0, Hc.real], [0, Hc.imag], 'C0', lw=.8)
ax.text(Hc.real - 0.2, Hc.imag - 0.25, 'PM = %.0f deg' % m['PM'], color='C0')
ax.set_xlim(-2.5, 1.5); ax.set_ylim(-2.5, 1.5); ax.set_aspect('equal')
ax.set_xlabel('Re L(jw)'); ax.set_ylabel('Im L(jw)'); ax.grid(alpha=.3); ax.legend(fontsize=8, loc='lower right')
ax.set_title('Nyquist of the identified loop with exact delay')
fig.tight_layout()
plt.show()

# %% [markdown]
# **Take-away.** The identified model (apparent dead time 18.6 s) is *conservative*: the true plant has only
# 4 s of transport delay and its sensor lag costs less phase near crossover, so the true PM is a little higher.
# The achievable crossover scales like $1/\theta$: halving the apparent delay roughly doubles the usable bandwidth.
