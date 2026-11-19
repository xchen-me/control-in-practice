# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 03, demo 2: position loop, P control on the root locus and a PD design
#
# ![scene](figures/scene.png)
#
# **Overview.** The arm must move 0.5 rad with overshoot at most 5 % and a 2 % settling time at most
# 0.15 s, without asking for more than the 12 V the PWM stage can give. Specs become a region of the
# s-plane (L16). P control cannot reach it, which the root locus shows at once (L13). A PD zero bends the
# locus into the region (L17), but the same zero appears in the closed loop and raises the overshoot (L9),
# so the first design is revised.

# %%
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.modules.pop('_model', None)          # run_all.py runs several labs in one process
import numpy as np
import matplotlib.pyplot as plt
import control as ct
from _model import identify_from_log, reduced_tf, pd_by_angle, step_specs

fit, _ = identify_from_log('03_motor_step')
K, tau = fit['K'], fit['tau']
G = reduced_tf(K, tau)
s = ct.tf('s')
r0, OS_max, Ts_max, V_max = 0.5, 5.0, 0.15, 12.0
print('design model G(s) = %.3f / (s (%.4f s + 1)), open-loop poles 0 and %.1f' % (K, tau, -1 / tau))

# %% [markdown]
# ## Specs as a region (L16)
# $OS = e^{-\pi\zeta/\sqrt{1-\zeta^2}} \le 5\% \Rightarrow \zeta \ge 0.69$;
# $T_s \approx 4/(\zeta\omega_n) \le 0.15$ s $\Rightarrow \sigma = \zeta\omega_n \ge 26.7$ 1/s.

# %%
zeta_min = -np.log(OS_max / 100) / np.sqrt(np.pi ** 2 + np.log(OS_max / 100) ** 2)
sigma_min = 4 / Ts_max
print('zeta >= %.3f, sigma >= %.1f 1/s' % (zeta_min, sigma_min))


def region(ax, ymax=60):
    """Shade the s-plane region that meets the specs."""
    th = np.arccos(zeta_min)
    x = np.array([-sigma_min, -ymax / np.tan(th) * 1.5, -ymax / np.tan(th) * 1.5])
    xs = np.linspace(-120, -sigma_min, 50)
    yb = np.minimum(-xs * np.tan(th), ymax)
    ax.fill_between(xs, -yb, yb, color='tab:green', alpha=0.12, label='meets specs')
    ax.axvline(-sigma_min, color='tab:green', lw=0.8)
    ax.plot([0, -ymax / np.tan(th)], [0, ymax], color='tab:green', lw=0.8)
    ax.plot([0, -ymax / np.tan(th)], [0, -ymax], color='tab:green', lw=0.8)


# %% [markdown]
# ## P control: $C(s)=K_p$ (L13)
# Closed loop: $\tau s^2 + s + KK_p = 0$. The two branches meet at $-1/(2\tau)$ and go straight up:
# $\sigma = 1/(2\tau)$ for every $K_p$, so $T_s \approx 8\tau$ cannot be shortened by gain.

# %%
gains = np.linspace(0, 40, 400)
loc = ct.root_locus_map(G, gains).loci
print('P control: sigma = 1/(2 tau) = %.1f 1/s -> Ts >= %.3f s (spec %.2f s): not reachable'
      % (1 / (2 * tau), 8 * tau, Ts_max))
Kp_z = 1 / (4 * K * tau * zeta_min ** 2)
print('largest Kp with zeta >= %.2f: Kp = %.2f V/rad' % (zeta_min, Kp_z))

# %% [markdown]
# ## PD by the angle condition (L17)
# Choose a target pole $s_d=-\sigma + j\sigma\tan(\cos^{-1}\zeta)$, then place the PD zero so that
# $\angle C(s_d)G(s_d) = 180^\circ$ and pick $K_d$ from $|C(s_d)G(s_d)|=1$.

# %%
designs = {}
for zeta in (0.70, 0.80):
    sigma = 27.0                     # just inside sigma_min, keeps Kp r near 12 V
    sd = complex(-sigma, sigma * np.tan(np.arccos(zeta)))
    Kd, z = pd_by_angle(K, tau, sd)
    C = Kd * (s + z)
    T = ct.feedback(C * G)
    tt = np.linspace(0, 0.4, 4001)
    _, y = ct.step_response(r0 * T, tt)
    _, y_nz = ct.step_response(r0 * ct.tf([T.den[0][0][-1]], T.den[0][0]), tt)   # same poles, no zero
    sp, sp_nz = step_specs(tt, y, r0), step_specs(tt, y_nz, r0)
    designs[zeta] = dict(Kd=Kd, z=z, sd=sd, T=T, y=y, y_nz=y_nz, t=tt)
    print('zeta=%.2f: sd = %.1f%+.1fj, zero at -%.1f, Kd = %.4f V s/rad, Kp = Kd z = %.2f V/rad'
          % (zeta, sd.real, sd.imag, z, Kd, Kd * z))
    print('   with the zero: OS = %.1f %%, Ts = %.3f s;  poles only: OS = %.1f %%, Ts = %.3f s;  proportional part Kp r = %.1f V'
          % (sp['OS'], sp['Ts'], sp_nz['OS'], sp_nz['Ts'], Kd * z * r0))

# %%
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
region(ax[0])
ax[0].plot(loc.real, loc.imag, 'b', lw=1.5)
ax[0].plot([0, -1 / tau], [0, 0], 'bx', ms=8, label='P: open-loop poles')
d = designs[0.80]
locpd = ct.root_locus_map(d['Kd'] * (s + d['z']) * G, np.linspace(0, 3, 600)).loci
ax[0].plot(locpd.real, locpd.imag, 'r', lw=1.5)
ax[0].plot(-d['z'], 0, 'ro', mfc='none', ms=8, label='PD zero (zeta=0.8)')
for zeta, m in ((0.70, 's'), (0.80, 'D')):
    sd = designs[zeta]['sd']
    ax[0].plot([sd.real] * 2, [sd.imag, -sd.imag], 'k' + m, mfc='none', label='target, zeta=%.2f' % zeta)
ax[0].set(xlim=(-120, 10), ylim=(-60, 60), xlabel='Re(s) (1/s)', ylabel='Im(s) (rad/s)',
          title='Root loci: P (blue) and PD (red)')
ax[0].grid(True)
ax[0].legend(fontsize=8, loc='lower left')
for zeta, col in ((0.70, 'tab:orange'), (0.80, 'tab:red')):
    d = designs[zeta]
    ax[1].plot(d['t'], d['y'], color=col, label='PD, zeta=%.2f' % zeta)
    ax[1].plot(d['t'], d['y_nz'], ':', color=col, label='same poles, no zero')
ax[1].axhline(r0 * (1 + OS_max / 100), color='k', lw=0.8, ls='--', label='5 % overshoot limit')
ax[1].axvline(Ts_max, color='0.5', lw=0.8, ls='--', label='Ts limit')
ax[1].set(xlabel='time (s)', ylabel='arm angle (rad)', title='0.5 rad step, design model')
ax[1].grid(True)
ax[1].legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# **Takeaways.** With zeta = 0.70 the poles alone would give 4.6 % overshoot, but the closed-loop zero at
# $-z$ (from $C(s)=K_d(s+z)$ acting on the error) pushes it above 5 %. Moving the target to zeta = 0.80
# moves the zero closer to the poles in angle but lowers the pole overshoot enough: the spec is met with
# about 3 % overshoot, and the initial voltage $K_p r$ stays under 12 V. Next (demo 3): a pure PD cannot
# be built, and differentiating encoder counts is noisy, so we replace it by a lead.
