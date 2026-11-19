# Author: Xu Chen, University of Washington
# %% [markdown]
# # Balancing robot 2: PD design on the root locus, and the Bode plot of an unstable plant
#
# ![scene](figures/scene.png)
#
# **Overview.** With the model identified in demo 1, $G(s) = b_0/(s^2 - p^2)$, we place the closed-loop
# poles with a PD controller using the root locus (L13-L17): the PD zero bends the locus of the two
# real poles into the left half plane. Then we look at the Bode plot of the unstable plant (L22 Part 3):
# its magnitude is identical to that of a stable plant, but its phase is not, so magnitude alone can
# not tell you about stability.

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import control as ct
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
import _model as M

p, b0, meta = M.identified_A()
Tf = meta['test']['Tf_s']
G = ct.tf([b0], [1, 0, -p ** 2])
print('identified model: G(s) = %.0f / (s^2 - %.1f), poles at +/- %.2f rad/s' % (b0, p ** 2, p))

# %% [markdown]
# ## Design target and PD gains
# Target: $\zeta = 0.7$ (about 5 % overshoot) and $\omega_n = 2p$ (the loop must be clearly faster than the
# fall, whose time constant is $1/p$). With an ideal PD, $C(s) = K_d(s + z)$, $z = K_p/K_d$, the
# characteristic equation is $s^2 + b_0K_d s + (b_0K_p - p^2) = s^2 + 2\zeta\omega_n s + \omega_n^2$:
#
# $$ K_d = \frac{2\zeta\omega_n}{b_0}, \qquad K_p = \frac{\omega_n^2 + p^2}{b_0} $$

# %%
zeta, wn = 0.7, 2 * p
Kp, Kd = M.design_pd(p, b0, zeta, wn)
z = Kp / Kd
sd = -zeta * wn + 1j * wn * np.sqrt(1 - zeta ** 2)
print('target poles s = %.1f +/- j%.1f rad/s' % (sd.real, sd.imag))
print('PD gains: Kp = %.3f N m/rad, Kd = %.4f N m s/rad, zero at s = -%.1f rad/s' % (Kp, Kd, z))
print('(firmware gains in the log: Kp = %.3f, Kd = %.4f)' % (meta['test']['Kp_Nm_per_rad'],
                                                             meta['test']['Kd_Nms_per_rad']))
# angle condition check (L13): angle of (s+z)/(s^2-p^2) at s_d must be 180 deg
ang = np.degrees(np.angle((sd + z) / (sd ** 2 - p ** 2)))
print('angle condition at s_d: %.1f deg; magnitude condition gives Kd = %.4f'
      % (ang, 1 / abs(b0 * (sd + z) / (sd ** 2 - p ** 2))))

# %% [markdown]
# ## Root locus of the PD loop (zero fixed, gain Kd varied)
# Open loop $K_d\,b_0 (s + z)/(s^2 - p^2)$: two real poles at $\pm p$, one zero at $-z$. The real axis
# segment between $+p$ and $-p$ belongs to the locus; the branches break away, circle around the zero and
# re-enter the real axis to the left of it. The firmware's derivative filter adds a pole at $-1/T_f$
# (L15): it pulls the far branch away but barely moves the dominant poles.

# %%
K = np.r_[0, np.logspace(-4, 0, 1500)]
ideal = np.array([np.roots(np.polyadd([1, 0, -p ** 2], k * b0 * np.array([1, z]))) for k in K])
den_f = np.polymul([Tf, 1], [1, 0, -p ** 2])
# filtered PD with the same zero ratio: Kd s/(Tf s+1) + Kp = Kd ((1+z Tf) s + z)/(Tf s + 1)
filt = np.array([np.roots(np.polyadd(den_f, k * b0 * np.array([1 + z * Tf, z]))) for k in K])
cl_design = np.roots(np.polyadd(den_f, Kd * b0 * np.array([1 + z * Tf, z])))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
for a, roots, title in ((ax[0], ideal, 'ideal PD'), (ax[1], filt, 'PD with derivative filter $T_f$')):
    for j in range(roots.shape[1]):
        a.plot(roots[:, j].real, roots[:, j].imag, '.', ms=1.5, color='C0')
    a.plot([p, -p], [0, 0], 'kx', ms=9, label='open-loop poles')
    a.plot(-z, 0, 'ko', mfc='none', ms=8, label='PD zero $-z$')
    a.plot([sd.real] * 2, [sd.imag, -sd.imag], 'r*', ms=11, label='target $\\zeta$ = 0.7, $\\omega_n$ = 2p')
    ww = np.linspace(0, 1.5 * wn, 50)
    a.plot(-zeta * ww, ww * np.sqrt(1 - zeta ** 2), 'r:', lw=0.8)
    a.plot(-zeta * ww, -ww * np.sqrt(1 - zeta ** 2), 'r:', lw=0.8)
    a.axvline(0, color='k', lw=0.6)
    a.set(xlabel='Re(s) (rad/s)', ylabel='Im(s) (rad/s)', title='Root locus, ' + title)
    a.grid(True, alpha=0.4)
ax[0].set_xlim(-4 * p, 2 * p)
ax[0].set_ylim(-3 * p, 3 * p)
ax[1].plot(-1 / Tf, 0, 'kx', ms=9)
ax[1].plot(cl_design.real, cl_design.imag, 'gs', mfc='none', ms=9, label='closed loop at design Kd')
ax[1].set_xlim(-1.2 / Tf, 2 * p)
ax[1].set_ylim(-3 * p, 3 * p)
ax[0].legend(fontsize=8, loc='upper left')
ax[1].legend(fontsize=8, loc='upper left')
fig.tight_layout()
plt.show()
print('closed-loop poles with the filter pole (Tf = %.3f s):' % Tf, np.round(cl_design, 1))

# %% [markdown]
# ## Bode plot of the unstable plant (L22 Part 3)
# Compare $G(s) = b_0/(s^2 - p^2) = b_0/((s - p)(s + p))$ with the stable $G_s(s) = b_0/(s + p)^2$.
# Since $|j\omega - p| = |j\omega + p|$, the magnitudes are identical. The phases are not: the RHP pole
# contributes $-180^\circ + \arctan(\omega/p)$ instead of $-\arctan(\omega/p)$, so $\angle G = -180^\circ$ at
# every frequency, while $\angle G_s$ goes from $0^\circ$ to $-180^\circ$. A Bode magnitude plot can not
# reveal an unstable pole; the phase (and later the Nyquist plot with P = 1) does.

# %%
w = np.logspace(0, 3, 600)
Gs = ct.tf([b0], [1, 2 * p, p ** 2])
Gj, Gsj = G(1j * w), Gs(1j * w)
C = M.pd_tf(Kp, Kd, Tf)
Lj = C(1j * w) * Gj
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
ax[0].semilogx(w, 20 * np.log10(abs(Gj)), lw=2.5, label='$G = b_0/(s^2-p^2)$ (unstable)')
ax[0].semilogx(w, 20 * np.log10(abs(Gsj)), 'C1--', label='$G_s = b_0/(s+p)^2$ (stable)')
ax[0].semilogx(w, 20 * np.log10(abs(Lj)), 'C2', label='loop $L = C G$ (design PD)')
ax[0].axhline(0, color='k', lw=0.6)
ax[0].set_ylabel('magnitude (dB)')
ax[1].semilogx(w, np.degrees(np.unwrap(np.angle(Gj))) - 360 * (np.angle(Gj[0]) > 0), lw=2.5, label='G')
ax[1].semilogx(w, np.degrees(np.unwrap(np.angle(Gsj))), 'C1--', label='$G_s$')
ax[1].semilogx(w, np.degrees(np.unwrap(np.angle(Lj))) - 360 * (np.angle(Lj[0]) > 0), 'C2', label='L')
ax[1].axhline(-180, color='k', lw=0.6)
ax[1].set(xlabel='frequency (rad/s)', ylabel='phase (deg)')
ax[1].set_yticks([-270, -180, -90, 0])
for a in ax:
    a.axvline(p, color='0.5', ls=':', lw=0.8)
    a.grid(True, which='both', alpha=0.4)
    a.legend(fontsize=8)
ax[0].set_title('Same magnitude, different phase: the unstable pitch dynamics')
fig.tight_layout()
plt.show()
print('|G(j0)| = |Gs(j0)| = %.2f rad/(N m) = %.1f dB; phase G = %.0f deg, Gs = %.0f deg at all w -> 0'
      % (abs(G(0)), 20 * np.log10(abs(G(0))), -180, 0))
