# Author: Xu Chen, University of Washington
# %% [markdown]
# # Belt drive 3: lead design on the Bode plot, margins checked on the measured response
#
# ![scene](figures/scene.png)
#
# **Overview.** We design a lead compensator on the low-order rigid-body model for a chosen crossover and
# phase margin (L21, L26). The design model ignores the resonance, so the decisive check is done on the
# MEASURED frequency response: phase and gain margins, the Nyquist plot with the resonance loop, and the
# sensitivity peak $M_s = \max|1/(1 + L)|$. Then we push the crossover toward the resonance and watch the
# margins, and the closed loop, fall apart.

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import control as ct
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
from ctrlprac import data, sim
import _model as M

df, meta = data.load('06_drive_prbs')
w, Gm, coh = M.measured_frf(df)
f = w / 2 / np.pi
Td = (meta['delay_samples'] + 0.5) * meta['sample_time_s']
band = (coh > 0.9) & (f >= 2) & (f <= 150)
fit = M.structured_fit(w[band], Gm[band], Td)
Glf, Gfull = fit['rigid'], fit['model']
print('design model G_lf = %.0f/(s(s + %.2f)); resonance at %.1f Hz, zeta %.3f; drive delay %.3f ms'
      % (fit['K'], fit['a'], fit['f_r'], fit['zeta'], 1e3 * Td))

# %% [markdown]
# ## Lead design on the Bode plot (L21)
# Target: crossover $f_c$ = 10 Hz (one sixth of the resonance) with PM = 50 deg. At $\omega_c$ the rigid body
# plus the delay has phase $-90^\circ - \arctan(\omega_c/a) - \omega_c T_d$; the lead must add the rest:
#
# $$ \phi_{max} = PM - 180^\circ - \angle G_{lf}(j\omega_c)e^{-j\omega_cT_d},\qquad
#    \alpha = \frac{1 - \sin\phi_{max}}{1 + \sin\phi_{max}},\quad z = \omega_c\sqrt\alpha,\ p = \omega_c/\sqrt\alpha $$
#
# and $K$ makes $|C(j\omega_c)G_{lf}(j\omega_c)| = 1$.

# %%
fc, PMt = 10.0, 50.0
C, dsg = M.design_lead(Glf, 2 * np.pi * fc, PMt, Td)
print('lead: phi_max = %.1f deg, zero %.1f Hz, pole %.1f Hz, K = %.4f A/mm'
      % (dsg['phi'], dsg['z'] / 2 / np.pi, dsg['p'] / 2 / np.pi, dsg['K']))
print('C(s) = %.4f (s/%.1f + 1)/(s/%.1f + 1)' % (dsg['K'], dsg['z'], dsg['p']))

# %% [markdown]
# ## Margins: design model vs MEASURED response (L26)
# $L_{meas}(j\omega) = C(j\omega)\,\hat G(j\omega)$ uses the ETFE directly (the drive delay is already in it).
# Coherence limits the trusted band to about 150 Hz; beyond it the fitted model (with delay) fills in
# for the Nyquist plot.

# %%
wm = w[band]
L_meas = C(1j * wm) * Gm[band]
wd = np.logspace(np.log10(2 * np.pi * 1), np.log10(np.pi * 4000), 6000)
L_lf = C(1j * wd) * Glf(1j * wd) * np.exp(-1j * wd * Td)
L_mod = C(1j * wd) * Gfull(1j * wd) * np.exp(-1j * wd * Td)


def report(tag, ww, L):
    m = M.frf_margins(ww, L)
    gm = ['%.1f dB at %.1f Hz' % (20 * np.log10(g), wx / 2 / np.pi) for wx, g in zip(m['w180'], m['GM'])]
    print('%-22s crossover %s Hz, PM %s deg | GM %s | Ms = %.2f (%.1f dB) at %.1f Hz'
          % (tag, np.round(np.array(m['wc']) / 2 / np.pi, 1), np.round(m['PM'], 1), ', '.join(gm) or 'none',
             m['Ms'], 20 * np.log10(m['Ms']), m['w_Ms'] / 2 / np.pi))
    return m


report('rigid-body model', wd, L_lf)
report('full fitted model', wd, L_mod)
mm = report('MEASURED response', wm, L_meas)

fig, ax = plt.subplots(2, 1, figsize=(8.5, 6.5), sharex=True)
fd = wd / 2 / np.pi
ax[0].semilogx(fd, 20 * np.log10(abs(L_lf)), 'C1--', label='design model (rigid body)')
ax[0].semilogx(f[band], 20 * np.log10(abs(L_meas)), 'C0', lw=1.5, label='measured G x C')
ax[0].axhline(0, color='k', lw=0.6)
ax[0].set(ylabel='|L| (dB)', title='Loop gain with the lead: design model vs measured', ylim=(-60, 50))
ph_m = np.degrees(np.unwrap(np.angle(L_meas)))
ph_m -= 360 * np.round((ph_m[0] + 120) / 360)
ph_l = np.degrees(np.unwrap(np.angle(L_lf)))
ph_l -= 360 * np.round((ph_l[0] + 120) / 360)
ax[1].semilogx(fd, ph_l, 'C1--', label='design model')
ax[1].semilogx(f[band], ph_m, 'C0', lw=1.5, label='measured')
ax[1].axhline(-180, color='k', lw=0.6)
ax[1].set(xlabel='frequency (Hz)', ylabel='phase (deg)', ylim=(-420, -90), xlim=(1, 300))
for x_, lab in ((mm['wc'][0] / 2 / np.pi, 'crossover'), (mm['w180'][0] / 2 / np.pi, 'phase crossover')):
    for a in ax:
        a.axvline(x_, color='0.5', ls=':', lw=0.8)
ax[1].annotate('PM = %.0f deg' % mm['PM'][0], xy=(mm['wc'][0] / 2 / np.pi, -180 + mm['PM'][0]),
               xytext=(3, -260), arrowprops=dict(arrowstyle='->'))
ax[0].annotate('GM = %.1f dB (resonance)' % (20 * np.log10(mm['GM'][0])),
               xy=(mm['w180'][0] / 2 / np.pi, -20 * np.log10(mm['GM'][0])), xytext=(70, 25),
               arrowprops=dict(arrowstyle='->'))
for a in ax:
    a.grid(True, which='both', alpha=0.4)
    a.legend(fontsize=8, loc='lower left')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Nyquist plot including the resonance (L24, L26)
# The open loop has no right-half-plane poles (P = 0; the integrator is indented around), so -1 must not be
# encircled. The resonance draws an extra loop: its size is the resonance peak times the loop gain
# there. The closest approach to -1 is $1/M_s$.

# %%
fig, ax = plt.subplots(1, 2, figsize=(11, 5))
for a, lim in ((ax[0], 3.0), (ax[1], 1.3)):
    a.plot(L_mod.real, L_mod.imag, color='0.6', lw=0.8, label='fitted model (beyond the data)')
    a.plot(L_meas.real, L_meas.imag, 'C0', lw=1.5, label='measured, w > 0')
    a.plot(L_meas.real, -L_meas.imag, 'C0--', lw=0.8, label='measured, w < 0')
    a.plot(-1, 0, 'r+', ms=14, mew=2)
    th = np.linspace(0, 2 * np.pi, 200)
    a.plot(-1 + np.cos(th) / mm['Ms'], np.sin(th) / mm['Ms'], 'r:', lw=0.8, label='circle radius 1/Ms')
    a.set(xlim=(-lim, 0.6 * lim), ylim=(-lim, lim), xlabel='Re L', ylabel='Im L')
    a.set_aspect('equal')
    a.grid(True, alpha=0.4)
ax[0].set_title('Nyquist, measured loop (lead, $f_c$ = 10 Hz)')
ax[1].set_title('Zoom: resonance loop near -1')
ax[0].legend(fontsize=8, loc='lower left')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Pushing the crossover toward the resonance
# The same procedure for 15 and 20 Hz: more lead gain at the resonance, so the resonance loop grows
# toward -1. Then the step responses of the sampled loop on the simulated stage (current limit, encoder
# quantization and noise, drive delay): `simulate_loop` with `LTIController(C)`. The simulator applies each
# command one sample late by itself, so the two samples of drive delay are `delay=(N_DELAY - 1) * TS`. The
# 0.01 A force ripple of the logs is left out; `_make_data.py` shows how to add it as the disturbance `d`.

# %%
res_all = {}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
for fcx, col in ((10, 'C0'), (15, 'C1'), (20, 'C3')):
    Cx, dx = M.design_lead(Glf, 2 * np.pi * fcx, PMt, Td)
    Lx = Cx(1j * wm) * Gm[band]
    mx = report('fc = %d Hz (measured)' % fcx, wm, Lx)
    Px = M.plant_ss(M.STAGE_A)
    cl = ct.poles(ct.feedback(Cx * M.plant_tf(M.STAGE_A) * ct.tf(*ct.pade(Td, 3)), 1))
    print('    closed-loop (true plant, Pade delay): max Re(pole) = %+.1f rad/s' % cl.real.max())
    r = sim.simulate_loop(Px, sim.LTIController(Cx, M.TS, u_lim=M.I_LIM), M.TS, 0.4,
                          lambda t: 1.0 if t >= 0.02 else 0.0, u_lim=M.I_LIM, noise_std=M.ENC_NOISE_MM,
                          quant=M.ENC_RES_MM, delay=(M.N_DELAY - 1) * M.TS)
    res_all[fcx] = (mx, r)
    stable = cl.real.max() < 0
    ax[0].plot(r['t'], r['y'], col, lw=1, label='$f_c$ = %d Hz: GM %.1f dB, %s'
               % (fcx, 20 * np.log10(min(mx['GM'])), 'Ms %.1f' % mx['Ms'] if stable else 'UNSTABLE'))
    ax[1].plot(Lx.real, Lx.imag, col, lw=1.2, label='$f_c$ = %d Hz' % fcx)
ax[0].set(xlabel='time (s)', ylabel='x (mm)', title='1 mm step, simulated stage', ylim=(-0.5, 2.5))
ax[1].plot(-1, 0, 'r+', ms=14, mew=2)
ax[1].set(xlim=(-3, 1), ylim=(-2.5, 2.5), xlabel='Re L', ylabel='Im L', title='Measured Nyquist, w > 0')
ax[1].set_aspect('equal')
for a in ax:
    a.grid(True, alpha=0.4)
    a.legend(fontsize=8, loc='upper right')
fig.tight_layout()
plt.show()

# %% [markdown]
# **Reading the result.** At 10 Hz the resonance loop stays well clear of -1 (GM about 10 dB, Ms < 2),
# and the step peaks at 1.31 A, the first output of the lead, well inside the 3 A limit. At 15 Hz the
# measured loop still has PM = 47 deg at its first crossover (16.1 Hz), but the resonance loop nearly
# touches -1 (GM under 2 dB, Ms above 10): the loop rings at the belt resonance. At 20 Hz the resonance loop
# encircles -1 and the stage oscillates at the resonance until the current limit holds it. Phase margin
# at the first crossover alone would have approved all three designs; the measured gain margin at the
# resonance and Ms do not. To go faster: a notch filter at the resonance, a stiffer belt, or a
# collocated (motor-side) feedback for the high-frequency part.
