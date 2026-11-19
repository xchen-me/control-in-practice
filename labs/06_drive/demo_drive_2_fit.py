# Author: Xu Chen, University of Washington
# %% [markdown]
# # Belt drive 2: fitting a transfer function to the measured response, and a low-order design model
#
# ![scene](figures/scene.png)
#
# **Overview.** The measured frequency response is the evidence; a transfer function is a compact
# summary we can design with. We fit one with `ident.fit_tf`, judge it (are its poles physical? where
# does it miss the data and why?), and then reduce it to the dominant low-frequency dynamics (L22):
# a rigid body $K/(s(s+a))$ that is accurate well below the resonance, which is where the position loop
# crossover must sit.

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import control as ct
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
from ctrlprac import data, ident
import _model as M

df, meta = data.load('06_drive_prbs')
w, G, coh = M.measured_frf(df)
f = w / 2 / np.pi
Td = (meta['delay_samples'] + 0.5) * meta['sample_time_s']      # known drive delay + ZOH
band = (coh > 0.9) & (f >= 2) & (f <= 150)
wb, Gb = w[band], G[band]
print('fitting %d ETFE points, %.0f-%.0f Hz, known delay %.3f ms removed' % (band.sum(), f[band].min(),
                                                                            f[band].max(), 1e3 * Td))

# %% [markdown]
# ## A first, unstructured fit
# Rigid body (2 poles) plus one resonance (2 poles) suggests $n_a = 4$ and, since the encoder is on the far
# side of the belt, no finite zeros in the band ($n_b = 0$). The fit is close in the band, but look at
# the poles.

# %%
F4 = ident.fit_tf(wb, Gb, 0, 4, delay=Td)
print('fit (nb=0, na=4) poles:', np.round(ct.poles(F4), 1))
print('-> a pole near s = 0 is pinned only by data well above it (the band starts at 12 rad/s); the fit')
print('   may put it in the right half plane, which is not physical for this stage.')

# %% [markdown]
# ## Put the known physics in: the integrator
# Position is the integral of velocity, so $G(s) = \frac{1}{s} G_v(s)$ exactly. Fit the velocity response
# $j\omega G(j\omega)$ with $n_b = 0$, $n_a = 3$:
#
# $$ G(s) = \frac{b_0}{s\,(s + a)\,(s^2 + 2\zeta\omega_r s + \omega_r^2)} $$

# %%
fit = M.structured_fit(wb, Gb, Td)
td = meta['true_derived']
Gfit = fit['model']
print('structured fit: a = %.2f rad/s (true %.2f), f_r = %.2f Hz (true %.2f), zeta = %.4f (true %.3f),'
      % (fit['a'], td['a_rigid_rad_s'], fit['f_r'], td['f_r_Hz'], fit['zeta'], td['zeta_r']))
print('                K = %.0f mm/(A s^2) (true %.0f, %+.1f %%)'
      % (fit['K'], td['K_rigid_mm_per_As2'], 100 * (fit['K'] / td['K_rigid_mm_per_As2'] - 1)))
Gf = Gfit(1j * w) * np.exp(-1j * w * Td)
err_db = 20 * np.log10(abs(Gf / G))
err_ph = np.degrees(np.angle(Gf / G))
for lo, hi in ((2, 30), (30, 90), (90, 150)):
    k = band & (f >= lo) & (f < hi)
    print('  %3d-%3d Hz: rms error %.2f dB, %.1f deg' % (lo, hi, np.sqrt(np.mean(err_db[k] ** 2)),
                                                       np.sqrt(np.mean(err_ph[k] ** 2))))

fig, ax = plt.subplots(3, 1, figsize=(8.5, 8), sharex=True, gridspec_kw=dict(height_ratios=[3, 3, 1.6]))
k = f > 1
ax[0].semilogx(f[band], 20 * np.log10(abs(G[band])), 'C0', lw=2, alpha=0.5, label='measured (coherence > 0.9)')
ax[0].semilogx(f[k], 20 * np.log10(abs(Gf[k])), 'k', lw=1, label='structured fit (with delay)')
ax[0].set(ylabel='|G| (dB, mm/A)', title='Fitted model vs measured frequency response', ylim=(-90, 45))
ph = np.degrees(np.unwrap(np.angle(G[band])))
ph -= 360 * np.round((ph[0] + 180) / 360)
phf = np.degrees(np.unwrap(np.angle(Gf[k])))
phf -= 360 * np.round((phf[np.argmin(abs(f[k] - f[band][0]))] - ph[0]) / 360)
ax[1].semilogx(f[band], ph, 'C0', lw=2, alpha=0.5, label='measured')
ax[1].semilogx(f[k], phf, 'k', lw=1, label='fit')
ax[1].set(ylabel='phase (deg)', ylim=(-480, -90))
ax[2].semilogx(f[band], err_db[band], 'C3', lw=0.8, label='magnitude error (dB)')
ax[2].semilogx(f[band], err_ph[band] / 10, 'C2', lw=0.8, label='phase error (deg/10)')
ax[2].axhline(0, color='k', lw=0.6)
ax[2].set(xlabel='frequency (Hz)', ylabel='fit - data', ylim=(-3, 3), xlim=(1, 300))
for a in ax:
    a.grid(True, which='both', alpha=0.4)
    a.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# **Judging the fit.** The errors are small and random below 90 Hz. Above, a slow phase drift appears:
# the 1 kHz current loop is not in the model; over the band it looks like a little extra delay. Both are
# far above any feasible crossover, so this model is good for design; the measured data remain the
# reference for the final margin check.
#
# ## Dominant-pole approximation (L22)
# Well below $\omega_r$ the resonance factor $\omega_r^2/(s^2 + 2\zeta\omega_r s + \omega_r^2) \approx 1$ (unit DC
# gain), so the dominant dynamics are the rigid body
#
# $$ G_{lf}(s) = \frac{K}{s\,(s + a)} $$
#
# The approximation error is $|1/(1 - (\omega/\omega_r)^2)|$: below 1 dB up to about $0.33 f_r$.

# %%
Glf = fit['rigid']
print('G_lf(s) = %.0f / (s (s + %.2f))   [x in mm, i in A]' % (fit['K'], fit['a']))
r = abs(Glf(1j * w) / Gfit(1j * w))
f1 = f[np.argmax(20 * np.log10(np.maximum(r, 1 / r)) > 1)]
print('G_lf within 1 dB of the full model up to %.1f Hz (= %.2f f_r)' % (f1, f1 / fit['f_r']))
fig, ax = plt.subplots(2, 1, figsize=(8, 5.5), sharex=True)
ax[0].semilogx(f[k], 20 * np.log10(abs(Gfit(1j * w[k]))), 'k', label='full fit')
ax[0].semilogx(f[k], 20 * np.log10(abs(Glf(1j * w[k]))), 'C1--', label='dominant poles $K/(s(s+a))$')
ax[0].axvline(f1, color='0.5', ls=':', label='1 dB apart (%.0f Hz)' % f1)
ax[0].set(ylabel='|G| (dB, mm/A)', title='Model approximation: keep the dominant (rigid-body) poles')
ax[1].semilogx(f[k], np.degrees(np.unwrap(np.angle(Gfit(1j * w[k])))), 'k', label='full fit (no delay)')
ax[1].semilogx(f[k], np.degrees(np.unwrap(np.angle(Glf(1j * w[k])))), 'C1--', label='dominant poles')
ax[1].set(xlabel='frequency (Hz)', ylabel='phase (deg)', xlim=(1, 300))
for a in ax:
    a.grid(True, which='both', alpha=0.4)
    a.legend(fontsize=8)
fig.tight_layout()
plt.show()
