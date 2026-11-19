# Author: Xu Chen, University of Washington
# %% [markdown]
# # Belt drive 1: measuring a frequency response and reading the Bode plot
#
# ![scene](figures/scene.png)
#
# **Overview.** A belt-driven linear stage behaves like a rigid mass at low frequency, but the elastic
# belt makes a lightly damped resonance in the tens of Hz, and that resonance, not the motor, limits how
# fast the position loop can be. Before designing anything we measure the frequency response from
# current to carriage position (L20): a PRBS test analysed with the empirical transfer function estimate
# (ETFE) and its coherence, cross-checked by a stepped-sine test. Then we read the Bode plot (L21):
# the -40 dB/dec rigid-body line, the resonance frequency and its damping from the peak height.

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
from ctrlprac import data, ident
import _model as M

# %% [markdown]
# ## The test log
# A double integrator drifts away under open-loop excitation, so the drive holds the carriage with a
# slow "parking" PD loop (crossover near 2 Hz) and adds the excitation current $i_{exc}$ to the
# controller output. The log has the excitation, the total current command and the encoder position,
# sampled at 4 kHz.

# %%
df, meta = data.load('06_drive_prbs')
t, d, u, x = (df[k].values for k in ('time_s', 'i_exc_A', 'i_cmd_A', 'x_mm'))
print('fs = %.0f Hz, %.1f s, encoder resolution %.0f um, known drive delay %d samples'
      % (1 / meta['sample_time_s'], t[-1], 1e3 * meta['encoder_resolution_mm'], meta['delay_samples']))
fig, ax = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
ax[0].plot(t, u, lw=0.4, label='current command $i_{cmd}$')
ax[0].plot(t, d, lw=0.4, alpha=0.6, label='excitation $i_{exc}$ (PRBS)')
ax[0].set_ylabel('current (A)')
ax[1].plot(t, x, lw=0.6, color='C2', label='carriage position x')
ax[1].set(xlabel='time (s)', ylabel='x (mm)')
ax[0].set_title('PRBS identification test inside the parking loop (simulated log)')
for a in ax:
    a.grid(True, alpha=0.4)
    a.legend(loc='upper right', fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## ETFE and coherence (L20)
# `ident.etfe(t, u, y)` returns $\hat H(j\omega) = S_{uy}/S_{uu}$ (Welch averages) and the coherence
# $\gamma^2 = |S_{uy}|^2/(S_{uu}S_{yy})$. In a closed loop the current $i_{cmd}$ also contains the
# controller's reaction to sensor noise, so we use the excitation as the reference signal:
#
# $$ \hat G(j\omega) = \frac{\hat H_{exc \to x}(j\omega)}{\hat H_{exc \to i_{cmd}}(j\omega)} $$
#
# Coherence near 1 means the output is linearly explained by the input; trust the estimate only there.

# %%
nps = 8192                              # 2 s windows: 0.49 Hz resolution, 6 averages (50 % overlap)
w, Hdx, coh = ident.etfe(t, d, x, nperseg=nps)
_, Hdu, _ = ident.etfe(t, d, u, nperseg=nps)
G_etfe = Hdx / Hdu
_, G_direct, _ = ident.etfe(t, u, x, nperseg=nps)
f = w / 2 / np.pi
good = (coh > 0.9) & (f >= 1.5) & (f <= 300)
print('coherence > 0.9 from %.1f to %.1f Hz' % (f[good].min(), f[good].max()))

# %% [markdown]
# ## Stepped sine (L20)
# One frequency at a time. The amplitude is faded in and the transient allowed to die out (longer near
# the lightly damped resonance, whose ringing lasts about $1/(\zeta\omega_r)$); in the analysis window
# (`measure = 1`, an integer number of periods) `ident.sine_response` projects the signals on sin/cos.
# Again the ratio of the responses to the excitation gives the plant.

# %%
ss, smeta = data.load('06_drive_steppedsine')
fs_list, G_ss = [], []
for seg, g in ss.groupby('segment'):
    fq = g['f_exc_Hz'].iloc[0]
    g = g[g['measure'] == 1]
    tt = g['time_s'].values
    Ax, phx = ident.sine_response(tt, g['i_exc_A'].values, g['x_mm'].values, 2 * np.pi * fq)
    Au, phu = ident.sine_response(tt, g['i_exc_A'].values, g['i_cmd_A'].values, 2 * np.pi * fq)
    fs_list.append(fq)
    G_ss.append(Ax / Au * np.exp(1j * np.radians(phx - phu)))
fs_arr, G_ss = np.array(fs_list), np.array(G_ss)


def phase_deg(H):
    p = np.degrees(np.unwrap(np.angle(H)))
    return p - 360 * np.round((p[0] + 90) / 360)        # start near -90 deg (integrator)


fig, ax = plt.subplots(3, 1, figsize=(8.5, 8), sharex=True, gridspec_kw=dict(height_ratios=[3, 3, 1.4]))
ax[0].semilogx(f[f > 1], 20 * np.log10(abs(G_direct[f > 1])), color='0.7', lw=0.8, label='direct $S_{ux}/S_{uu}$')
ax[0].semilogx(f[f > 1], 20 * np.log10(abs(G_etfe[f > 1])), 'C0', lw=1, label='ETFE via excitation')
ax[0].semilogx(fs_arr, 20 * np.log10(abs(G_ss)), 'C3o', ms=5, mfc='none', label='stepped sine')
ax[0].set(ylabel='|G| (dB, mm/A)', title='Measured frequency response, current -> carriage position')
ax[1].semilogx(f[good], phase_deg(G_etfe[good]), 'C0', lw=1, label='ETFE (coherence > 0.9)')
ax[1].semilogx(fs_arr, phase_deg(G_ss), 'C3o', ms=5, mfc='none', label='stepped sine')
ax[1].set(ylabel='phase (deg)', ylim=(-450, 0))
ax[1].set_yticks([-450, -360, -270, -180, -90, 0])
ax[2].semilogx(f[f > 1], coh[f > 1], 'C0', lw=0.8)
ax[2].axhline(0.9, color='k', ls=':', lw=0.8)
ax[2].set(xlabel='frequency (Hz)', ylabel='coherence', ylim=(0, 1.05), xlim=(1, 1000))
for a in ax:
    a.grid(True, which='both', alpha=0.4)
ax[0].legend(fontsize=8)
ax[1].legend(fontsize=8)
fig.tight_layout()
plt.show()
for fq in (5, 20, 60, 100):
    k = np.argmin(abs(f - fq))
    print('  %3d Hz: ETFE %.4f mm/A, %.0f deg, coherence %.3f' % (fq, abs(G_etfe[k]),
                                                                 np.degrees(np.angle(G_etfe[k])), coh[k]))

# %% [markdown]
# ## Reading the Bode plot (L21)
# * **Rigid body.** Below the resonance $|G| \approx K/\omega^2$ (slope -40 dB/dec, phase near -180 deg):
#   fit $K$ on 4-12 Hz (above the friction corner near 1 Hz, below the resonance's quasi-static
#   rise, which is already $1/(1 - (f/f_r)^2)$ = +4 % at 12 Hz). Since $K = K_t r/J$, it gives the total inertia if $K_t$ and $r$ are known.
# * **Resonance.** A second-order mode multiplies the rigid-body line by
#   $\omega_r^2/(s^2 + 2\zeta\omega_r s + \omega_r^2)$, whose peak height at $\omega_r$ is $Q = 1/(2\zeta)$. So
#   $\zeta \approx 1/(2Q)$ with $Q = |G(j\omega_r)|/(K/\omega_r^2)$.
# * Above the resonance the slope becomes -80 dB/dec and the phase falls by another 180 deg (the encoder
#   sits on the far side of the belt); the drive delay bends the phase down further.

# %%
band = good & (f >= 4) & (f <= 12)
K_rb = np.exp(np.mean(np.log(abs(G_etfe[band]) * w[band] ** 2)))
slope = np.polyfit(np.log10(w[band]), 20 * np.log10(abs(G_etfe[band])), 1)[0]
rng = good & (f > 30) & (f < 120)
kpk = np.argmax(abs(G_etfe[rng]))
f_r, Gpk = f[rng][kpk], abs(G_etfe[rng][kpk])
Q = Gpk / (K_rb / (2 * np.pi * f_r) ** 2)
zeta_r = 1 / (2 * Q)
# finer: refine the peak with the stepped-sine points around it
near = fs_arr > 30
kss = np.flatnonzero(near)[np.argmax(abs(G_ss[near]))]
Q_ss = abs(G_ss[kss]) / (K_rb / (2 * np.pi * fs_arr[kss]) ** 2)
td = meta['true_derived']
J_est = meta['parameters']['Kt'] * meta['parameters']['r'] * 1e3 / K_rb
print('rigid-body slope %.1f dB/dec (ideal -40); K = %.0f mm/(A s^2) (true %.0f, %+.1f %%)'
      % (slope, K_rb, td['K_rigid_mm_per_As2'], 100 * (K_rb / td['K_rigid_mm_per_As2'] - 1)))
print('total inertia J = Kt r/K = %.3g kg m^2 (true %.3g)' % (J_est, td['J_total_kgm2']))
print('resonance f_r = %.1f Hz (true %.1f), peak %.1f dB above the rigid-body line: Q = %.1f, zeta = %.3f (true %.3f)'
      % (f_r, td['f_r_Hz'], 20 * np.log10(Q), Q, zeta_r, td['zeta_r']))
print('stepped sine: peak at %.0f Hz, Q = %.1f -> zeta = %.3f (2 Hz grid near the peak: a coarse grid misses the top)'
      % (fs_arr[kss], Q_ss, 1 / (2 * Q_ss)))

fig, ax = plt.subplots(figsize=(8, 4.2))
ax.loglog(f[good], abs(G_etfe[good]), 'C0', lw=1, label='measured (ETFE)')
ff = np.logspace(0, np.log10(300), 200)
ax.loglog(ff, K_rb / (2 * np.pi * ff) ** 2, 'k--', lw=1, label='rigid body $K/\\omega^2$ (-40 dB/dec)')
ax.loglog(f_r, Gpk, 'r^', ms=8, label='peak: $f_r$ = %.1f Hz, Q = %.1f, $\\zeta$ = %.3f' % (f_r, Q, zeta_r))
ax.annotate('', xy=(f_r * 1.08, Gpk), xytext=(f_r * 1.08, K_rb / (2 * np.pi * f_r) ** 2),
            arrowprops=dict(arrowstyle='<->', color='r'))
ax.text(f_r * 1.15, np.sqrt(Gpk * K_rb / (2 * np.pi * f_r) ** 2), 'Q = 1/(2$\\zeta$)', color='r', fontsize=9)
ax.loglog(ff[ff > f_r], K_rb / (2 * np.pi * f_r) ** 2 * (f_r / ff[ff > f_r]) ** 4, 'k:', lw=1,
          label='-80 dB/dec above $f_r$')
ax.set(xlabel='frequency (Hz)', ylabel='|G| (mm/A)', title='Reading the Bode magnitude')
ax.grid(True, which='both', alpha=0.4)
ax.legend(fontsize=8)
fig.tight_layout()
plt.show()
