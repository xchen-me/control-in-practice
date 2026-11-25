# Author: Xu Chen, University of Washington
# %% [markdown]
# # Disk drive 1: identifying the VCM actuator from a measured frequency response
#
# ![scene](figures/scene.png)
#
# **Overview.** A hard disk drive keeps its read/write head on a data track, and the tracks are only 60 nm
# apart (the track pitch). The head sits at the end of an arm that a voice-coil motor (VCM) rotates about a
# pivot. At low frequency the arm is a rigid inertia, a double integrator from coil current to head position;
# at a few kHz the arm, the pivot bearing and the suspension add lightly damped resonances. The servo
# samples the head position once per servo sector (36 kHz here), so the zero-order hold and the computation
# delay also count.
# The drive's servo analyzer measured the frequency response from VCM current to head position with a
# swept sine (L20). From that measurement we read the delay from the phase slope, the rigid-body gain
# from the -40 dB/dec line, and every resonance from its peak (L21), and then refine all of them with a
# least-squares fit of the modal model.

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
import _model as M

# %% [markdown]
# ## The measurement
# The analyzer adds a sine to the current command while the firmware controller keeps the head on the
# track, steps it through 400 frequencies from 100 Hz to 17 kHz, and estimates
# $\hat G(j\omega) = H_{exc \to y}/H_{exc \to u}$ with the excitation as reference (as in Lab 06). It exports
# the magnitude in dB re 1 track/A, the phase wrapped to $\pm 180^\circ$ and the coherence $\gamma^2$.

# %%
f, G, coh, meta = M.load_frf('07_hdd_frf')
Ts = meta['sample_time_s']
print('drive %s: %d rpm, %d servo sectors -> Ts = %.2f us (%.1f kHz), track pitch %.0f nm'
      % (meta['drive'], meta['spindle_rpm'], meta['servo_sectors'], 1e6 * Ts, 1e-3 / Ts, meta['track_pitch_nm']))
ph = np.degrees(np.unwrap(np.angle(G)))
ph -= 360 * np.round((ph[0] + 180) / 360)               # rigid body: start near -180 deg
good = coh >= 0.9
print('coherence >= 0.9 on %d of %d points; below 0.9 near %s Hz'
      % (good.sum(), len(f), ', '.join('%.0f' % x for x in f[~good][np.r_[True, np.diff(f[~good]) > 200]])))

fig, ax = plt.subplots(3, 1, figsize=(8.5, 8), sharex=True, gridspec_kw=dict(height_ratios=[3, 3, 1.4]))
ax[0].semilogx(f, 20 * np.log10(abs(G)), color='0.75', lw=0.8)
ax[0].semilogx(f[good], 20 * np.log10(abs(G[good])), 'C0.', ms=2.5, label='coherence >= 0.9')
ax[0].set(ylabel='|G| (dB re 1 track/A)', title='Measured frequency response, VCM current -> head position')
ax[1].semilogx(f, ph, color='0.75', lw=0.8)
ax[1].semilogx(f[good], ph[good], 'C0.', ms=2.5)
ax[1].set(ylabel='phase (deg, unwrapped)')
ax[2].semilogx(f, coh, 'C0', lw=0.8)
ax[2].axhline(0.9, color='k', ls=':', lw=0.8)
ax[2].set(xlabel='frequency (Hz)', ylabel='coherence', ylim=(0, 1.05), xlim=(100, 17000))
for a in ax:
    a.grid(True, which='both', alpha=0.4)
ax[0].legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# The coherence is close to 1 almost everywhere. It dips where the response is small compared with the
# disturbances: at the disk-flutter frequencies near 0.9 and 1.5 kHz (non-repeatable runout that the
# excitation competes with), at the antiresonance near 5 kHz, and above 12 kHz, where the plant gain has
# fallen by another 20 dB.
#
# ## Delay from the phase slope (L20)
# Below the resonances the plant is $K_p/s^2$, phase $-180^\circ$. The zero-order hold acts like a delay
# of $T_s/2$ and the computation adds $\tau$, so the measured phase is a straight line in linear frequency:
#
# $$ \angle G(j\omega) \approx -180^\circ - 360^\circ f\,T_d, \qquad T_d = \frac{T_s}{2} + \tau . $$
#
# A free straight-line fit also estimates the intercept; since physics fixes it at $-180^\circ$, the line
# through $-180^\circ$ uses the data only for the slope and is less sensitive to the phase noise near the
# disk-flutter frequencies.

# %%
lo = good & (coh >= 0.95) & (f >= 150) & (f <= 1500)
slope, icpt = np.polyfit(f[lo], ph[lo], 1)
Td = -np.sum(f[lo] * (ph[lo] + 180)) / (360 * np.sum(f[lo] ** 2))     # line through -180 deg at f = 0
td = meta['true_derived']
print('phase fit on %.0f-%.0f Hz: free line %.1f deg at f = 0 (ideal -180), slope -> Td = %.2f us'
      % (f[lo].min(), f[lo].max(), icpt, -1e6 * slope / 360))
print('line through -180 deg: Td = %.2f us (true %.2f us)' % (1e6 * Td, td['Td_us']))
print('-> computation delay tau = Td - Ts/2 = %.2f us; phase lag of the delay at 1 kHz: %.1f deg, at 6 kHz: %.0f deg'
      % (1e6 * (Td - Ts / 2), 360 * 1e3 * Td, 360 * 6e3 * Td))

# %% [markdown]
# ## Rigid-body gain from the -40 dB/dec line (L21)
# Below the modes $|G| \approx K_p/\omega^2$, so $|G|\,\omega^2$ is flat and equal to $K_p$. The modes pull it down as
# the frequency rises (by 0.8 % at 600 Hz and about 6 % at 1.5 kHz), so $K_p$ is the median of $|G|\,\omega^2$ on
# 150-600 Hz. With the VCM torque constant $K_t$, the pivot-to-head length $l$ and the track pitch TP,
# $K_p = K_t l/(J\,\mathrm{TP})$ gives the actuator inertia $J$.

# %%
lo2 = lo & (f <= 600)
Kp = np.median(abs(G[lo2]) * (2 * np.pi * f[lo2]) ** 2)
p = meta['parameters']
J = p['Kt'] * p['l_arm'] / (Kp * meta['track_pitch_nm'] * 1e-9)
print('Kp = %.4g track/(A s^2) (true %.4g, %+.1f %%); J = Kt l/(Kp TP) = %.3g kg m^2 (true %.3g)'
      % (Kp, td['Kp_track_per_As2'], 100 * (Kp / td['Kp_track_per_As2'] - 1), J, p['J']))
print('acceleration: %.0f m/s^2 per A (%.0f g/A); 1 mA at 1 kHz moves the head %.2f track'
      % (Kp * meta['track_pitch_nm'] * 1e-9, Kp * meta['track_pitch_nm'] * 1e-9 / 9.81, 1e-3 * Kp / (2 * np.pi * 1e3) ** 2))

# %% [markdown]
# ## Resonances from the peaks (L21)
# Dividing the measurement by the rigid-body line $-K_p/\omega^2\,e^{-j\omega T_d}$ leaves the effect of
# the modes. In the modal model
#
# $$ G(s) = K_p\left[\frac{1}{s^2} + \sum_i \frac{r_i}{s^2 + 2\zeta_i\omega_i s + \omega_i^2}\right] $$
#
# mode $i$ peaks near $\omega_i$, and its half-power width gives $\zeta_i \approx \Delta f/(2 f_i)$ (for a peak well
# above the line; the small pivot-bearing peak has no clear half-power points, so we start from
# $\zeta = 0.05$). Each mode also makes an antiresonance where the rigid-body and modal terms cancel,
# $\omega_{z,i} = \omega_i/\sqrt{1 + r_i}$, so $r_i = (f_i/f_{z,i})^2 - 1$: an in-phase mode ($r_i > 0$) has its
# antiresonance below the resonance, an out-of-phase mode ($r_i < 0$) above it.

# %%
R = G / (-Kp / (2 * np.pi * f) ** 2 * np.exp(-2j * np.pi * f * Td))
RdB = 20 * np.log10(abs(R))
band = np.flatnonzero((f > 2000) & (f < 12000))            # peaks standing 3 dB above their surroundings
pk = band[signal.find_peaks(RdB[band], prominence=3)[0]]
pk = pk[coh[pk] >= 0.9]
names = ['pivot-bearing mode', 'butterfly mode', 'suspension sway mode']
print('peaks found at', np.round(f[pk]))
modes0 = []
for i, k in enumerate(pk):
    half = RdB[k] - 3
    k1, k2 = k, k
    while k1 > 0 and RdB[k1] > half:
        k1 -= 1
    while k2 < len(f) - 1 and RdB[k2] > half:
        k2 += 1
    zeta = (f[k2] - f[k1]) / (2 * f[k]) if RdB[k] > 10 else 0.05
    near = np.flatnonzero((f > 0.75 * f[k]) & (f < 1.3 * f[k]) & (coh >= 0.3))
    kz = near[np.argmin(RdB[near])]
    r = (f[k] / f[kz]) ** 2 - 1
    modes0.append(dict(name=names[i], f=f[k], zeta=zeta, r=r))
    print('%-21s peak %6.0f Hz, %4.1f dB above the line, zeta %.3f %s; antiresonance %6.0f Hz -> r = %+.2f'
          % (names[i], f[k], RdB[k], zeta, '(half power)' if RdB[k] > 10 else '(guess)  ', f[kz], r))

fig, ax = plt.subplots(2, 1, figsize=(8.5, 6), sharex=True)
ax[0].semilogx(f[good], RdB[good], 'C0', lw=1)
ax[0].semilogx(f[pk], RdB[pk], 'r^', ms=7, label='peaks')
for m, (dy, ha) in zip(modes0, ((6, 'center'), (4, 'right'), (5, 'center'))):
    k = np.argmin(abs(f - m['f']))
    ax[0].text(m['f'] * (0.93 if ha == 'right' else 1.0), RdB[k] + dy, '%s\n%.0f Hz' % (m['name'], m['f']),
               fontsize=8, ha=ha, va='bottom')
ax[0].axhline(0, color='k', lw=0.6)
ax[0].set(ylabel='|G| / rigid line (dB)', ylim=(-30, 42), title='Measured response divided by the rigid-body line')
ax[1].semilogx(f[good], np.degrees(np.angle(R[good])), 'C0', lw=1)
ax[1].set(xlabel='frequency (Hz)', ylabel='phase (deg)', xlim=(100, 17000), ylim=(-200, 200))
ax[1].set_yticks([-180, -90, 0, 90, 180])
for a in ax:
    a.grid(True, which='both', alpha=0.4)
ax[0].legend(fontsize=8, loc='upper left')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Least-squares fit of the modal model
# The readings above are starting values. `M.fit_modal` adjusts $K_p$, $T_d$ and every $(f_i, \zeta_i, r_i)$
# together so that the sampled model (exact zero-order hold and delay) matches the measured points with
# coherence $\ge 0.8$, weighting the log-magnitude and phase errors by the coherence.

# %%
fit = M.fit_modal(f, G, coh, Ts, Kp, Td, modes0)
true_modes = meta['parameters']['modes']
print('fit on %d points: rms error %.2f dB, %.2f deg' % (fit['n_points'], fit['rms_dB'], fit['rms_deg']))
print('Kp = %.4g track/(A s^2) (true %.4g), Td = %.2f us (true %.2f)'
      % (fit['Kp'], td['Kp_track_per_As2'], 1e6 * fit['Td'], td['Td_us']))
print('%-21s %16s %16s %16s' % ('mode', 'f (Hz)', 'zeta', 'r'))
for m, mt in zip(fit['modes'], true_modes):
    print('%-21s %7.0f (%6.0f) %7.4f (%6.4f) %+7.3f (%+6.3f)' % (m['name'], m['f'], mt['f'], m['zeta'], mt['zeta'],
                                                                m['r'], mt['r']))
print('(true values in parentheses, from the demo sidecar)')
Gfit = M.frf(M.plant_dt(M.fitted_drive(fit['Kp'], fit['Td'], fit['modes'], Ts)), f)

fig, ax = plt.subplots(2, 1, figsize=(8.5, 6.5), sharex=True)
ax[0].semilogx(f, 20 * np.log10(abs(G)), 'C0.', ms=3, alpha=0.6, label='measured')
ax[0].semilogx(f, 20 * np.log10(abs(Gfit)), 'k', lw=1, label='fitted modal model (sampled, with delay)')
ax[0].semilogx(f, 20 * np.log10(fit['Kp'] / (2 * np.pi * f) ** 2), 'C1--', lw=1, label='rigid body $K_p/\\omega^2$')
ax[0].set(ylabel='|G| (dB re 1 track/A)', title='Modal model fitted to the measured response')
phf = np.degrees(np.unwrap(np.angle(Gfit)))
phf -= 360 * np.round((phf[0] + 180) / 360)
ax[1].semilogx(f, ph, 'C0.', ms=3, alpha=0.6, label='measured')
ax[1].semilogx(f, phf, 'k', lw=1, label='fit')
ax[1].semilogx(f, -180 - 360 * f * fit['Td'], 'C1--', lw=1, label='rigid body with delay $T_d$')
ax[1].set(xlabel='frequency (Hz)', ylabel='phase (deg)', xlim=(100, 17000))
for a in ax:
    a.grid(True, which='both', alpha=0.4)
    a.legend(fontsize=8, loc='lower left')
fig.tight_layout()
plt.show()

# %% [markdown]
# **Reading the result.** Three numbers describe the low-frequency plant: $K_p$, the delay $T_d$ (zero-order
# hold plus computation, about one sample in total) and the sample rate. The butterfly mode stands about
# 25 dB above the rigid-body line, and at 6 kHz the delay alone costs about 55 degrees of phase, so the
# controller can not rely on the phase of this mode. Demo 2 designs the track-following controller on the rigid-body
# model and then checks it against this measured response, with and without a notch at the butterfly mode.
