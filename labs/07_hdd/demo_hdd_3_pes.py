# Author: Xu Chen, University of Washington
# %% [markdown]
# # Disk drive 3: the position error signal, runout, and what the sensitivity function does to it
#
# ![scene](figures/scene.png)
#
# **Overview.** How well a drive follows its tracks is judged by the position error signal (PES) logged
# during track following, summarized as $3\sigma$ in % of the track pitch. The PES has a repeatable part
# (RRO), locked to the disk rotation, and a non-repeatable part (NRRO): disk flutter, airflow on the arm,
# measurement noise. The loop shapes both through the sensitivity function $S = 1/(1 + L)$ from runout to
# PES (L26): runout below the crossover is followed, runout above it passes through or is amplified where
# $|S| > 1$. Pushing $|S|$ down in one band pushes it up in another (the waterbed effect). We analyse a PES
# log recorded with the firmware controller, predict the PES with the controller of demo 2, and add a peak
# filter at the spindle frequency for the largest runout component.

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
from ctrlprac import data
import _model as M

dp, pm = data.load('07_hdd_pes')
Ts, ns = pm['sample_time_s'], pm['servo_sectors']
f_rot = pm['spindle_rpm'] / 60.0
pes = dp['pes_pct'].values
fw = pm['firmware_controller']
C0, _ = M.controller(Ts, fw['K'], fw['fz'], fw['fp'], fw['fi'], tuple(fw['notch']))
print('%d revolutions at %.0f rpm (spindle frequency %.0f Hz), %d sectors per revolution, Ts = %.2f us'
      % (dp['rev'].max() + 1, pm['spindle_rpm'], f_rot, ns, 1e6 * Ts))
print('firmware controller C0: K = %.4g A/track, lead %.0f/%.0f Hz, integral %.0f Hz, notch %.0f Hz'
      % (fw['K'], fw['fz'], fw['fp'], fw['fi'], fw['notch'][0]))

# plant model and controller of demos 1 and 2
f, G, coh, _ = M.load_frf('07_hdd_frf')
idf = M.identify(f, G, coh, Ts)
P_fit = M.plant_dt(M.fitted_drive(idf['Kp'], idf['Td'], idf['modes'], Ts))
notch = (round(idf['modes'][1]['f']), 20.0, 0.2)
d1 = M.design_lead(idf['Kp'], idf['Td'], 1200.0, 40.0, Ts, fi=120.0, notch=notch)
C1, _ = M.controller(Ts, d1['K'], d1['fz'], d1['fp'], 120.0, notch)
S_of = lambda C: (lambda ff: 1 / (1 + M.frf(C, ff) * M.frf(P_fit, ff)))
S0, S1 = S_of(C0), S_of(C1)
m0 = M.frf_margins(f, M.frf(C0, f) * G)
print('C0 on the measured response: crossover %.0f Hz, PM %.1f deg, GM %.1f dB, Ms %.2f'
      % (m0['fc'][0], m0['PM'][0], 20 * np.log10(m0['GM_up'][1]), m0['Ms']))

# %% [markdown]
# ## RRO and NRRO by averaging over revolutions
# The RRO repeats every revolution, so the average of the PES at each sector over the $N$ logged revolutions
# estimates it; what is left is the NRRO. The two are uncorrelated, so
# $\sigma^2_{PES} \approx \sigma^2_{RRO} + \sigma^2_{NRRO}$. The average still contains NRRO divided by $\sqrt{N}$,
# which is why logs are long.

# %%
st0 = M.pes_stats(pes, ns)
print('PES 3 sigma = %.1f %% TP: RRO %.1f %%, NRRO %.1f %% (RSS %.1f %%); %.1f nm 3 sigma at %.0f nm track pitch'
      % (st0['three_sigma'], st0['rro3'], st0['nrro3'], np.hypot(st0['rro3'], st0['nrro3']),
         st0['three_sigma'] * pm['track_pitch_nm'] / 100, pm['track_pitch_nm']))
print('NRRO left in the RRO estimate: about %.1f %% 3 sigma (NRRO/sqrt(%d))' % (st0['nrro3'] / np.sqrt(st0['n_rev']), st0['n_rev']))
print('mean current %.1f mA (the integral term holds the flex-cable bias); mean PES %.2f %% TP'
      % (dp['u_mA'].mean(), pes.mean()))

fig, ax = plt.subplots(2, 1, figsize=(9, 6))
k5 = 5 * ns
t_ms = 1e3 * dp['time_s'].values[:k5]
ax[0].plot(t_ms, pes[:k5], 'C0', lw=0.6, label='PES')
ax[0].plot(t_ms, np.tile(st0['rro'], 5), 'C3', lw=1.4, label='RRO (average over %d revolutions)' % st0['n_rev'])
ax[0].set(xlabel='time (ms)', ylabel='PES (% of track pitch)',
          title='PES with the firmware controller, first 5 revolutions (simulated log)')
sec = np.arange(ns)
nr = st0['nrro'].reshape(-1, ns)
ax[1].plot(sec, nr.T, color='0.8', lw=0.3)
ax[1].plot(sec, st0['rro'], 'C3', lw=1.4, label='RRO: %.1f %% 3 sigma' % st0['rro3'])
ax[1].fill_between(sec, -st0['nrro3'] / 1, st0['nrro3'], color='C0', alpha=0.15, label='NRRO: +/- %.1f %% (3 sigma)' % st0['nrro3'])
ax[1].set(xlabel='servo sector', ylabel='% of track pitch', xlim=(0, ns - 1), title='One revolution: RRO and the spread of the NRRO (grey: every revolution)')
for a in ax:
    a.grid(True, alpha=0.4)
    a.legend(fontsize=8, loc='upper right')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Which runout is followed and which is amplified (L26)
# The RRO in the PES is the track runout seen through the loop, $\mathrm{RRO}_{PES}(k f_{rot}) =
# S(e^{j\omega T_s})\,r(k f_{rot})$. Dividing each harmonic of the PES by $|S_0|$ at that frequency gives the runout
# of the track itself. The NRRO spectrum divided by $|S_0|^2$ estimates the spectrum of the non-repeatable
# disturbance in the same way.

# %%
h = M.rro_harmonics(st0['rro'])                       # % TP, harmonics 1x, 2x, ...
fk = f_rot * np.arange(1, len(h) + 1)
run = h / abs(S0(fk))
true_rro = pm['disturbances']['rro_track']
print('harmonic   in PES   |S0|    track runout (estimate, true)')
for k in range(5):
    print('  %dx %4.0f Hz %5.2f %%  %.3f   %6.1f %% TP (%5.1f %%)' % (k + 1, fk[k], h[k], abs(S0(fk[k:k + 1]))[0], run[k], 100 * true_rro[k]))
amp = np.flatnonzero(abs(S0(fk)) > 1) + 1               # harmonic numbers with |S0| > 1
runs = np.split(amp, np.flatnonzero(np.diff(amp) > 1) + 1)
print('|S0| > 1 from %.0f Hz on: harmonics %s are amplified, not attenuated'
      % (amp[0] * f_rot, ' and '.join('%dx-%dx' % (r[0], r[-1]) for r in runs)))

fn, Pn = signal.welch(st0['nrro'], 1 / Ts, nperseg=10 * ns)
S0n = abs(S0(fn[1:]))
fig, ax = plt.subplots(3, 1, figsize=(8.5, 8.5), sharex=True, gridspec_kw=dict(height_ratios=[1.3, 2, 2]))
ax[0].semilogx(fn[1:], 20 * np.log10(S0n), 'k', lw=1)
ax[0].axhline(0, color='k', lw=0.6, ls=':')
ax[0].set(ylabel='|S0| (dB)', ylim=(-60, 10), title='Firmware loop: sensitivity, RRO harmonics and NRRO spectrum')
ax[1].semilogx(fk, run, 'C1s', ms=3.5, mfc='none', label='track runout (PES harmonic / |S0|)')
ax[1].vlines(fk, 1e-4, h, color='C3', lw=1.2)
ax[1].semilogx(fk, h, 'C3o', ms=3, label='RRO harmonics in the PES')
ax[1].set_yscale('log')
ax[1].set(ylabel='amplitude (% TP)', ylim=(1e-2, 300))
ax[2].loglog(fn[1:], Pn[1:], 'C0', lw=0.9, label='NRRO PSD')
ax[2].loglog(fn[1:], Pn[1:] / S0n ** 2, 'C2', lw=0.9, label='NRRO PSD / |S0|^2 (disturbance estimate)')
for fl, _, _ in pm['disturbances']['flutter']:
    ax[2].axvline(fl, color='0.5', ls=':', lw=0.8)
ax[2].text(1160, 0.12, 'disk flutter', fontsize=8, ha='center', color='0.3')
ax[2].set(xlabel='frequency (Hz)', ylabel='PSD ((% TP)^2/Hz)', xlim=(f_rot * 0.8, 18000), ylim=(1e-5, 1))
for a in ax:
    a.grid(True, which='both', alpha=0.4)
for a in ax[1:]:
    a.legend(fontsize=8, loc='upper right')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## The waterbed effect
# For a stable loop with no unstable open-loop poles the area under $\ln|S|$ over the frequency axis is zero
# (Bode's sensitivity integral; for a sampled loop over $0 \le f \le f_s/2$):
#
# $$ \int_0^{f_s/2} \ln\left|S(e^{j2\pi f T_s})\right|\,df = 0 . $$
#
# Attenuation below the crossover has to be paid for with amplification above it. Raising the crossover from
# the firmware's 0.8 kHz to 1.2 kHz deepens the attenuation at the low spindle harmonics and moves the
# amplification to higher frequencies, where some written-in runout and the windage-driven arm motion sit.

# %%
fl = np.linspace(0.01, 0.5 / Ts, 400000)
for tag, S in (('C0 (firmware)', S0), ('C1 (demo 2)', S1)):
    lnS = np.log(abs(S(fl)))
    df_ = fl[1] - fl[0]
    print('%-14s integral of ln|S| over 0..fs/2: %+.1f Np Hz (below 0 dB %.0f, above %+.0f), max |S| %.2f'
          % (tag, lnS.sum() * df_, lnS[lnS < 0].sum() * df_, lnS[lnS > 0].sum() * df_, np.exp(lnS.max())))
fig, ax = plt.subplots(figsize=(9, 4))
for S, col, lab in ((S0, 'C7', 'C0 firmware, crossover 0.8 kHz'), (S1, 'C0', 'C1 demo 2, crossover 1.2 kHz')):
    v = 20 * np.log10(abs(S(fl)))
    ax.plot(fl / 1e3, v, col, lw=1.2, label=lab)
    ax.fill_between(fl / 1e3, v, 0, where=v < 0, color=col, alpha=0.12)
    ax.fill_between(fl / 1e3, v, 0, where=v > 0, color=col, alpha=0.35)
ax.axhline(0, color='k', lw=0.6)
ax.set(xlabel='frequency (kHz, linear axis)', ylabel='|S| (dB)', xlim=(0, 18), ylim=(-30, 8),
       title='Waterbed: equal areas of ln|S| below and above 0 dB')
ax.grid(True, alpha=0.4)
ax.legend(fontsize=8, loc='lower right')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Predicting the PES with the new controller
# For a linear loop with the same disturbances, the PES with controller $C_1$ is the logged PES filtered by
# $S_1/S_0$ (`M.predict_pes`, applied on the FFT of the whole-revolution log). This uses only the log and
# the identified plant. As a check we also simulate the drive with the disturbance model of the demo
# sidecar and the same random sequences (`M.simulate_pes`), which the homework drive does not provide. The
# simulation runs the loop once over the sequence before the logged pass, so that the slow closed-loop poles
# that the peak filter adds near 120 Hz have settled.

# %%
peak = (f_rot, 20.0, 0.1)
C1p, _ = M.controller(Ts, d1['K'], d1['fz'], d1['fp'], 120.0, notch, peak)
m1p = M.frf_margins(f, M.frf(C1p, f) * G)
m1 = M.frf_margins(f, M.frf(C1, f) * G)
print('C1          on the measured response: PM %.1f deg, GM %.1f dB, Ms %.2f' % (m1['PM'][0], 20 * np.log10(m1['GM_up'][1]), m1['Ms']))
print('C1 + peak   on the measured response: PM %.1f deg, GM %.1f dB, Ms %.2f (peak filter %.0f Hz, +20 dB, zeta_z 0.1)'
      % (m1p['PM'][0], 20 * np.log10(m1p['GM_up'][1]), m1p['Ms'], peak[0]))
print('|S| at %.0f Hz: C0 %.4f, C1 %.4f, C1 + peak %.5f' % (f_rot, abs(S0([f_rot]))[0], abs(S1([f_rot]))[0],
                                                          abs(S_of(C1p)([f_rot]))[0]))
rows = [('C0 firmware (log)', pes, None)]
for tag, C in (('C1', C1), ('C1 + peak filter', C1p)):
    pr = M.predict_pes(pes, Ts, S0, S_of(C))
    sim = 100 * M.simulate_pes(M.DRIVE_A, C, st0['n_rev'], pm['seed'])['pes']
    rows += [(tag + ' (predicted)', pr, None), (tag + ' (simulated)', sim, None)]
print('%-30s %8s %8s %8s %10s' % ('PES, % of track pitch', '3 sigma', 'RRO', 'NRRO', '1x RRO amp'))
stats = {}
for tag, x, _ in rows:
    s = M.pes_stats(x, ns)
    stats[tag] = (s, x)
    print('%-30s %8.1f %8.1f %8.1f %10.2f' % (tag, s['three_sigma'], s['rro3'], s['nrro3'], M.rro_harmonics(s['rro'])[0]))

# %% [markdown]
# ## Where the PES variance comes from
# The cumulative $3\sigma$ adds up the PES spectrum from low to high frequency; its final value is the $3\sigma$
# of the whole log. Steps mark single components (the spindle harmonics), ramps mark broadband parts.

# %%
def cum3s(x):
    X = np.fft.rfft(x - x.mean())
    p = np.abs(X) ** 2
    p[1:-1] *= 2
    return np.fft.rfftfreq(len(x), Ts), 3 * np.sqrt(np.cumsum(p) / len(x) ** 2)


fig, ax = plt.subplots(1, 2, figsize=(11, 4.3))
for tag, col in (('C0 firmware (log)', 'C7'), ('C1 (predicted)', 'C0'), ('C1 + peak filter (predicted)', 'C3')):
    ff, c = cum3s(stats[tag][1])
    ax[0].semilogx(ff[1:], c[1:], col, lw=1.3, label='%s: %.1f %%' % (tag, stats[tag][0]['three_sigma']))
ax[0].set(xlabel='frequency (Hz)', ylabel='cumulative 3 sigma (% TP)', xlim=(50, 18000),
          title='Cumulative PES 3 sigma')
ax[0].legend(fontsize=7.5, loc='upper left')
fs_ = np.logspace(np.log10(20), np.log10(0.4999 / Ts), 5000)
for S, col, lab in ((S0, 'C7', 'C0'), (S1, 'C0', 'C1'), (S_of(C1p), 'C3', 'C1 + peak filter')):
    ax[1].semilogx(fs_, 20 * np.log10(abs(S(fs_))), col, lw=1.2, label=lab)
ax[1].axhline(0, color='k', lw=0.6)
ax[1].set(xlabel='frequency (Hz)', ylabel='|S| (dB)', xlim=(20, 18000), ylim=(-70, 10), title='Sensitivity functions')
ax[1].legend(fontsize=8, loc='lower right')
for a in ax:
    a.grid(True, which='both', alpha=0.4)
fig.tight_layout()
plt.show()

# %% [markdown]
# **Reading the result.** With the firmware controller the PES is dominated by the 1x runout, which the loop
# attenuates only by a factor of about 15 at 120 Hz. The 1.2 kHz design of demo 2 attenuates it more than twice as
# much and cuts the PES $3\sigma$ roughly in half, while the waterbed shows up as more amplification between 1 and
# 3 kHz. The peak filter adds 20 dB of loop gain in a narrow band at 120 Hz and removes most of the remaining 1x
# runout at almost no cost in margins. The prediction from the log and the simulation agree to a few tenths
# of a percent, so the log alone was enough to evaluate the new controller.
