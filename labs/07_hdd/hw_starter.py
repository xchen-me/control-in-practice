# Author: Xu Chen, University of Washington
# %% [markdown]
# # Homework starter, Lab 07: disk drive D2
#
# ![scene](figures/scene_hw.png)
#
# Loads the two D2 logs, builds the firmware controller from its `.json` and plots the raw data. Fill in the
# TODO steps (see homework.md).

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import control as ct
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
from ctrlprac import data
import _model as M          # M.controller, M.design_lead, M.frf, M.frf_margins, M.fit_modal, M.pes_stats, ...

f, G, coh, meta = M.load_frf('07_hw_frf')
dp, pmeta = data.load('07_hw_pes')
Ts, ns = meta['sample_time_s'], meta['servo_sectors']
f_rot = meta['spindle_rpm'] / 60.0
print('drive %s: %d rpm (%.0f Hz), %d sectors, Ts = %.2f us, track pitch %.0f nm; Kt = %.3f N m/A, l_arm = %.3f m'
      % (meta['drive'], meta['spindle_rpm'], f_rot, ns, 1e6 * Ts, meta['track_pitch_nm'],
         meta['parameters']['Kt'], meta['parameters']['l_arm']))
fw = pmeta['firmware_controller']
C0, _ = M.controller(Ts, fw['K'], fw['fz'], fw['fp'], fw['fi'], tuple(fw['notch']))
print('firmware controller C0:', {k: fw[k] for k in ('K', 'fz', 'fp', 'fi', 'notch')})

fig, ax = plt.subplots(3, 1, figsize=(8.5, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[3, 2, 1.3]))
ax[0].semilogx(f, 20 * np.log10(abs(G)), 'C0.', ms=3)
ax[0].set(ylabel='|G| (dB re 1 track/A)', title='07_hw_frf (analyzer export)')
ax[1].semilogx(f, np.degrees(np.angle(G)), 'C0.', ms=3)
ax[1].set(ylabel='phase (deg, wrapped)')
ax[2].semilogx(f, coh, 'C0', lw=0.8)
ax[2].set(xlabel='frequency (Hz)', ylabel='coherence', ylim=(0, 1.05))
for a in ax:
    a.grid(True, which='both', alpha=0.4)
fig.tight_layout()
plt.show()

fig, ax = plt.subplots(2, 1, figsize=(8.5, 5.5))
k3 = 3 * ns
ax[0].plot(1e3 * dp['time_s'][:k3], dp['pes_pct'][:k3], 'C0', lw=0.6)
ax[0].set(xlabel='time (ms)', ylabel='PES (% TP)', title='07_hw_pes, first 3 revolutions')
ax[1].plot(1e3 * dp['time_s'][:k3], dp['u_mA'][:k3], 'C1', lw=0.6)
ax[1].set(xlabel='time (ms)', ylabel='VCM current (mA)')
for a in ax:
    a.grid(True, alpha=0.4)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 1. Frequency response and coherence
# TODO: unwrap the phase (np.unwrap), mark the points you trust, explain the coherence dips

# %% [markdown]
# ## 2. Delay and rigid-body gain
# TODO: Td from the phase slope (line through -180 deg), tau = Td - Ts/2; Kp from |G| w^2; J = Kt l_arm/(Kp TP)

# %% [markdown]
# ## 3. Resonances and the modal fit
# TODO: divide G by -Kp/w^2 e^{-jwTd}; peaks, half-power zeta, r from the antiresonance;
# fit = M.fit_modal(f, G, coh, Ts, Kp, Td, modes0); P_fit = M.plant_dt(M.fitted_drive(fit['Kp'], fit['Td'], fit['modes'], Ts))

# %% [markdown]
# ## 4-5. Lead + integral, then the notch and the highest crossover
# TODO: d = M.design_lead(Kp, Td, fc, 40, Ts, fi=fc/10, notch=...); C, _ = M.controller(Ts, d['K'], d['fz'], d['fp'], fc/10, notch)
# a single lead stage supplies less than 90 deg: check d['phi'] (M.design_lead warns when it would need more)
# L = M.frf(C, f) * G; m = M.frf_margins(f, L); stability: abs(M.cl_poles(P_fit, C)).max() < 1
# check: ct.stability_margins(ct.frd(L, 2 * np.pi * f))

# %% [markdown]
# ## 6. PES: RRO, NRRO, runout through |S0|
# TODO: st = M.pes_stats(dp['pes_pct'].values, ns); M.rro_harmonics(st['rro']); S0 = 1/(1 + C0 G)

# %% [markdown]
# ## 7. Prediction and the peak filter
# TODO: M.predict_pes(pes, Ts, S0_fun, S1_fun) with S_fun = lambda ff: 1 / (1 + M.frf(C, ff) * M.frf(P_fit, ff))
