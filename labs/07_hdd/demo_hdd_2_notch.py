# Author: Xu Chen, University of Washington
# %% [markdown]
# # Disk drive 2: lead and integral design, a notch for the butterfly mode, margins on the measured response
#
# ![scene](figures/scene.png)
#
# **Overview.** The track-following controller has to make the head follow the track runout and reject
# disturbances up to about 1 kHz, so the loop crossover should be as high as the plant allows. We design
# the controller on the rigid-body model from demo 1, $K_p e^{-sT_d}/s^2$: a lead for the phase margin at
# the crossover (L18, L21) and an integral term for the bias force of the flex cable (L12). The rigid-body
# model ignores the resonances, so the decisive check is on the MEASURED frequency response: Bode and
# Nyquist plots (L24), gain, phase and delay margins and the sensitivity peak $M_s$ (L26). Without a notch
# the butterfly mode destabilizes the loop; a notch filter at that mode restores the margins.

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import control as ct
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
import _model as M

f, G, coh, meta = M.load_frf('07_hdd_frf')
Ts = meta['sample_time_s']
idf = M.identify(f, G, coh, Ts)                       # the identification of demo 1
Kp, Td = idf['Kp'], idf['Td']
P_fit = M.plant_dt(M.fitted_drive(Kp, Td, idf['modes'], Ts))
f_bfly = idf['modes'][1]['f']
print('identified: Kp = %.4g track/(A s^2), Td = %.2f us, butterfly mode %.0f Hz (zeta %.4f)'
      % (Kp, 1e6 * Td, f_bfly, idf['modes'][1]['zeta']))
fd = np.logspace(1, np.log10(0.4999 / Ts), 20000)    # dense grid up to the Nyquist frequency (model)
Gd = M.frf(P_fit, fd)

# %% [markdown]
# ## Design targets
# A crossover of $f_c$ = 1.2 kHz, about one fifth of the butterfly frequency, with PM = 40 deg; typical
# track-following specifications ask for PM $\ge$ 35-40 deg, GM $\ge$ 6 dB and $M_s \le 2$ (6 dB). The integral
# term $(s + \omega_i)/s$ with $f_i = f_c/10$ = 120 Hz removes the steady offset from the bias force and adds
# low-frequency gain; it costs $\arctan(f_i/f_c) \approx 6^\circ$ at the crossover, which the lead makes up.
#
# ## Lead on the rigid-body model (L18, L21)
# At $\omega_c$ the rigid body with the delay has the phase $-180^\circ - \omega_c T_d$, so the lead must add
#
# $$ \phi_{max} = PM - 180^\circ - \angle\left[\frac{K_p}{(j\omega_c)^2}e^{-j\omega_c T_d}\,C_{rest}(j\omega_c)\right],\qquad
#    \alpha = \frac{1 - \sin\phi_{max}}{1 + \sin\phi_{max}},\quad f_z = f_c\sqrt\alpha,\ f_p = f_c/\sqrt\alpha , $$
#
# where $C_{rest}$ holds the integral term (and later the notch), and $K$ makes $|L(j\omega_c)| = 1$. The
# firmware runs at the sector rate, so every block is discretized with the Tustin rule, prewarped at its
# own key frequency (`M.controller`); the loop is evaluated at $z = e^{j\omega T_s}$.

# %%
fc, PM, fi = 1200.0, 40.0, 120.0
dA = M.design_lead(Kp, Td, fc, PM, Ts, fi=fi)
CA, _ = M.controller(Ts, dA['K'], dA['fz'], dA['fp'], fi)
print('lead + integral: phi_max = %.1f deg, zero %.0f Hz, pole %.0f Hz, K = %.4g A/track'
      % (dA['phi'], dA['fz'], dA['fp'], dA['K']))


def report(tag, ff, L):
    m = M.frf_margins(ff, L)
    gm = m['GM_up']
    print('%-28s crossover %s Hz, PM %s deg | GM %.1f dB at %.0f Hz | Ms = %.2f (%.1f dB) at %.0f Hz'
          % (tag, np.round(m['fc']).astype(int), np.round(m['PM'], 1), 20 * np.log10(gm[1]), gm[0],
             m['Ms'], 20 * np.log10(m['Ms']), m['f_Ms']))
    return m


LA_meas = M.frf(CA, f) * G
mA = report('no notch, MEASURED', f, LA_meas)
zA = M.cl_poles(P_fit, CA)
print('closed loop on the fitted model: max |z| = %.4f -> %s; unstable pole near %.0f Hz'
      % (abs(zA).max(), 'stable' if abs(zA).max() < 1 else 'UNSTABLE',
         abs(np.angle(zA[np.argmax(abs(zA))])) / (2 * np.pi * Ts)))

# %% [markdown]
# On the design model this loop has PM = 40 deg, but the measured loop gain rises above 0 dB again at the
# butterfly mode, with two more crossovers near 6 kHz. There the delay alone costs about 55 deg and the
# controller can not keep the phase of the mode under control. The Nyquist plot shows what happens: the
# resonance draws a loop large enough to encircle $-1$. The open loop has no poles outside the unit circle
# (the integrators sit on it and are indented around), so any encirclement means an unstable closed loop.
#
# ## A notch at the butterfly mode (L21)
# $N(s) = (s^2 + 2\zeta_z\omega_0 s + \omega_0^2)/(s^2 + 2\zeta_p\omega_0 s + \omega_0^2)$ with $\zeta_z/\zeta_p$ setting the
# depth: 20 dB at $f_0$ = the identified butterfly frequency, width $\zeta_p$ = 0.2. In its Tustin form the notch
# attenuates by at least 10 dB from 5.1 % below to 5.2 % above $f_0$, room for part of the spread of this mode
# between drives and with temperature. The notch lowers the loop gain at the mode, so the loop is
# gain-stabilized there; it also costs a few degrees at the crossover, so the lead is designed again with the
# notch included.

# %%
notch = (round(f_bfly), 20.0, 0.2)
d1 = M.design_lead(Kp, Td, fc, PM, Ts, fi=fi, notch=notch)
C1, _ = M.controller(Ts, d1['K'], d1['fz'], d1['fp'], fi, notch)
print('lead + integral + notch: notch %.0f Hz, 20 dB, zeta_p 0.2; phi_max = %.1f deg, zero %.0f Hz, pole %.0f Hz, '
      'K = %.4g A/track' % (notch[0], d1['phi'], d1['fz'], d1['fp'], d1['K']))
phase_notch = np.degrees(np.angle(M.frf(ct.c2d(M.notch_ct(*notch), Ts, 'tustin', prewarp_frequency=2 * np.pi * notch[0]), fc)[0]))
print('the notch costs %.1f deg at fc; the lead high-frequency gain is fp/fz = %.1f (%.0f dB)'
      % (-phase_notch, d1['fp'] / d1['fz'], 20 * np.log10(d1['fp'] / d1['fz'])))
L1_meas = M.frf(C1, f) * G
m1 = report('with notch, MEASURED', f, L1_meas)
m1f = report('with notch, fitted model', fd, M.frf(C1, fd) * Gd)
z1 = M.cl_poles(P_fit, C1)
print('closed loop on the fitted model: max |z| = %.4f -> %s' % (abs(z1).max(), 'stable' if abs(z1).max() < 1 else 'UNSTABLE'))
low = [(x, g) for x, g in zip(m1['f180'], m1['GM']) if g < 1]
print('lower gain margin (integral action): %.1f dB at %.0f Hz; the loop stays stable if the gain drops by less'
      % (20 * np.log10(low[0][1]), low[0][0]))
print('delay margin PM/wc = %.1f us = %.1f samples (the loop already has Td = %.1f us)'
      % (1e6 * m1['DM'][0], m1['DM'][0] / Ts, 1e6 * Td))

# %% [markdown]
# ## An independent check of the margins
# `python-control` computes the margins of the same loop from frequency-response data
# (`ct.frd`, `ct.stability_margins`), with its own crossover search and interpolation.

# %%
for tag, ff, L in (('measured', f, L1_meas), ('fitted model', fd, M.frf(C1, fd) * Gd)):
    gm, pm, sm, wpc, wgc, wms = ct.stability_margins(ct.frd(L, 2 * np.pi * ff))
    print('python-control, %-12s GM %.2f dB at %.0f Hz, PM %.2f deg at %.0f Hz, Ms = 1/sm = %.3f'
          % (tag, 20 * np.log10(gm), wpc / 2 / np.pi, pm, wgc / 2 / np.pi, 1 / sm))
print('frf_margins,    measured     GM %.2f dB at %.0f Hz, PM %.2f deg at %.0f Hz, Ms = %.3f'
      % (20 * np.log10(m1['GM_up'][1]), m1['GM_up'][0], m1['PM'][0], m1['fc'][0], m1['Ms']))

# %% [markdown]
# ## Bode and Nyquist plots on the measured response (L24, L26)

# %%
fig, ax = plt.subplots(2, 1, figsize=(8.5, 6.8), sharex=True)
for L, col, lab in ((LA_meas, 'C3', 'lead + integral, no notch'), (L1_meas, 'C0', 'with notch at %.0f Hz' % notch[0])):
    ax[0].semilogx(f, 20 * np.log10(abs(L)), col, lw=1.2, label=lab)
    p = np.degrees(np.unwrap(np.angle(L)))
    p -= 360 * np.round((p[np.argmin(abs(f - fc))] + 140) / 360)
    ax[1].semilogx(f, p, col, lw=1.2, label=lab)
ax[0].axhline(0, color='k', lw=0.6)
ax[0].set(ylabel='|L| (dB)', ylim=(-40, 50), title='Loop gain L = C G on the measured response')
ax[1].axhline(-180, color='k', lw=0.6)
ax[1].set(xlabel='frequency (Hz)', ylabel='phase (deg)', ylim=(-540, -60), xlim=(100, 17000))
ax[1].set_yticks([-540, -450, -360, -270, -180, -90])
ax[1].annotate('PM = %.0f deg' % m1['PM'][0], xy=(m1['fc'][0], -180 + m1['PM'][0]), xytext=(250, -100),
               arrowprops=dict(arrowstyle='->'), fontsize=8)
g = m1['GM_up']
ax[0].annotate('GM = %.1f dB' % (20 * np.log10(g[1])), xy=(g[0], -20 * np.log10(g[1])), xytext=(1800, -32),
               arrowprops=dict(arrowstyle='->'), fontsize=8)
ax[0].annotate('butterfly mode:\n|L| > 1 again', xy=(6200, 10), xytext=(1700, 25), color='C3',
               arrowprops=dict(arrowstyle='->', color='C3'), fontsize=8)
for a in ax:
    a.grid(True, which='both', alpha=0.4)
    a.legend(fontsize=8, loc='lower left')
fig.tight_layout()
plt.show()

fig, ax = plt.subplots(1, 2, figsize=(11, 5))
th = np.linspace(0, 2 * np.pi, 200)
for a, C_, L, col, ttl in ((ax[0], CA, LA_meas, 'C3', 'no notch: the resonance loop encircles -1'),
                           (ax[1], C1, L1_meas, 'C0', 'with notch: PM %.0f deg, GM %.1f dB, Ms %.2f'
                            % (m1['PM'][0], 20 * np.log10(m1['GM_up'][1]), m1['Ms']))):
    Lm = M.frf(C_, fd) * Gd
    a.plot(Lm.real, Lm.imag, col, lw=1.2, label='fitted model, f > 0')
    a.plot(Lm.real, -Lm.imag, col, ls='--', lw=0.6, alpha=0.6, label='fitted model, f < 0')
    a.plot(L.real, L.imag, 'k.', ms=3, label='measured points')
    a.plot(-1, 0, 'r+', ms=14, mew=2)
    a.plot(-1 + np.cos(th) / m1['Ms'], np.sin(th) / m1['Ms'], 'k:', lw=0.8, label='circle 1/Ms of the notched loop')
    a.set(xlim=(-3, 2), ylim=(-2.5, 2.5), xlabel='Re L', ylabel='Im L', title=ttl)
    a.set_aspect('equal')
    a.grid(True, alpha=0.4)
    a.legend(fontsize=7.5, loc='lower right')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Sensitivity and complementary sensitivity
# $S = 1/(1 + L)$ is the transfer function from track runout to PES; $T = L/(1 + L)$ from PES noise to the head
# motion. $|S|$ is small below the crossover (runout is followed), peaks at $M_s$ above it and returns to 1.
# Demo 3 uses this to explain the PES.

# %%
S1, T1 = 1 / (1 + M.frf(C1, fd) * Gd), 1 - 1 / (1 + M.frf(C1, fd) * Gd)
fig, ax = plt.subplots(figsize=(8.5, 4.3))
ax.semilogx(fd, 20 * np.log10(abs(S1)), 'C0', lw=1.3, label='|S|, fitted model')
ax.semilogx(f, 20 * np.log10(abs(1 / (1 + L1_meas))), 'C0.', ms=2.5, label='|S|, measured')
ax.semilogx(fd, 20 * np.log10(abs(T1)), 'C1', lw=1, label='|T|, fitted model')
ax.plot(m1['f_Ms'], 20 * np.log10(m1['Ms']), 'rv', ms=7, label='Ms = %.2f (%.1f dB) at %.0f Hz'
        % (m1['Ms'], 20 * np.log10(m1['Ms']), m1['f_Ms']))
ax.axhline(20 * np.log10(2), color='k', ls=':', lw=0.8, label='Ms = 2 (6 dB)')
ax.axhline(0, color='k', lw=0.6)
ax.set(xlabel='frequency (Hz)', ylabel='magnitude (dB)', xlim=(20, 18000), ylim=(-60, 15),
       title='Sensitivity S and complementary sensitivity T, lead + integral + notch')
ax.grid(True, which='both', alpha=0.4)
ax.legend(fontsize=8, loc='lower right')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## How robust is the notch?
# The butterfly frequency varies between drives and with temperature. Moving the mode of the fitted model by
# $\pm$5 % and $\pm$10 % with the controller unchanged shows how much of that the 20 dB, $\zeta_p$ = 0.2 notch tolerates.

# %%
for shift in (-0.10, -0.05, 0.05, 0.10):
    modes = [dict(m) for m in idf['modes']]
    modes[1]['f'] *= 1 + shift
    Ps = M.plant_dt(M.fitted_drive(Kp, Td, modes, Ts))
    ms = M.frf_margins(fd, M.frf(C1, fd) * M.frf(Ps, fd))
    zs = M.cl_poles(Ps, C1)
    if abs(zs).max() < 1:
        print('butterfly at %+3.0f %% (%.0f Hz): stable, GM %.1f dB at %.0f Hz, Ms %.2f'
              % (100 * shift, modes[1]['f'], 20 * np.log10(ms['GM_up'][1]), ms['GM_up'][0], ms['Ms']))
    else:
        print('butterfly at %+3.0f %% (%.0f Hz): UNSTABLE (max |z| = %.4f)' % (100 * shift, modes[1]['f'], abs(zs).max()))

# %% [markdown]
# A downward shift is the dangerous direction. Below $f_0$ the notch attenuates less and adds phase lag, which
# turns the resonance loop of the Nyquist plot toward $-1$: 5 % lower leaves only a few dB of gain margin and
# 10 % lower is unstable, while an upward shift of 10 % is harmless. If the butterfly frequency spreads by
# $\pm$5 % over a production lot, the notch should sit a few percent below the nominal mode or be wider, at
# the price of more phase lag at the crossover.
#
# **Reading the result.** The design model promised PM = 40 deg for both controllers. Without the notch the
# measured loop has three crossovers and is unstable, which only the resonance part of the Bode plot, or the
# Nyquist plot, reveals. With the notch the measured margins are PM about 40 deg, GM about 12 dB (now set at
# 3.7 kHz, between the pivot-bearing and butterfly modes) and $M_s \approx 1.7$, and `python-control`
# reproduces them from the same data. The delay margin of about 95 us, a little over three samples, says how much extra
# computation time the firmware could afford. Demo 3 asks what this loop does to the PES.
