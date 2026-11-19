# Author: Xu Chen, University of Washington
# %% [markdown]
# # Balancing robot 3: Nyquist with an unstable open loop, two gain margins, and the delay limit
#
# ![scene](figures/scene.png)
#
# **Overview.** The open loop $L = C G$ has one right-half-plane pole (P = 1), so the Nyquist criterion
# (L23-L24) asks for ONE counter-clockwise encirclement of -1, not zero. That turns the gain margin into
# a range: too little gain loses the encirclement (a LOWER limit), too much gain with the loop delay adds
# a clockwise one (an UPPER limit). The phase margin converts into the largest loop delay the robot can
# tolerate (L26). We check that number against simulations and against the logged tip-over.

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

p, b0, meta = M.identified_A()
Kp, Kd = meta['test']['Kp_Nm_per_rad'], meta['test']['Kd_Nms_per_rad']   # firmware gains in the logs
Tf, Ts, nd = meta['test']['Tf_s'], M.TS, meta['test']['loop_delay_samples']
G = ct.tf([b0], [1, 0, -p ** 2])
Td = nd * Ts + Ts / 2          # one sample of computation + half a sample for the zero-order hold
print('G(s) = %.0f/(s^2 - %.1f); firmware PD Kp = %.3f, Kd = %.4f, Tf = %.3f s; loop delay Td = %.1f ms'
      % (b0, p ** 2, Kp, Kd, Tf, 1e3 * Td))

# %% [markdown]
# ## The loop as the firmware runs it
# $L(j\omega) = C(e^{j\omega T_s})\,G(j\omega)\,e^{-j\omega T_d}$ with $T_d = n T_s + T_s/2$ (n samples of
# sensor + computation delay, half a sample for the zero-order hold). The firmware derivative is a
# filtered backward difference, $C(z) = K_p + K_d\frac{(1-a)(1 - z^{-1})}{T_s(1 - a z^{-1})}$,
# $a = T_f/(T_f+T_s)$; at 200 Hz it costs several degrees of phase at crossover compared with
# $K_p + K_d s/(T_f s + 1)$, so we evaluate it exactly (`M.pd_frf`).
#
# ## Nyquist plot, P = 1 (L23-L24)
# $Z = N + P$ with $N$ the clockwise encirclements of -1. For $Z = 0$ we need $N = -1$: one
# counter-clockwise encirclement. The plot starts at $L(0) = -K_p b_0/p^2$ (on the negative real axis,
# left of -1), while the D-term phase lead lifts the phase above $-180^\circ$ (so for $\omega > 0$ it runs through
# the lower half plane), and it spirals into the origin because of the delay.

# %%
w = np.logspace(-2, np.log10(np.pi / Ts), 40000)               # up to the Nyquist frequency
Lj = M.loop_frf(G, Kp, Kd, Tf, nd, w)
mg = M.loop_margins(Lj, w)
print('net counter-clockwise encirclements of -1: %.0f (P = 1 -> need 1)' % M.winding_ccw(Lj))


def nyq(ax, L, label, **kw):
    ax.plot(L.real, L.imag, label=label, **kw)
    ax.plot(L.real, -L.imag, ls='--', lw=0.8, color=kw.get('color'))


fig, ax = plt.subplots(1, 2, figsize=(11, 4.8))
for a, lim in ((ax[0], 1.15 * abs(Lj[0])), (ax[1], 2.0)):
    nyq(a, Lj, 'L(jw), firmware gains, Td = %.1f ms' % (1e3 * Td), color='C0')
    a.plot(-1, 0, 'r+', ms=14, mew=2)
    k = np.searchsorted(w, mg['wc'])
    a.annotate('', xy=(Lj[k + 40].real, Lj[k + 40].imag), xytext=(Lj[k].real, Lj[k].imag),
               arrowprops=dict(arrowstyle='->', color='C0'))
    th = np.linspace(0, 2 * np.pi, 200)
    a.plot(np.cos(th), np.sin(th), 'k:', lw=0.6)
    a.set(xlim=(-lim, 0.4 * lim), ylim=(-0.7 * lim, 0.7 * lim), xlabel='Re L', ylabel='Im L')
    a.set_aspect('equal')
    a.grid(True, alpha=0.4)
ax[0].set_title('Nyquist plot (solid w > 0, dashed w < 0)')
ax[1].set_title('Zoom near -1: unit circle dotted')
for wx, gm in mg['gm_list']:
    ax[0].plot(-1 / gm, 0, 'ks', ms=5)
ax[0].legend(fontsize=8, loc='lower left')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Gain margins (two of them) and phase margin
# The real-axis crossings of $L(j\omega)$ at $-1/\mathrm{GM}$ mark where a gain change would move a crossing
# across -1. At $\omega = 0$, $L(0) = -K_pb_0/p^2$: if the gain drops by more than the factor
# $p^2/(K_pb_0)$, the plot no longer encircles -1 (the Routh limit $K_p > p^2/b_0$ again). The second
# crossing, made by the delay, gives the upper limit.

# %%
lo = [g for wx, g in mg['gm_list'] if g < 1]
hi = [(wx, g) for wx, g in mg['gm_list'] if g > 1]
GM_lo = max(lo)
w180, GM_hi = min(hi, key=lambda x: x[1])
print('lower gain margin (w = 0): GM = %.3f = %.1f dB -> gain may drop to %.0f %% before it falls'
      % (GM_lo, 20 * np.log10(GM_lo), 100 * GM_lo))
print('upper gain margin (w = %.0f rad/s): GM = %.2f = %+.1f dB' % (w180, GM_hi, 20 * np.log10(GM_hi)))
print('phase margin PM = %.1f deg at crossover wc = %.1f rad/s (%.1f Hz)'
      % (mg['PM'], mg['wc'], mg['wc'] / 2 / np.pi))
# check the gain range with the exact sampled-data closed loop (stable iff all |z| < 1)
for f in (0.95 * GM_lo, 1.05 * GM_lo, 0.95 * GM_hi, 1.05 * GM_hi):
    zmax = abs(M.discrete_cl_poles(G, Kp, Kd, Tf, nd, gain=f)).max()
    print('  gain x %.3f: max |closed-loop pole z| = %.3f -> %s' % (f, zmax, 'stable' if zmax < 1 else 'UNSTABLE'))

# %% [markdown]
# ## Maximum loop delay (L26)
# A pure delay $T$ keeps $|L|$ and rotates the phase by $-\omega T$. At crossover the loop can absorb
# $\omega_c T_{extra} = \mathrm{PM}$ (in radians) before the plot touches -1:
#
# $$ T_{extra} = \frac{\mathrm{PM}\,[\mathrm{rad}]}{\omega_c}, \qquad T_{max} = T_d + T_{extra} $$

# %%
T_extra = mg['Tdm']
T_max = Td + T_extra
n_max = int(np.floor((T_max - Ts / 2) / Ts))
print('delay margin T_extra = %.1f ms -> total loop delay T_max = %.1f ms' % (1e3 * T_extra, 1e3 * T_max))
print('in whole samples of the 200 Hz loop: stable up to %d samples (%.0f ms), unstable at %d samples'
      % (n_max, 1e3 * n_max * Ts, n_max + 1))
for n in (n_max, n_max + 1):
    print('  exact sampled loop with %d samples: max |z| = %.4f' % (n, abs(M.discrete_cl_poles(G, Kp, Kd, Tf, n)).max()))
print('at the limit the loop oscillates at about wc = %.1f rad/s, period %.0f ms'
      % (mg['wc'], 1e3 * 2 * np.pi / mg['wc']))

# %% [markdown]
# ## Check by simulation, and the logged tip-over
# `simulate_loop` on the linear model (its controller output is applied one sample later, so
# `delay = (n - 1) Ts` gives n samples), then the nonlinear robot (sin theta, torque limit, IMU) from
# the same 2 deg tilt. Just past the limit the torque limit turns the growing oscillation into a large
# limit cycle; at 6 samples the robot falls. Finally the firmware-v2 log, where blocking telemetry made
# the delay 6 samples.

# %%
Gtrue = M.plant_tf(M.ROBOT_A)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for n, col in ((nd, 'C0'), (n_max, 'C2'), (n_max + 1, 'C3')):
    pid = sim.PID(Kp, Kd=Kd, Ts=Ts, Tf=Tf)
    res = sim.simulate_loop(Gtrue, pid, Ts, 1.5, 0.0, d=lambda t: -0.05 if t < 0.05 else 0.0,
                            delay=(n - 1) * Ts, noise_std=np.radians(0.02))
    ax[0].plot(res['t'], np.degrees(res['y']), col, label='%d samples (%.0f ms)' % (n, n * Ts * 1e3))
    o = M.simulate(M.ROBOT_A, sim.PID(Kp, Kd=Kd, Ts=Ts, Tf=Tf), 3.0, theta0_deg=2.0, n_delay=n, seed=7)
    ax[1].plot(o['t'], np.degrees(o['theta']), col, label='%d samples%s' % (n, ', FELL' if
                                                                         abs(o['theta'][-1]) > 0.7 else ''))
# well beyond the limit (the firmware-v2 delay): the oscillation grows past what the torque can catch
o = M.simulate(M.ROBOT_A, sim.PID(Kp, Kd=Kd, Ts=Ts, Tf=Tf), 3.0, theta0_deg=2.0, n_delay=6, seed=7)
ax[1].plot(o['t'], np.degrees(o['theta']), 'k', label='6 samples%s' % (', FELL' if abs(o['theta'][-1]) > 0.7 else ''))
ax[0].set(xlabel='time (s)', ylabel='theta (deg)', title='Linear model, torque tap at t = 0', ylim=(-8, 8))
ax[1].set(xlabel='time (s)', ylabel='theta (deg)', title='Nonlinear robot from 2 deg tilt', ylim=(-50, 50))
for a in ax:
    a.grid(True, alpha=0.4)
    a.legend(fontsize=8, title='loop delay')
fig.tight_layout()
plt.show()

tip, tm = data.load('05_balance_tipover')
tt, thl = tip['time_s'].values, tip['theta_deg'].values
# oscillation frequency of the growing oscillation (zero crossings of the logged pitch)
big = np.abs(thl) > 0.3                                          # above the noise
zc = tt[np.where(np.diff(np.sign(thl)) != 0)[0]]
zc = zc[zc > tt[np.argmax(big)]]
w_log = np.pi / np.median(np.diff(zc))
nbad = tm['test']['loop_delay_samples']
zc_ = M.discrete_cl_poles(G, Kp, Kd, Tf, nbad)
cl_bad = np.log(zc_[np.argmax(abs(zc_))]) / Ts                   # z = e^{sTs}
Lbad = M.loop_frf(G, Kp, Kd, Tf, nbad, w)
mbad = M.loop_margins(Lbad, w)
print('tip-over log: delay %d samples = %.0f ms > T_max = %.1f ms; fell at t = %.2f s'
      % (tm['test']['loop_delay_samples'], tm['test']['loop_delay_ms'], 1e3 * T_max, tt[-1]))
print('  model with this delay: PM = %.0f deg, CCW encirclements = %.0f (need 1) -> unstable pole '
      's = %.1f +/- j%.1f rad/s' % (mbad['PM'], M.winding_ccw(Lbad), cl_bad.real, abs(cl_bad.imag)))
print('  growing oscillation in the log: about %.0f rad/s (zero crossings, before saturation dominates)'
      % w_log)
fig, ax = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
ax[0].plot(tt, thl, lw=0.9, label='theta, firmware v2 (6-sample delay)')
ax[0].axhline(45, color='k', ls=':', lw=0.8, label='tilt cutoff')
ax[0].set(ylabel='theta (deg)', title='Logged tip-over (simulated log): oscillation grows, robot falls')
ax[1].plot(tt, tip['tau_cmd_Nm'], 'C1', lw=0.9, label='torque command')
ax[1].axhline(M.ROBOT_A['tau_max'], color='k', ls=':', lw=0.8, label='torque limit')
ax[1].axhline(-M.ROBOT_A['tau_max'], color='k', ls=':', lw=0.8)
ax[1].set(xlabel='time (s)', ylabel='tau (N m)')
for a in ax:
    a.grid(True, alpha=0.4)
    a.legend(fontsize=8, loc='upper left')
fig.tight_layout()
plt.show()
