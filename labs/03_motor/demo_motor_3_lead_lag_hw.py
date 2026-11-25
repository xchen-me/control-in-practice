# Author: Xu Chen, University of Washington
# %% [markdown]
# # Lab 03, demo 3: lead as the realizable PD, lag/PI for steady-state error, hardware test
#
# ![scene](figures/scene.png)
#
# **Overview.** A PD needs the derivative of the angle, but the angle comes from an encoder in steps of
# one count. Differentiating a staircase at 1 kHz turns every count into a voltage spike. A lead
# $C(s)=K_c(s+z)/(s+p)$ is the realizable PD: it adds phase where needed and caps the high-frequency
# gain at $K_c$ (L18). Then a lag (or PI) raises the low-frequency gain for a ramp reference and a
# constant load torque (L12, L18). Every design is tested in `simulate_loop` on the third-order plant
# with PWM saturation (12 V), a one-sample PWM delay and encoder quantization.

# %%
import os
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.modules.pop('_model', None)          # run_all.py runs several labs in one process
import numpy as np
import matplotlib.pyplot as plt
import control as ct
from ctrlprac.sim import PID, LTIController, simulate_loop
from _model import MOTOR_A, identify_from_log, reduced_tf, full_tf, counts_per_rev, pd_by_angle, \
    lead_by_angle, step_specs

p = MOTOR_A
fit, _ = identify_from_log('03_motor_step')
K, tau = fit['K'], fit['tau']
G, G3 = reduced_tf(K, tau), full_tf(p)          # design model, "hardware" model
s = ct.tf('s')
Ts, q, ulim = p['Ts'], 2 * np.pi / counts_per_rev(p), (-p['V_max'], p['V_max'])
r0, zeta, sigma = 0.5, 0.85, 28.0     # margin over demo 2 (zeta 0.80, sigma 27): see the text below
sd = complex(-sigma, sigma * np.tan(np.arccos(zeta)))
print('encoder resolution q = %.2e rad; a one-count step differentiated over Ts = %.2f rad/s' % (q, q / Ts))

# %% [markdown]
# ## Lead by the angle condition, zero on the plant pole
# Target with margin: on the hardware the 1 ms PWM delay, the electrical pole and the discretization add
# about 1 % of overshoot. For zeta = 0.80 the overshoot then exceeds the 2 % band, and the settling time
# jumps from the first entry into the band to the end of the overshoot (about 0.17 s, spec failed). So the
# lead targets zeta = 0.85 and sigma = 28 1/s.
#
# Put the lead zero on the mechanical pole, $z = 1/\tau$, then the angle condition at $s_d$ gives the pole
# $p$ and the magnitude condition gives $K_c$. With the cancellation the loop is
# $K_cK/(\tau s(s+p))$, so $p = 2\sigma$ exactly.

# %%
Kd, zpd = pd_by_angle(K, tau, sd)
z_l = 1 / tau
Kc, p_l = lead_by_angle(K, tau, sd, z_l)
C_lead = Kc * (s + z_l) / (s + p_l)
print('PD:   Kd = %.4f V s/rad, zero %.1f  (Kp = %.2f V/rad)' % (Kd, zpd, Kd * zpd))
print('lead: C(s) = %.2f (s + %.1f)/(s + %.1f);  HF gain %.1f V/rad, DC gain %.2f V/rad'
      % (Kc, z_l, p_l, Kc, Kc * z_l / p_l))

# %% Hardware-like test of PD (raw derivative), PD with a derivative filter, and the lead
T_end = 0.5
runs = {
    'PD, raw derivative of counts': PID(Kd * zpd, 0, Kd, Ts, Tf=Ts / 50),
    'PD, derivative filter Tf = 5 ms': PID(Kd * zpd, 0, Kd, Ts, Tf=5e-3),
    'lead (Tustin)': LTIController(C_lead, Ts),
}
res = {}
for name, ctrl in runs.items():
    # simulate_loop applies u(t_k) from t_(k+1): this is the 1 ms PWM update delay, so no extra delay
    out = simulate_loop(G3, ctrl, Ts, T_end, r0, u_lim=ulim, quant=q)
    sp = step_specs(out['t'], out['y'], r0)
    du = np.diff(out['u'][out['t'] > 0.02])           # skip the first 20 ms (reference step)
    res[name] = out
    print('%-34s OS = %4.1f %%, Ts = %.3f s, rms sample-to-sample du = %.2f V, max|u| = %.1f V'
          % (name, sp['OS'], sp['Ts'], np.sqrt(np.mean(du ** 2)), np.max(np.abs(out['u']))))

fig, ax = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
for name, out in res.items():
    ax[0].plot(out['t'], out['y'], label=name)
    ax[1].plot(out['t'], out['u'], lw=0.8, label=name)
ax[0].axhline(r0 * 1.05, color='k', ls='--', lw=0.8, label='5 % overshoot')
ax[0].set(ylabel='arm angle (rad)', title='0.5 rad step on the 3rd-order plant: 12 V limit, 1 ms delay, encoder')
ax[1].set(xlabel='time (s)', ylabel='PWM voltage (V)')
for a in ax:
    a.grid(True)
    a.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# The raw derivative turns each encoder count into a $K_d q/T_s$ jump of the voltage (about 1.8 V):
# audible PWM chatter and heating in a real motor. Filtering the derivative or using a lead removes it;
# a filtered PD **is** a lead. (The `PID` class differentiates the measurement, not the error, so its
# closed loop has no zero and overshoots less than the forward-path PD of demo 2.)
#
# ## Steady-state error: ramp reference and constant load torque (L12)
# The plant has an integrator (type 1): zero error to a step, but with a ramp $r=vt$ the error is
# $v/K_v$, $K_v = K\,C(0)$, and a constant load torque $T_L$ at the output enters like an input voltage
# $d = -T_L R/(N k_t)$, leaving $e_{ss} = |d|/C(0)$. A lag $(s+z_g)/(s+z_g/\beta)$ multiplies $C(0)$ by
# $\beta$; a PI factor $(s+z_g)/s$ makes $C(0)$ infinite and removes the load error entirely.

# %%
v, TL = 2.0, 0.3                           # ramp speed (rad/s), load torque at the output (N m)
d_load = -TL * p['R'] / (p['N'] * p['kt'])
e_ramp_spec = 0.01                         # rad
Kv = K * Kc * z_l / p_l
beta = np.ceil(v / (Kv * e_ramp_spec))
z_g = sigma / 10                           # lag zero a decade below the dominant poles
C_lag = C_lead * (s + z_g) / (s + z_g / beta)
C_pi = C_lead * (s + z_g) / s
print('lead only: Kv = %.1f 1/s, ramp error %.3f rad; load %.1f N m = %.2f V -> error %.3f rad'
      % (Kv, v / Kv, TL, d_load, abs(d_load) / (Kc * z_l / p_l)))
print('lag: beta = %.0f, zero %.1f, pole %.2f -> predicted ramp error %.4f rad, load error %.4f rad'
      % (beta, z_g, z_g / beta, v / (Kv * beta), abs(d_load) / (Kc * z_l / p_l * beta)))
print('closed-loop poles with lag:', np.round(ct.feedback(C_lag * G).poles(), 2))

r_ramp = lambda t: v * t
d_step = lambda t: d_load if t >= 0.5 else 0.0
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
for name, C in (('lead', C_lead), ('lead + lag', C_lag), ('lead + PI', C_pi)):
    o1 = simulate_loop(G3, LTIController(C, Ts), Ts, 3.0, r_ramp, u_lim=ulim, quant=q)
    o2 = simulate_loop(G3, LTIController(C, Ts), Ts, 3.0, 0.0, d=d_step, u_lim=ulim, quant=q)
    o3 = simulate_loop(G3, LTIController(C, Ts), Ts, 2.0, r0, u_lim=ulim, quant=q)
    e1, e2, sp = o1['r'] - o1['y'], o2['r'] - o2['y'], step_specs(o3['t'], o3['y'], r0)
    print('%-11s ramp error at 3 s = %.4f rad; load error at 3 s = %+.4f rad; 0.5 rad step: OS = %.1f %%, '
          'Ts = %.2f s' % (name, e1[-1], e2[-1], sp['OS'], sp['Ts']))
    ax[0].plot(o1['t'], e1, label=name)
    ax[1].plot(o2['t'], e2, label=name)
ax[0].axhline(e_ramp_spec, color='k', ls='--', lw=0.8, label='spec 0.01 rad')
ax[0].set(xlabel='time (s)', ylabel='tracking error (rad)', title='Ramp reference, %.0f rad/s' % v)
ax[1].set(xlabel='time (s)', ylabel='position error (rad)',
          title='Load torque %.1f N m from t = 0.5 s (hold at 0)' % TL)
for a in ax:
    a.grid(True)
    a.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# **Price of the lag/PI.** The step now overshoots about 13-14 % and settles in about 0.7 s: the lag adds a
# slow closed-loop pole near $-z_g$ that does not quite cancel the lag zero, and their small residue
# decays slowly (L9: a zero near a pole). A lower $z_g$ shrinks the overshoot but slows the recovery
# from the load (try $z_g = \sigma/40$). This is the classic lag trade-off.

# %% [markdown]
# ## Large move: the 12 V limit (saturation and integrator windup)
# For a 2 rad move the lead asks for $K_c r = 44$ V; the PWM gives 12 V. The loop is open while
# saturated. With an integrator (PI) the controller state keeps integrating the large error (windup)
# and the arm overshoots. The lag pole at $z_g/\beta \approx 0.25$ 1/s is almost an integrator on this
# time scale, so it winds up nearly as much.

# %%
fig, ax = plt.subplots(2, 1, figsize=(9, 5.5), sharex=True)
for name, C in (('lead', C_lead), ('lead + lag', C_lag), ('lead + PI', C_pi)):
    o = simulate_loop(G3, LTIController(C, Ts), Ts, 1.5, 2.0, u_lim=ulim, quant=q)
    sp = step_specs(o['t'], o['y'], 2.0)
    print('2 rad move, %-11s OS = %4.1f %%, Ts = %.3f s, time at the 12 V limit = %.0f ms'
          % (name, sp['OS'], sp['Ts'], 1e3 * Ts * np.sum(np.abs(o['u']) >= p['V_max'] - 1e-9)))
    ax[0].plot(o['t'], o['y'], label=name)
    ax[1].plot(o['t'], o['u'], label=name)
ax[0].set(ylabel='arm angle (rad)', title='2 rad move with the 12 V PWM limit')
ax[1].set(xlabel='time (s)', ylabel='PWM voltage (V)')
for a in ax:
    a.grid(True)
    a.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# **Takeaways.** The lag meets the ramp spec at the price of a slow tail; the PI drives the load-torque
# error to zero, to within about one encoder count (3.3e-3 rad), the quantization floor. The test with
# saturation, delay and quantization confirms the linear predictions for small moves; for large moves the
# 12 V limit dominates, and an integrator (or a slow lag) needs anti-windup, e.g. the clamping in
# `ctrlprac.sim.PID`.
