# Author: Xu Chen, University of Washington
# %% [markdown]
# # Balancing robot 1: an unstable plant, identified from logs; why P fails, Routh range for PD
#
# ![scene](figures/scene.png)
#
# **Overview.** A two-wheeled robot balances by driving its wheels under its centre of mass. Its pitch
# dynamics near upright are an inverted pendulum: one pole in the right half plane (L3, L10). We
# (1) write the linearized model, (2) identify the unstable pole p and the input gain b0 from two logs
# (a motors-off drop test and a closed-loop push test), (3) show with the root locus that P control can
# never stabilize it, and (4) use Routh-Hurwitz to find the stabilizing PD gains (L10-L11).

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import control as ct
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
from ctrlprac import data
import _model as M

# %% [markdown]
# ## Mechanism and model
# Wheel torque tau rolls the wheels and reacts on the body. Eliminating the wheel position from the
# two linearized equations of motion (wheels free to roll) leaves the pitch transfer function
#
# $$ G(s) = \frac{\Theta(s)}{T(s)} = \frac{b_0}{s^2 - p^2}, \qquad
#    p^2 = \frac{a b g}{a c - b^2},\quad b_0 = \frac{a + b/r}{a c - b^2} $$
#
# with $a = m_w + I_w/r^2 + m_b$, $b = m_b l$, $c = I_b + m_b l^2$. Poles at $\pm p$: one is unstable.

# %%
p_true, b0_true = M.pitch_params(M.ROBOT_A)
print('robot R1 (representative values): m_b = %.2f kg, l = %.3f m, I_b = %.4f kg m^2, r = %.3f m'
      % tuple(M.ROBOT_A[k] for k in ('m_b', 'l', 'I_b', 'r')))
print('model from parameters: p = %.2f rad/s, b0 = %.0f rad/s^2/(N m)' % (p_true, b0_true))

# %% [markdown]
# ## Identify p from the drop test
# With the motors off, $\ddot\theta = p^2\theta$, so $\theta(t) = A e^{pt} + B e^{-pt}$ (released from
# rest, $\theta = \theta_0\cosh pt$). For a trial $p$ this is LINEAR in $(A, B)$: solve by least squares,
# keep the $p$ with the smallest residual (a one-dimensional search). Use $|\theta| < 10$ deg, where
# $\sin\theta \approx \theta$. (A straight line through $\ln|\theta|$ would ignore the $e^{-pt}$ term and
# underestimate $p$.)

# %%
drop, meta = data.load('05_balance_drop')
t, th = drop['time_s'].values, drop['theta_deg'].values
p_drop, coef, m = M.fit_drop(t, th)
slope = np.polyfit(t[m & (np.abs(th) > 3)], np.log(np.abs(th[m & (np.abs(th) > 3)])), 1)[0]
print('drop test: p = %.2f rad/s (true %.2f, error %+.1f %%); straight line on ln|theta|: %.2f rad/s'
      % (p_drop, p_true, 100 * (p_drop / p_true - 1), slope))
print('time to double the tilt = ln2/p = %.0f ms' % (1e3 * np.log(2) / p_drop))

fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
ax[0].plot(t, th, '.', ms=4, label='log: theta (complementary filter)')
tt = np.linspace(0, t[m][-1], 200)
fit = coef[0] * np.exp(p_drop * tt) + coef[1] * np.exp(-p_drop * tt) + coef[2]
ax[0].plot(tt, fit, 'r', label='fit $Ae^{pt}+Be^{-pt}$, p = %.1f 1/s' % p_drop)
ax[0].set(xlabel='time (s)', ylabel='theta (deg)', title='Drop test, motors off')
ax[1].semilogy(t, np.abs(th), '.', ms=4, label='|theta| log')
ax[1].semilogy(tt, np.abs(fit), 'r', label='fit')
ax[1].axhline(10, color='k', ls=':', lw=0.8, label='fit limit 10 deg')
ax[1].set(xlabel='time (s)', ylabel='|theta| (deg)', title='Exponential growth (log scale)')
for a in ax:
    a.grid(True, which='both', alpha=0.4)
    a.legend(fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Identify b0 (and p again) from the closed-loop push log
# An unstable plant cannot be step-tested open loop; we use the balancing log. Integrating
# $\ddot\theta = p^2\theta + b_0\tau(t - T_d)$ over a window that excludes the pushes gives a linear
# least-squares problem in $(p^2, b_0)$ without differentiating noisy data:
#
# $$ \omega_{gyro}(t) - \omega_{gyro}(t_0) = p^2\!\int_{t_0}^{t}\!\theta\,dt + b_0\!\int_{t_0}^{t}\!\tau\,dt + c\,(t - t_0) $$
#
# ($c$ absorbs the small angle offset from gyro bias). The torque is the logged command shifted by the
# loop delay (meta: one sample). In closed loop $\tau \approx -K_p\theta - K_d\dot\theta$ is almost a copy of
# $\theta$, so the two regressors are nearly collinear; the firmware adds a small PRBS dither to the torque
# during this test to make $p^2$ and $b_0$ separable (standard closed-loop identification practice).

# %%
push, pm = data.load('05_balance_push')
t = push['time_s'].values
nd = pm['test']['loop_delay_samples']
p_id, b0_id, A = M.fit_closed_loop(push, nd, M.PUSH_FREE_A)        # see _model.py: 12 lines of numpy
print('push log: p = %.2f rad/s (true %.2f, %+.1f %%), b0 = %.0f (true %.0f, %+.1f %%)'
      % (p_id, p_true, 100 * (p_id / p_true - 1), b0_id, b0_true, 100 * (b0_id / b0_true - 1)))
print('regressor correlation theta-integral vs tau-integral: %.2f' % np.corrcoef(A[:, 0], A[:, 1])[0, 1])
G = ct.tf([b0_id], [1, 0, -p_id ** 2])
print('identified G(s) = %.0f / (s^2 - %.1f)' % (b0_id, p_id ** 2))

fig, ax = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
ax[0].plot(t, push['theta_deg'], lw=0.8, label='theta (deg)')
ax[0].axvspan(3.0, 3.1, color='C3', alpha=0.3, label='push forward')
ax[0].axvspan(6.5, 6.6, color='C2', alpha=0.3, label='push backward')
ax[0].set_ylabel('theta (deg)')
ax[1].plot(t, push['tau_cmd_Nm'], lw=0.6, color='C1', label='torque command (PD + dither)')
ax[1].axhline(M.ROBOT_A['tau_max'], color='k', ls=':', lw=0.8, label='torque limit')
ax[1].axhline(-M.ROBOT_A['tau_max'], color='k', ls=':', lw=0.8)
ax[1].set(xlabel='time (s)', ylabel='tau (N m)')
ax[0].set_title('Logged balancing run with two pushes (simulated log, 200 Hz)')
for a in ax:
    a.grid(True, alpha=0.4)
    a.legend(loc='upper right', fontsize=8)
fig.tight_layout()
plt.show()

# %% [markdown]
# ## Why P alone fails (root locus, L13)
# With $C = K$ the characteristic equation is $s^2 + (K b_0 - p^2) = 0$: the two poles at $\pm p$ meet at
# the origin when $K = p^2/b_0$ and then move up and down the imaginary axis. No gain puts both in the
# open left half plane: at best an undamped oscillation, and any delay tips it over.

# %%
fig = plt.figure(figsize=(5.5, 4.5))
Ks = np.linspace(0, 4 * p_id ** 2 / b0_id, 400)
for K in Ks[::8]:
    r = np.roots([1, 0, K * b0_id - p_id ** 2])
    plt.plot(r.real, r.imag, 'b.', ms=3)
plt.plot([p_id, -p_id], [0, 0], 'kx', ms=10, label='open-loop poles $\\pm p$')
plt.axvline(0, color='k', lw=0.6)
plt.xlabel('Re(s) (rad/s)')
plt.ylabel('Im(s) (rad/s)')
plt.title('Root locus with P control: never in the open LHP')
plt.grid(True, alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()
print('P control: K = p^2/b0 = %.3f N m/rad puts both poles at s = 0; larger K gives s = +/- j sqrt(K b0 - p^2)'
      % (p_id ** 2 / b0_id))

# %% [markdown]
# ## Routh range for PD (L10, L11)
# The firmware PD is $C(s) = K_p + K_d s/(T_f s + 1)$. Closed-loop characteristic polynomial:
#
# $$ T_f s^3 + s^2 + \big(b_0(K_d + K_p T_f) - p^2 T_f\big)s + (b_0 K_p - p^2) = 0 $$
#
# Routh rows: $s^3: T_f,\ b_0(K_d+K_pT_f)-p^2T_f$; $s^2: 1,\ b_0K_p - p^2$;
# $s^1: b_0(K_d+K_pT_f)-p^2T_f - T_f(b_0K_p-p^2) = b_0 K_d$; $s^0: b_0K_p - p^2$.
# Stable iff $K_p > p^2/b_0$ and $K_d > 0$: the D term is what makes balancing possible. There is no
# upper limit here because delay is not yet in the model; demo 3 adds it.

# %%
Kp_min = p_id ** 2 / b0_id
print('Routh: stable iff Kp > p^2/b0 = %.3f N m/rad and Kd > 0' % Kp_min)
Tf = pm['test']['Tf_s']
for Kp, Kd in ((0.8 * Kp_min, 0.05), (1.38, 0.0), (1.38, 0.0585)):
    poles = np.roots([Tf, 1, b0_id * (Kd + Kp * Tf) - p_id ** 2 * Tf, b0_id * Kp - p_id ** 2])
    print('  Kp = %.3f, Kd = %.4f: max Re(pole) = %+.2f rad/s' % (Kp, Kd, poles.real.max()))
