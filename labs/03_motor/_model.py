# Author: Xu Chen, University of Washington
"""DC gearmotor with encoder: physics model, reduced model, log simulator and small helpers for Lab 03.

Angles and speeds are on the OUTPUT shaft (after the gearbox) unless a name says _m (motor side)."""
import numpy as np
import control as ct

# Motor A (demos): representative values for a 12 V brushed gearmotor class (37 mm can, 30:1 spur gearbox)
MOTOR_A = dict(R=2.4, L=2.0e-3, kt=0.019, Jm=3.3e-6, JL=2.0e-3, N=30, b=2.0e-6, Tc=2.0e-3,
               cpr_m=64, V_max=12.0, Ts=1e-3)
# Motor B (homework): a different unit, 50:1 gearbox, heavier load, coarser encoder. Datasheet-class values
# only: its Coulomb friction is not listed (estimate it from the log).
MOTOR_B = dict(R=1.9, L=1.1e-3, kt=0.0235, Jm=4.5e-6, JL=6.0e-3, N=50, b=3.0e-6,
               cpr_m=48, V_max=12.0, Ts=1e-3)
# R [ohm], L [H], kt = ke [N m/A = V s/rad], Jm rotor inertia and JL load inertia at the output [kg m^2],
# N gear ratio, b viscous friction at the motor [N m s/rad], Tc Coulomb friction at the motor [N m],
# cpr_m encoder counts per motor revolution (quadrature), V_max PWM voltage limit [V], Ts sample time [s]


def J_total(p):
    """Inertia seen by the motor shaft: rotor + load reflected through the gearbox (JL / N^2)."""
    return p['Jm'] + p['JL'] / p['N'] ** 2


def counts_per_rev(p):
    """Encoder counts per OUTPUT revolution."""
    return p['cpr_m'] * p['N']


def full_tf(p):
    """theta_out(s)/V(s) = (kt/N) / ( s [ (L s + R)(J s + b) + kt ke ] ), third order."""
    J = J_total(p)
    den = np.polymul([p['L'], p['R']], [J, p['b']])
    den = np.polyadd(den, [p['kt'] ** 2])
    return ct.tf([p['kt'] / p['N']], np.polymul(den, [1, 0]))


def reduced_params(p):
    """K/(s(tau s + 1)) by neglecting L: K [rad/s per V at output], tau [s]."""
    J, den = J_total(p), p['R'] * p['b'] + p['kt'] ** 2
    return p['kt'] / (p['N'] * den), p['R'] * J / den


def reduced_tf(K, tau):
    return ct.tf([K], [tau, 1, 0])


def simulate_log(p, V_of_t, T, seed, count0=0, sub=10):
    """Nonlinear motor (L, Coulomb friction, PWM limit, one-sample PWM update delay) sampled at Ts.

    Returns t, V_cmd (what the controller wrote), encoder counts (integer), and the true output angle."""
    rng = np.random.default_rng(seed)
    Ts, J = p['Ts'], J_total(p)
    n = int(round(T / Ts)) + 1
    t = np.arange(n) * Ts
    V_cmd = np.clip(np.array([V_of_t(tk) for tk in t]), -p['V_max'], p['V_max'])
    # small PWM amplitude ripple of the H-bridge supply (0.5 %), constant within a sample
    V_app = np.concatenate([[0.0], V_cmd[:-1]]) * (1 + 0.005 * rng.standard_normal(n))
    x = np.zeros(3)                                           # i [A], omega_m [rad/s], theta_m [rad]

    def f(x, V):
        i, w, _ = x
        # Coulomb friction, smoothed over 0.5 rad/s of motor speed so the integrator stays well posed
        dw = (p['kt'] * i - p['b'] * w - p['Tc'] * np.tanh(w / 0.5)) / J
        return np.array([(V - p['R'] * i - p['kt'] * w) / p['L'], dw, w])

    theta = np.zeros(n)
    h = Ts / sub
    for k in range(n):
        theta[k] = x[2] / p['N']
        for _ in range(sub):                                  # RK4 sub-steps (electrical pole is fast)
            k1 = f(x, V_app[k]); k2 = f(x + h / 2 * k1, V_app[k])
            k3 = f(x + h / 2 * k2, V_app[k]); k4 = f(x + h * k3, V_app[k])
            x = x + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    counts = count0 + np.floor(theta * counts_per_rev(p) / (2 * np.pi)).astype(int)
    return t, V_cmd, counts, theta


def speed_from_counts(t, counts, cpr, window=10, causal=False):
    """Output speed [rad/s] from encoder counts by a difference over `window` samples.

    causal=False: centered difference (offline, no lag). causal=True: backward difference (what a
    controller can compute in real time; lag of window/2 samples)."""
    th = 2 * np.pi * np.asarray(counts, float) / cpr
    Ts = np.median(np.diff(t))
    w = np.full_like(th, np.nan)
    if causal:
        w[window:] = (th[window:] - th[:-window]) / (window * Ts)
    else:
        h = window // 2
        w[h:-h] = (th[2 * h:] - th[:-2 * h]) / (2 * h * Ts)
    return w


def pd_by_angle(K, tau, sd):
    """PD C(s) = Kd (s + z) placing a closed-loop pole of K/(s(tau s+1)) at sd (angle + magnitude)."""
    G = reduced_tf(K, tau)
    phiG = np.angle(G(sd))                                    # angle of G at sd
    phi_z = np.pi - phiG                                      # zero must add this angle (mod 2 pi)
    phi_z = np.mod(phi_z + np.pi, 2 * np.pi) - np.pi
    z = -sd.real + sd.imag / np.tan(phi_z)
    Kd = 1 / abs((sd + z) * G(sd))
    return Kd, z


def lead_by_angle(K, tau, sd, z):
    """Lead C(s) = Kc (s + z)/(s + p): for a chosen zero z, the pole p that meets the angle condition."""
    G = reduced_tf(K, tau)
    phi_p = np.angle(sd + z) + np.angle(G(sd)) + np.pi        # angle the pole must subtract
    phi_p = np.mod(phi_p, 2 * np.pi)
    p = -sd.real + sd.imag / np.tan(phi_p)
    Kc = abs((sd + p) / ((sd + z) * G(sd)))
    return Kc, p


def step_specs(t, y, r_final, band=0.02):
    """Overshoot [%], 2 % settling time [s] and peak from a sampled step response to r_final."""
    y = np.asarray(y, float)
    os_ = 100 * max(np.max(y) / r_final - 1, 0.0)
    out = np.abs(y - r_final) > band * abs(r_final)
    ts = t[np.nonzero(out)[0][-1] + 1] if out.any() and np.nonzero(out)[0][-1] + 1 < len(t) else (
        t[0] if not out.any() else np.nan)
    return dict(OS=os_, Ts=ts, peak=float(np.max(y)))


def identify_from_log(name, window=10):
    """K [rad/s/V] and tau [s] of K/(s(tau s+1)) from a voltage-step log data/<name> (fit_fopdt on speed)."""
    from ctrlprac import data, ident
    df, meta = data.load(name)
    t, V, c = df['t (s)'].values, df['V_cmd (V)'].values, df['encoder (counts)'].values
    k_step = int(np.argmax(np.abs(np.diff(V)) > 0)) + 1
    w = speed_from_counts(t, c, meta['counts_per_rev'], window)
    m = ~np.isnan(w)
    fit = ident.fit_fopdt(t[m], w[m], V[k_step] - V[k_step - 1], t_step=t[k_step])
    return fit, (t, V, c, w, k_step, meta)
