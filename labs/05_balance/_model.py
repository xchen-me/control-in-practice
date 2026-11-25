# Author: Xu Chen, University of Washington
"""Two-wheeled self-balancing robot, pitch loop only. Physics, nonlinear simulator with IMU and
complementary filter. Representative hobby-robot values (not a specific product).

Sign convention: theta > 0 when the body leans BACKWARD (toward -x). Then a positive wheel torque
tau (which drives the wheels toward +x) tips the body further back, and the linearized pitch model is

    theta(s) / tau(s) = G(s) = b0 / (s^2 - p^2),   p^2 = a b g / (a c - b^2),  b0 = (a + b/r) / (a c - b^2)

with a = m_w + I_w/r^2 + m_b,  b = m_b l,  c = I_b + m_b l^2 (wheels free to roll, no slip).
The usual negative feedback tau = -C(s) theta balances it."""
import numpy as np
import control as ct

G0 = 9.81
TS = 0.005                    # 200 Hz control loop

# demo robot "R1": 1.0 kg body, centre of mass 10 cm above the axle, 68 mm wheels. I_b stays fixed when
# m_b or l is changed (a change of m_b is mass at the centre of mass). n_delay: samples of sensor +
# computation delay; with the zero-order hold the total loop delay is Td = (n_delay + 1/2) TS.
ROBOT_A = dict(name='R1', m_b=1.0, l=0.10, I_b=1.0 * (0.22 ** 2 + 0.06 ** 2) / 12, m_w=0.10, r=0.034,
               tau_max=0.50, h_push=0.20, n_delay=1)

IMU = dict(gyro_noise_dps=0.10, gyro_bias_dps=0.4, acc_noise_deg=0.30, alpha=0.98)
LOG_RES_DEG = 0.01            # logged angle resolution
TAU_RES = 1.0 / 255           # 8-bit PWM torque command, as a fraction of tau_max


def coeffs(rb):
    Iw = 0.5 * rb['m_w'] * rb['r'] ** 2
    a = rb['m_w'] + Iw / rb['r'] ** 2 + rb['m_b']
    b = rb['m_b'] * rb['l']
    c = rb['I_b'] + rb['m_b'] * rb['l'] ** 2
    return a, b, c


def pitch_params(rb):
    """Return p (rad/s, unstable pole) and b0 (rad/s^2 per N m) of G(s) = b0/(s^2 - p^2)."""
    a, b, c = coeffs(rb)
    d = a * c - b * b
    return np.sqrt(a * b * G0 / d), (a + b / rb['r']) / d


def plant_tf(rb):
    p, b0 = pitch_params(rb)
    return ct.tf([b0], [1, 0, -p ** 2])


def pd_tf(Kp, Kd, Tf):
    """PD with first-order filtered derivative, Kp + Kd s/(Tf s + 1), as implemented by ctrlprac.sim.PID."""
    return ct.tf([Kd + Kp * Tf, Kp], [Tf, 1])


def _deriv(rb, s, tau, F):
    """Nonlinear wheel + body equations, s = [x, xdot, theta, thetadot]; F = horizontal push at h_push."""
    a, b, c = coeffs(rb)
    x, xd, th, thd = s
    # [a, -b cos; -b cos, c] [xdd; thdd] = [tau/r - b sin thd^2 + F; tau + b g sin - F h cos]
    M = np.array([[a, -b * np.cos(th)], [-b * np.cos(th), c]])
    f = np.array([tau / rb['r'] - b * np.sin(th) * thd ** 2 + F,
                  tau + b * G0 * np.sin(th) - F * rb['h_push'] * np.cos(th)])
    xdd, thdd = np.linalg.solve(M, f)
    return np.array([xd, xdd, thd, thdd]), xdd


def simulate(rb, controller, T, push=None, theta0_deg=0.0, n_delay=None, motors_on=True,
             stop_deg=None, seed=0, fall_deg=45.0):
    """Nonlinear robot under a discrete controller(r, theta_hat) -> tau, with IMU noise, complementary
    filter, torque limit and PWM quantization, n_delay samples of sensor+computation delay (with the
    zero-order hold, Td = (n_delay + 1/2) TS).

    push(t) -> horizontal force (N) at h_push. The run stops when |theta| > fall_deg (the firmware's
    tilt cutoff, a fall; the logs end there) or, if given, when |theta| > stop_deg (the stand of the
    drop test). Returns dict of arrays."""
    rng = np.random.default_rng(seed)
    n_delay = rb['n_delay'] if n_delay is None else n_delay
    n = int(round(T / TS)) + 1
    sub = 10
    h = TS / sub
    s = np.array([0.0, 0.0, np.radians(theta0_deg), 0.0])
    th_hat = s[2]
    bias = np.radians(IMU['gyro_bias_dps'])
    buf = [0.0] * n_delay
    out = {k: [] for k in ('t', 'theta', 'theta_hat', 'gyro', 'acc_tilt', 'tau_cmd', 'tau', 'x', 'push')}
    for k in range(n):
        t = k * TS
        F = push(t) if push else 0.0
        # IMU: gyro (rate + bias + noise) and accelerometer tilt (noise; the error from forward
        # acceleration is left out here, the complementary filter weights it by only 1 - alpha)
        gyro = s[3] + bias + np.radians(IMU['gyro_noise_dps']) * rng.standard_normal()
        acc = s[2] + np.radians(IMU['acc_noise_deg']) * rng.standard_normal()
        al = IMU['alpha']
        # complementary filter; trapezoidal gyro integration keeps it free of phase error
        th_hat = al * (th_hat + 0.5 * (gyro + g_prev) * TS) + (1 - al) * acc if k else acc
        g_prev = gyro
        u = controller(0.0, th_hat) if motors_on else 0.0
        u = float(np.clip(u, -rb['tau_max'], rb['tau_max']))
        q = TAU_RES * rb['tau_max']
        u = q * np.round(u / q)
        buf.append(u)
        ua = buf.pop(0)
        for key, v in zip(out, (t, s[2], th_hat, gyro, acc, u, ua, s[0], F)):
            out[key].append(v)
        if abs(s[2]) > np.radians(fall_deg) or (stop_deg and abs(s[2]) > np.radians(stop_deg)):
            break
        for _ in range(sub):            # RK4 with the torque held (zero-order hold)
            k1, _ = _deriv(rb, s, ua, F)
            k2, _ = _deriv(rb, s + h / 2 * k1, ua, F)
            k3, _ = _deriv(rb, s + h / 2 * k2, ua, F)
            k4, _ = _deriv(rb, s + h * k3, ua, F)
            s = s + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return {k: np.array(v) for k, v in out.items()}


def loop_margins(L_jw, w):
    """Margins from samples of L(jw) on a dense grid w (rad/s), delay included in L_jw.

    Returns dict(wc, PM, Tdm, gm_list) where gm_list holds (w, GM) at every -180 deg crossing
    (w = 0 included when L(0) is negative real): GM < 1 means the gain may only be REDUCED by that
    factor (a lower gain limit), GM > 1 an upper limit. Tdm = PM/wc is the extra delay the loop tolerates."""
    mag = np.abs(L_jw)
    ph = np.unwrap(np.angle(L_jw))
    if ph[0] > 0:
        ph -= 2 * np.pi
    i = np.where(np.diff(np.sign(mag - 1)) != 0)[0]
    wc = [np.interp(1, mag[[k + 1, k]], w[[k + 1, k]]) if mag[k] > mag[k + 1] else
          np.interp(1, mag[[k, k + 1]], w[[k, k + 1]]) for k in i]
    pm = [180 + np.degrees(np.interp(x, w, ph)) for x in wc]
    gms = []
    if abs(np.degrees(ph[0]) + 180) < 1:                       # L(0) on the negative real axis
        gms.append((0.0, 1 / mag[0]))
    j = np.where(np.diff(np.sign(np.degrees(ph) + 180)) != 0)[0]
    for k in j:
        if abs(np.degrees(ph[k]) + 180) < 1 and k < 3:
            continue
        wx = np.interp(-np.pi, ph[[k, k + 1]] if ph[k] < ph[k + 1] else ph[[k + 1, k]],
                       w[[k, k + 1]] if ph[k] < ph[k + 1] else w[[k + 1, k]])
        gms.append((wx, 1 / np.interp(wx, w, mag)))
    k = int(np.argmin(pm)) if pm else None
    return dict(wc=wc[k] if pm else np.nan, PM=pm[k] if pm else np.nan,
                Tdm=np.radians(pm[k]) / wc[k] if pm else np.nan, gm_list=gms)


def winding_ccw(L_jw):
    """Net counter-clockwise encirclements of -1 by L(jw), w from -inf to +inf (conjugate mirror).
    The grid ends at a finite w (pi/Ts in the demos), where |L| is small but not zero, so the curve is
    closed from L(j w_max) back to its mirror image before the angle is counted; the result is an integer."""
    z = np.concatenate([np.conj(L_jw[::-1]), L_jw, np.conj(L_jw[-1:])]) + 1
    ang = np.unwrap(np.angle(z))
    return float(np.round((ang[-1] - ang[0]) / (2 * np.pi))) + 0.0      # no "-0"


def fit_drop(t, th_deg, max_deg=10.0):
    """Drop test fit theta = A e^{pt} + B e^{-pt} + c on |theta| < max_deg: linear least squares in
    (A, B, c) for each trial p, keep the p with the smallest residual. Returns p, (A, B, c), mask."""
    m = np.abs(th_deg) < max_deg

    def fit(p):
        A = np.column_stack([np.exp(p * t[m]), np.exp(-p * t[m]), np.ones(m.sum())])
        x = np.linalg.lstsq(A, th_deg[m], rcond=None)[0]
        return np.sum((A @ x - th_deg[m]) ** 2), x

    ps = np.linspace(2, 40, 7601)
    p = ps[np.argmin([fit(q)[0] for q in ps])]
    return p, fit(p)[1], m


def fit_closed_loop(df, n_delay, windows):
    """Integrated equation-error fit of theta'' = p^2 theta + b0 tau(t - Td) on push-free windows:
    gyro(t) - gyro(t0) = p^2 int theta + b0 int tau + c (t - t0). Returns p, b0, regressor matrix."""
    t = df['time_s'].values
    th = np.radians(df['theta_deg'].values)
    gy = np.radians(df['gyro_dps'].values)
    tau = np.r_[np.zeros(n_delay), df['tau_cmd_Nm'].values[:-n_delay]]   # torque actually applied
    integ = lambda v: np.r_[0, np.cumsum(v[:-1])] * TS                   # integral up to sample k
    rows, rhs = [], []
    for t0, t1 in windows:
        k = (t >= t0) & (t <= t1)
        rows.append(np.column_stack([integ(th[k]), integ(tau[k]), t[k] - t0]))
        rhs.append(gy[k] - gy[k][0])
    A = np.vstack(rows)
    p2, b0, _ = np.linalg.lstsq(A, np.concatenate(rhs), rcond=None)[0]
    return np.sqrt(p2), b0, A


PUSH_FREE_A = ((0.0, 2.95), (3.2, 6.45), (6.7, 10.0))


def identified_A():
    """The model demo 1 identifies from the R1 push log (used by demos 2 and 3)."""
    from ctrlprac import data
    df, meta = data.load('05_balance_push')
    p, b0, _ = fit_closed_loop(df, meta['test']['loop_delay_samples'], PUSH_FREE_A)
    return p, b0, meta


def pd_frf(Kp, Kd, Tf, Ts, w):
    """Frequency response of the firmware PD exactly as ctrlprac.sim.PID computes it (backward-difference
    derivative with a first-order filter, a = Tf/(Tf + Ts)), acting on e = -theta:
    C(z) = Kp + Kd (1 - a)(1 - z^-1) / (Ts (1 - a z^-1)),  evaluated at z = e^{jwTs}."""
    a = Tf / (Tf + Ts)
    zi = np.exp(-1j * w * Ts)
    return Kp + Kd * (1 - a) * (1 - zi) / (Ts * (1 - a * zi))


def loop_frf(G, Kp, Kd, Tf, n_delay, w, Ts=TS):
    """L(jw) of the sampled loop: firmware PD, n_delay samples of delay, zero-order hold (Ts/2 delay),
    so the total loop delay is Td = (n_delay + 1/2) Ts."""
    return pd_frf(Kp, Kd, Tf, Ts, w) * G(1j * w) * np.exp(-1j * w * (n_delay + 0.5) * Ts)


def discrete_cl_poles(G, Kp, Kd, Tf, n_delay, Ts=TS, gain=1.0):
    """Exact closed-loop poles (z-plane) of the sampled loop: ZOH plant, PD, n_delay samples of delay."""
    a = Tf / (Tf + Ts)
    Cz = Kp + Kd * (1 - a) / Ts * ct.tf([1, -1], [1, -a], Ts)
    L = gain * Cz * ct.c2d(G, Ts, 'zoh') * ct.tf([1], [1] + [0] * n_delay, Ts)
    return ct.poles(ct.feedback(L, 1))


def design_pd(p, b0, zeta=0.7, wn=None):
    """Ideal-PD pole placement for b0/(s^2 - p^2): s^2 + b0 Kd s + (b0 Kp - p^2) = s^2 + 2 zeta wn s + wn^2."""
    wn = 2 * p if wn is None else wn
    return (wn ** 2 + p ** 2) / b0, 2 * zeta * wn / b0
