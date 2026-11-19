# Author: Xu Chen, University of Washington
"""Sampled-data loops with the non-idealities of real hardware: saturation, noise, quantization, delay."""
import numpy as np
import control as ct


class PID:
    """Discrete PID, u = Kp e + Ki int(e) + Kd d/dt(-y), derivative on the measurement with a
    first-order filter (time constant Tf), and clamping anti-windup when u_lim = (lo, hi)."""

    def __init__(self, Kp, Ki=0.0, Kd=0.0, Ts=0.01, Tf=None, u_lim=None, antiwindup=True, u0=0.0):
        self.Kp, self.Ki, self.Kd, self.Ts = Kp, Ki, Kd, Ts
        self.Tf = Tf if Tf is not None else 10 * Ts
        self.u_lim, self.antiwindup, self.u0 = u_lim, antiwindup, u0
        self.reset()

    def reset(self):
        self.I, self.D, self.y_prev = 0.0, 0.0, None

    def __call__(self, r, y):
        e = r - y
        if self.y_prev is None:
            self.y_prev = y
        a = self.Tf / (self.Tf + self.Ts)
        self.D = a * self.D - self.Kd * (1 - a) * (y - self.y_prev) / self.Ts
        self.y_prev = y
        u_unsat = self.u0 + self.Kp * e + self.I + self.D
        u = u_unsat if self.u_lim is None else float(np.clip(u_unsat, *self.u_lim))
        # integrate unless saturated and the error pushes further into saturation
        if not (self.antiwindup and u != u_unsat and np.sign(e) == np.sign(u_unsat - u)):
            self.I += self.Ki * self.Ts * e
        return u


class LTIController:
    """Any continuous-time controller C(s) (e.g. a lead), discretized with Tustin; acts on e = r - y."""

    def __init__(self, C, Ts, u_lim=None, u0=0.0):
        self.sys = ct.ss(ct.c2d(ct.tf(C), Ts, 'tustin'))
        self.u_lim, self.u0 = u_lim, u0
        self.reset()

    def reset(self):
        self.x = np.zeros(self.sys.nstates)

    def __call__(self, r, y):
        e = r - y
        A, B, Cm, D = (np.asarray(m) for m in (self.sys.A, self.sys.B, self.sys.C, self.sys.D))
        u = float((Cm @ self.x).ravel()[0] + D.ravel()[0] * e) + self.u0
        self.x = A @ self.x + B.ravel() * e
        return u if self.u_lim is None else float(np.clip(u, *self.u_lim))


def simulate_loop(plant, controller, Ts, T, r, d=None, u_lim=None, noise_std=0.0, quant=None,
                  delay=0.0, x0=None, y0=0.0, seed=0):
    """Simulate plant (continuous LTI, input u + d) under a discrete controller(r, y_meas) -> u.

    r, d: callables of t or constants. u_lim: actuator limits (lo, hi). noise_std, quant: sensor noise
    and resolution. Timing as on hardware: u computed from the sample at t_k is applied from t_{k+1}
    (one sample of computation delay); delay (seconds, rounded to samples) adds transport delay on top.
    The plant works in deviation variables; y0 is added to the
    plant output (an operating point, e.g. ambient temperature). Returns dict of arrays t, r, y, y_meas, u."""
    rng = np.random.default_rng(seed)
    P = ct.ss(ct.c2d(ct.ss(plant), Ts, 'zoh'))
    A, B, C, D = (np.asarray(m) for m in (P.A, P.B, P.C, P.D))
    n = int(round(T / Ts)) + 1
    t = np.arange(n) * Ts
    f = lambda v, tk: v(tk) if callable(v) else (0.0 if v is None else v)
    x = np.zeros(P.nstates) if x0 is None else np.asarray(x0, float)
    buf = [0.0] * (int(round(delay / Ts)) + 1)
    out = {k: np.zeros(n) for k in ('r', 'y', 'y_meas', 'u')}
    for k in range(n):
        y = float((C @ x).ravel()[0] + D.ravel()[0] * buf[0]) + y0
        ym = y + noise_std * rng.standard_normal()
        if quant:
            ym = quant * np.round(ym / quant)
        rk = f(r, t[k])
        u = controller(rk, ym)
        if u_lim is not None:
            u = float(np.clip(u, *u_lim))
        buf.append(u)
        ua = buf.pop(0)
        x = A @ x + B.ravel() * (ua + f(d, t[k]))
        out['r'][k], out['y'][k], out['y_meas'][k], out['u'][k] = rk, y, ym, u
    out['t'] = t
    return out


def open_loop(plant, Ts, u, noise_std=0.0, quant=None, delay=0.0, y0=0.0, seed=0):
    """Response to a sampled input sequence u (array), with sensor noise/quantization and input delay."""
    u = np.asarray(u, float)
    seq = iter(u)
    res = simulate_loop(plant, lambda r, y: next(seq), Ts, (len(u) - 1) * Ts, 0.0,
                        noise_std=noise_std, quant=quant, delay=delay, y0=y0, seed=seed)
    return res['t'], res['y_meas'], res['y']
