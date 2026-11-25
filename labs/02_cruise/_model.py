# Author: Xu Chen, University of Washington
"""Lab 02 helpers: vehicle parameters, the nonlinear longitudinal model that generates the logs,
its linearization, and root-locus / margin utilities."""
import numpy as np
import control as ct

G0, RHO = 9.81, 1.2          # m/s^2, kg/m^3
TS = 0.1                     # s, cruise-control loop rate
U_LIM = (0.0, 100.0)         # % throttle (engine braking and brakes are not used by this cruise controller)
WHEEL_Q, WHEEL_N = 0.1 / 3.6, 0.02     # wheel-speed resolution (0.1 km/h) and noise std (m/s)
GPS_N = 0.12                            # GPS Doppler speed noise std (m/s)

# Representative values (not a specific vehicle's datasheet). Fmax is the tractive force at full throttle in
# top gear near highway speed (no downshift): engine torque x overall ratio x efficiency / wheel radius.
CAR_A = dict(name='mid-size sedan (demo)', m=1550.0, m_eff=1610.0, Cd=0.29, A=2.30, Crr=0.010,
             Fmax=1600.0, tau_e=0.6, delay=0.1, v0=27.0)


def road_load(c, v, grade=0.0):
    """Aerodynamic drag + rolling resistance + grade force (N) at speed v and grade (rise/run)."""
    th = np.arctan(grade)
    return 0.5 * RHO * c['Cd'] * c['A'] * v * abs(v) + c['Crr'] * c['m'] * G0 * np.cos(th) \
        + c['m'] * G0 * np.sin(th)


def trim_throttle(c, v):
    """Throttle (%) that holds speed v on a flat road."""
    return 100 * road_load(c, v) / c['Fmax']


def linearize(c, v0=None):
    """Linear model about v0 on a flat road. Returns dict(b, K, tau, G) with
    G(s) = (Fmax/100) / ((m_eff s + b)(tau_e s + 1)) in (m/s) per % throttle (delay not included)."""
    v0 = c['v0'] if v0 is None else v0
    b = RHO * c['Cd'] * c['A'] * v0                 # d(drag)/dv; rolling resistance is constant
    K, tau = c['Fmax'] / 100 / b, c['m_eff'] / b
    G = ct.tf([K], np.polymul([tau, 1], [c['tau_e'], 1]))
    return dict(b=b, K=K, tau=tau, G=G, u0=trim_throttle(c, v0))


def grade_to_pct(c, grade):
    """Grade (rise/run) expressed as the equivalent throttle (%) that it removes at the plant input."""
    return -100 * c['m'] * G0 * np.sin(np.arctan(grade)) / c['Fmax']


def drive(c, T, controller=None, u_seq=None, r=None, grade=0.0, v_init=None, seed=0, wind_std=15.0):
    """Nonlinear longitudinal simulation, integrated at 1 ms, logged at TS.

    Either u_seq (throttle array, open loop) or controller(r, v_meas) -> throttle (closed loop on the
    wheel-speed signal). grade: constant or callable of t. Wind gusts: low-pass noise force (N std).
    Timing as in sim.simulate_loop: u computed from the sample at t_k reaches the engine one sample later,
    after the delay. Returns dict t, r, u, v (true), v_wheel, v_gps, grade, all at the sample times t."""
    rng = np.random.default_rng(seed)
    n, sub = int(round(T / TS)) + 1, 100
    dt = TS / sub
    gr = grade if callable(grade) else (lambda t: grade)
    rr = r if callable(r) else (lambda t: r)
    v = c['v0'] if v_init is None else v_init
    F = road_load(c, v)                         # engine force state starts at trim
    nd = int(round(c['delay'] / TS))
    buf = [100 * F / c['Fmax']] * (nd + 1)
    w, a_w = 0.0, np.exp(-TS / 5.0)             # gust force, 5 s correlation time
    out = {k: np.zeros(n) for k in ('t', 'r', 'u', 'v', 'v_wheel', 'v_gps', 'grade')}
    for k in range(n):
        t, vk = k * TS, v
        vw = WHEEL_Q * np.round((v + WHEEL_N * rng.standard_normal()) / WHEEL_Q)
        vg = v + GPS_N * rng.standard_normal()
        rk = rr(t) if controller is not None else np.nan
        u = float(np.clip(controller(rk, vw) if controller is not None else u_seq[k], *U_LIM))
        buf.append(u)
        ua = buf.pop(0)
        w = a_w * w + wind_std * np.sqrt(1 - a_w ** 2) * rng.standard_normal()
        g = gr(t)
        for _ in range(sub):
            F += dt / c['tau_e'] * (c['Fmax'] * ua / 100 - F)
            v += dt / c['m_eff'] * (F - road_load(c, v, g) - w)
        for key, val in (('t', t), ('r', rk), ('u', u), ('v', vk), ('v_wheel', vw), ('v_gps', vg), ('grade', g)):
            out[key][k] = val
    return out


def v_steady(c, u):
    """Flat-road steady speed for throttle u (%)."""
    Fa = c['Fmax'] * u / 100 - c['Crr'] * c['m'] * G0
    return np.sqrt(Fa / (0.5 * RHO * c['Cd'] * c['A']))


def dominant(poles):
    poles = np.asarray(poles)
    return poles[np.argmax(poles.real)]


def gain_for_zeta(L, zeta, kmax=1e4):
    """Smallest root-locus gain k for which the slowest complex closed-loop pair of 1 + k L(s) has damping
    zeta. Returns (k, closed-loop poles)."""
    num, den = np.ravel(L.num[0][0]), np.ravel(L.den[0][0])

    def z_of(k):
        p = np.roots(np.polyadd(den, k * num))
        p = p[np.abs(p.imag) > 1e-9]
        return 1.0 if not len(p) else -dominant(p).real / abs(dominant(p))

    ks = np.logspace(-4, np.log10(kmax), 4000)
    zs = np.array([z_of(k) for k in ks])
    i = np.argmax(zs <= zeta)
    if zs[i] > zeta:
        raise ValueError('damping %.2f not reached' % zeta)
    lo, hi = ks[max(i - 1, 0)], ks[i]
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        lo, hi = (mid, hi) if z_of(mid) > zeta else (lo, mid)
    return hi, np.roots(np.polyadd(den, hi * num))


def margins(L, delay=0.0, w=None):
    """wc, PM (deg), wpc, GM (abs), delay margin (s) of L(s) e^{-s delay} from a dense frequency grid."""
    w = np.logspace(-4, 2, 20000) if w is None else w
    H = np.asarray(L(1j * w)).ravel() * np.exp(-1j * w * delay)
    mag, ph = np.abs(H), np.degrees(np.unwrap(np.angle(H)))
    out = dict(wc=np.nan, PM=np.nan, wpc=np.nan, GM=np.inf, DM=np.nan)
    i = np.where((mag[:-1] >= 1) & (mag[1:] < 1))[0]
    if len(i):
        i = i[0]
        out['wc'] = np.interp(0, -np.log(mag[i:i + 2]), w[i:i + 2])
        out['PM'] = 180 + np.interp(out['wc'], w[i:i + 2], ph[i:i + 2])
        out['DM'] = np.radians(out['PM']) / out['wc']
    j = np.where((ph[:-1] > -180) & (ph[1:] <= -180))[0]
    if len(j):
        j = j[0]
        out['wpc'] = np.interp(-180, ph[j:j + 2][::-1], w[j:j + 2][::-1])
        out['GM'] = 1 / np.interp(out['wpc'], w[j:j + 2], mag[j:j + 2])
    return out


def recovery(t, v, v_ref, t0, band=0.2):
    """Max speed drop after t0 and the time (s, from t0) after which |v - v_ref| stays within band."""
    m = t >= t0
    tt, e = t[m] - t0, v[m] - v_ref
    out = np.abs(e) > band
    last = np.where(out)[0]
    t_rec = 0.0 if not len(last) else (tt[last[-1] + 1] if last[-1] + 1 < len(tt) else np.nan)
    return -e.min(), t_rec
