# Author: Xu Chen, University of Washington
"""Hard disk drive track-following servo: a voice-coil motor (VCM) rotates the actuator arm about its pivot,
and the head at the end of the suspension follows a data track. Representative values for a 3.5-inch,
7200 rpm class drive (not a specific product).

    VCM current u (A) -> rigid body + lightly damped arm/suspension modes (modal sum) -> head position y (tracks)
    G(s) = Kp [ 1/s^2 + sum_i r_i / (s^2 + 2 zeta_i w_i s + w_i^2) ],   Kp = Kt l_arm / (J TP)

The servo reads the position error signal (PES) from the servo bursts once per servo sector,
Ts = 60/(rpm * sectors). The current computed from PES(k) reaches the coil tau seconds later (computation
delay) and is held until the next update (zero-order hold); at low frequency this looks like a delay of
Ts/2 + tau. PES is logged in % of the track pitch TP (100 % = one track)."""
import numpy as np
import control as ct
from scipy import optimize, signal
from scipy.linalg import expm

# demo drive D1: 7200 rpm (120 Hz), 300 servo sectors per revolution (36 kHz), 60 nm track pitch
DRIVE_A = dict(
    name='D1', rpm=7200, sectors=300, track_pitch_nm=60.0, tau_s=12e-6,
    Kt=0.11, l_arm=0.048, J=1.9e-5,             # VCM torque constant (N m/A), pivot-to-head (m), inertia (kg m^2)
    modes=[dict(name='pivot-bearing mode', f=3300.0, zeta=0.050, r=0.08),
           dict(name='butterfly mode', f=6200.0, zeta=0.015, r=0.60),
           dict(name='suspension sway mode', f=11500.0, zeta=0.010, r=-0.25)],
    rro_track=[1.20, 0.12, 0.05, 0.03, 0.02],   # runout at 1x..5x the spindle frequency (after runout feedforward)
    rro_wi_track=0.003,                         # written-in runout per harmonic at 6x, falling as 1/sqrt(k) to 60x
    flutter=[(900.0, 0.03, 0.012), (1500.0, 0.03, 0.010)],   # disk flutter: f (Hz), zeta, rms (track)
    windage_A=0.15e-3, windage_fc=3000.0,       # airflow torque on the arm: current-equivalent rms (A), bandwidth
    noise_track=0.012,                          # PES demodulation noise, rms (track)
    pes_res_track=1 / 256,                      # PES resolution (track)
    bias_A=0.012,                               # flex-cable bias force, current-equivalent (A)
    i_lim_A=(-1.0, 1.0),
)


def sample_time(dr):
    return 60.0 / (dr['rpm'] * dr['sectors'])


def rigid_gain(dr):
    """Kp = Kt l_arm / (J TP) in tracks/(A s^2), or dr['Kp'] when given directly (identified models)."""
    if 'Kp' in dr:
        return dr['Kp']
    return dr['Kt'] * dr['l_arm'] / dr['J'] / (dr['track_pitch_nm'] * 1e-9)


def plant_ct(dr):
    """Continuous modal state space: input VCM current (A), output head position (tracks)."""
    Kp = rigid_gain(dr)
    n = 2 + 2 * len(dr['modes'])
    A, B, C = np.zeros((n, n)), np.zeros((n, 1)), np.zeros((1, n))
    A[0, 1], B[1, 0], C[0, 0] = 1.0, Kp, 1.0
    for i, m in enumerate(dr['modes']):
        w = 2 * np.pi * m['f']
        k = 2 + 2 * i
        A[k, k + 1] = 1.0
        A[k + 1, k], A[k + 1, k + 1] = -w ** 2, -2 * m['zeta'] * w
        B[k + 1, 0], C[0, k] = Kp * m['r'], 1.0
    return ct.ss(A, B, C, 0)


def _zoh_parts(A, B, t):
    """e^{At} and int_0^t e^{As} ds B."""
    n, m = A.shape[0], B.shape[1]
    Mx = np.zeros((n + m, n + m))
    Mx[:n, :n], Mx[:n, n:] = A, B
    E = expm(Mx * t)
    return E[:n, :n], E[:n, n:]


def plant_dt(dr, tau=None):
    """Exact sampled plant with the zero-order hold and a computation delay tau >= 0 (default dr['tau_s']).
    Input: current computed at sample k; output: head position at sample k (tracks).
    States: plant x (first), the previous held current, then n_int whole-sample delay states."""
    Ts = sample_time(dr)
    tau = dr['tau_s'] if tau is None else tau
    n_int = int(np.floor(tau / Ts + 1e-9))
    tf_ = tau - n_int * Ts
    P = plant_ct(dr)
    A, B, C = (np.asarray(m) for m in (P.A, P.B, P.C))
    n = A.shape[0]
    Phi, _ = _zoh_parts(A, B, Ts)
    _, G1 = _zoh_parts(A, B, tf_)               # previous current acts during [0, tau_frac)
    E0, G0 = _zoh_parts(A, B, Ts - tf_)         # new current during [tau_frac, Ts)
    N = n + 1 + n_int
    Ad, Bd = np.zeros((N, N)), np.zeros((N, 1))
    Cv, Dv = np.zeros(N), 1.0                    # current reaching the hold at sample k: v = Cv X + Dv u
    if n_int:
        Cv[N - 1], Dv = 1.0, 0.0
        Bd[n + 1, 0] = 1.0                       # buffer: newest first
        for i in range(n + 2, N):
            Ad[i, i - 1] = 1.0
    Ad[:n, :n], Ad[:n, n] = Phi, (E0 @ G1)[:, 0]
    Ad[:n, :] += np.outer(G0[:, 0], Cv)
    Bd[:n, 0] += G0[:, 0] * Dv
    Ad[n, :] += Cv
    Bd[n, 0] += Dv
    Cd = np.zeros((1, N))
    Cd[0, :n] = C[0]
    return ct.ss(Ad, Bd, Cd, 0, Ts)


def frf(sys_d, f):
    """Frequency response of a discrete LTI system at f (Hz): sys(e^{j 2 pi f Ts}), batched."""
    f = np.atleast_1d(np.asarray(f, float))
    S = ct.ss(sys_d)
    A, B, C, D = (np.asarray(m, complex) for m in (S.A, S.B, S.C, S.D))
    if A.shape[0] == 0:
        return np.full(len(f), D[0, 0])
    z = np.exp(2j * np.pi * f * S.dt)
    X = np.linalg.solve(z[:, None, None] * np.eye(A.shape[0]) - A, np.broadcast_to(B, (len(f),) + B.shape))
    return (C @ X)[:, 0, 0] + D[0, 0]


# ---------------------------------------------------------------------------------------- controller
def lead_ct(K, fz, fp):
    """K (s/wz + 1)/(s/wp + 1), corner frequencies in Hz."""
    return ct.tf([K / (2 * np.pi * fz), K], [1 / (2 * np.pi * fp), 1])


def notch_ct(f0, depth_dB, zeta_p):
    """(s^2 + 2 zeta_z w0 s + w0^2)/(s^2 + 2 zeta_p w0 s + w0^2), zeta_z = zeta_p 10^(-depth/20): unit gain
    away from f0, -depth_dB at f0; zeta_p sets the width (-3 dB band about 2 zeta_p f0 for a deep notch)."""
    w = 2 * np.pi * f0
    zz = zeta_p * 10 ** (-depth_dB / 20)
    return ct.tf([1, 2 * zz * w, w ** 2], [1, 2 * zeta_p * w, w ** 2])


def peak_ct(f0, gain_dB, zeta_z):
    """Peak filter (an inverted notch): +gain_dB at f0, unit gain away from it; zeta_z sets the width."""
    w = 2 * np.pi * f0
    zp = zeta_z * 10 ** (-gain_dB / 20)
    return ct.tf([1, 2 * zeta_z * w, w ** 2], [1, 2 * zp * w, w ** 2])


def controller(Ts, K, fz, fp, fi=None, notch=None, peak=None):
    """Discrete track-following controller C(z) = lead x integral x notch x peak filter (A/track).
    Lead K (s/wz + 1)/(s/wp + 1); integral (s + wi)/s with fi in Hz; notch = (f0, depth_dB, zeta_p);
    peak = (f0, gain_dB, zeta_z). Each block is discretized with the Tustin rule prewarped at its own key
    frequency (lead and integral at sqrt(fz fp), notch and peak at f0), as servo firmware does.
    Returns C(z) and the list of discrete blocks."""
    wc = 2 * np.pi * np.sqrt(fz * fp)
    blocks = [ct.c2d(lead_ct(K, fz, fp), Ts, 'tustin', prewarp_frequency=wc)]
    if fi:
        blocks.append(ct.c2d(ct.tf([1, 2 * np.pi * fi], [1, 0]), Ts, 'tustin', prewarp_frequency=wc))
    if notch:
        blocks.append(ct.c2d(notch_ct(*notch), Ts, 'tustin', prewarp_frequency=2 * np.pi * notch[0]))
    if peak:
        blocks.append(ct.c2d(peak_ct(*peak), Ts, 'tustin', prewarp_frequency=2 * np.pi * peak[0]))
    C = blocks[0]
    for b in blocks[1:]:
        C = C * b
    return C, blocks


def design_lead(Kp, Td, fc, PM, Ts, fi=None, notch=None, peak=None):
    """Lead for the loop crossover fc (Hz) and phase margin PM (deg) on the rigid-body design model
    Kp/s^2 e^{-s Td}. The integral, notch and peak blocks are chosen first and their phase at fc is
    included; the lead supplies the rest: phi = PM - 180 - angle(G C_rest), alpha = (1 - sin phi)/(1 + sin phi),
    fz = fc sqrt(alpha), fp = fc/sqrt(alpha), and K makes |L| = 1 at fc."""
    wc = 2 * np.pi * fc
    _, bl = controller(Ts, 1.0, fc / 2, fc * 2, fi, notch, peak)
    rest = np.prod([frf(b, fc)[0] for b in bl[1:]]) if len(bl) > 1 else 1.0 + 0j
    G = -Kp / wc ** 2 * np.exp(-1j * wc * Td) * rest
    phi = PM - 180 - (np.degrees(np.angle(G)) - 360 * (np.angle(G) > 0))
    s = np.sin(np.radians(phi))
    alpha = (1 - s) / (1 + s)
    fz, fp = fc * np.sqrt(alpha), fc / np.sqrt(alpha)
    C1, _ = controller(Ts, 1.0, fz, fp, fi, notch, peak)
    K = 1 / abs(frf(C1, fc)[0] * Kp / wc ** 2)
    return dict(K=K, fz=fz, fp=fp, phi=phi, fi=fi, notch=notch, peak=peak)


# ---------------------------------------------------------------------------------------- margins
def frf_margins(f, L):
    """Margins from samples of L(e^{jwTs}) at f (Hz), measured or model.
    fc, PM: every gain crossover and its phase margin; DM = PM/wc: the delay margin (s).
    f180, GM: every crossing of the negative real axis and 1/|L| there. Crossings below the first gain
    crossover come from the integral action and are lower limits (the gain may drop by that factor);
    GM_up = (f, GM) of the smallest GM among the crossings above it (GM < 1, i.e. negative in dB, means the
    loop gain already exceeds 1 there). Ms = max |1/(1 + L)| at f_Ms."""
    f, L = np.asarray(f, float), np.asarray(L, complex)
    mag = np.abs(L)
    out = dict(fc=[], PM=[], DM=[], f180=[], GM=[])
    for k in np.where(np.diff(np.sign(mag - 1)) != 0)[0]:
        a = (1 - mag[k]) / (mag[k + 1] - mag[k])
        Lx = L[k] + a * (L[k + 1] - L[k])
        fx = f[k] + a * (f[k + 1] - f[k])
        pm = np.degrees(np.angle(Lx)) % 360 - 180
        out['fc'].append(fx)
        out['PM'].append(pm)
        out['DM'].append(np.radians(pm) / (2 * np.pi * fx))
    im = L.imag
    for k in np.where(np.diff(np.sign(im)) != 0)[0]:
        a = im[k] / (im[k] - im[k + 1])
        Lx = L[k] + a * (L[k + 1] - L[k])
        if Lx.real < 0:
            out['f180'].append(f[k] + a * (f[k + 1] - f[k]))
            out['GM'].append(1 / abs(Lx))
    S = 1 / (1 + L)
    k = int(np.argmax(np.abs(S)))
    out['Ms'], out['f_Ms'] = float(np.abs(S[k])), float(f[k])
    f1 = out['fc'][0] if out['fc'] else 0.0
    up = [(fx, g) for fx, g in zip(out['f180'], out['GM']) if fx > f1]
    out['GM_up'] = min(up, key=lambda x: x[1]) if up else (np.nan, np.inf)
    return out


def cl_poles(P_d, C_d):
    """Closed-loop poles of the sampled loop with negative feedback; stable if all |z| < 1."""
    return ct.poles(ct.feedback(ct.series(ct.ss(C_d), P_d), 1))


# ---------------------------------------------------------------------------------------- data, model fit
def load_frf(name):
    """Swept-sine analyzer export: f (Hz), complex G (track/A) from magnitude and wrapped phase, coherence."""
    from ctrlprac import data
    df, meta = data.load(name)
    G = 10 ** (df['mag_dB'].values / 20) * np.exp(1j * np.radians(df['phase_deg'].values))
    return df['f_Hz'].values, G, df['coherence'].values, meta


def fitted_drive(Kp, Td, modes, Ts):
    """Drive dict for plant_dt from identified values: Kp, total delay Td = Ts/2 + tau, modes."""
    return dict(rpm=60.0 / Ts, sectors=1, Kp=Kp, tau_s=max(Td - Ts / 2, 0.0), modes=modes)


def fit_modal(f, G, coh, Ts, Kp0, Td0, modes0, coh_min=0.8):
    """Least-squares fit of Kp, Td and every mode (f_i, zeta_i, r_i) of the sampled modal model to the
    measured response, on the points with coherence >= coh_min, starting from values read off the Bode
    plot. Residual: log(G_model/G_meas) (magnitude in nepers, phase in rad), weighted by sqrt(coherence)."""
    k = coh >= coh_min
    fk, Gk, wk = f[k], G[k], np.sqrt(coh[k])

    def unpack(p):
        modes = [dict(name=m0.get('name', ''), f=p[2 + 3 * i], zeta=np.exp(p[3 + 3 * i]), r=p[4 + 3 * i])
                 for i, m0 in enumerate(modes0)]
        return np.exp(p[0]), p[1] * 1e-6, modes

    def res(p):
        Kp, Td, modes = unpack(p)
        e = np.log(frf(plant_dt(fitted_drive(Kp, Td, modes, Ts)), fk) / Gk)
        return np.concatenate([wk * e.real, wk * ((e.imag + np.pi) % (2 * np.pi) - np.pi)])

    p0 = [np.log(Kp0), Td0 * 1e6] + sum([[m['f'], np.log(m['zeta']), m['r']] for m in modes0], [])
    sol = optimize.least_squares(res, p0, x_scale='jac')
    Kp, Td, modes = unpack(sol.x)
    e = res(sol.x)
    n = len(e) // 2
    return dict(Kp=Kp, Td=Td, modes=modes, n_points=int(k.sum()),
                rms_dB=20 / np.log(10) * np.sqrt(np.mean((e[:n] / wk) ** 2)),
                rms_deg=np.degrees(np.sqrt(np.mean((e[n:] / wk) ** 2))))


def identify(f, G, coh, Ts):
    """The identification of demo 1 in one call: delay from the phase slope (line through -180 deg on
    150-1500 Hz), Kp from |G| w^2 on 150-600 Hz, modes from the peaks of G/rigid line between 2 and 12 kHz
    (half-power zeta, r from the antiresonance), then fit_modal. Returns the fit_modal dict."""
    good = coh >= 0.9
    ph = np.degrees(np.unwrap(np.angle(G)))
    ph -= 360 * np.round((ph[0] + 180) / 360)
    lo = good & (coh >= 0.95) & (f >= 150) & (f <= 1500)
    Td = -np.sum(f[lo] * (ph[lo] + 180)) / (360 * np.sum(f[lo] ** 2))
    lo2 = lo & (f <= 600)
    Kp = np.median(abs(G[lo2]) * (2 * np.pi * f[lo2]) ** 2)
    R = G / (-Kp / (2 * np.pi * f) ** 2 * np.exp(-2j * np.pi * f * Td))
    RdB = 20 * np.log10(abs(R))
    band = np.flatnonzero((f > 2000) & (f < 12000))            # peaks standing 3 dB above their surroundings
    pk = band[signal.find_peaks(RdB[band], prominence=3)[0]]
    pk = pk[coh[pk] >= 0.9]
    modes0 = []
    for i, k in enumerate(pk):
        half = RdB[k] - 3
        k1, k2 = k, k
        while k1 > 0 and RdB[k1] > half:
            k1 -= 1
        while k2 < len(f) - 1 and RdB[k2] > half:
            k2 += 1
        zeta = (f[k2] - f[k1]) / (2 * f[k]) if RdB[k] > 10 else 0.05
        near = np.flatnonzero((f > 0.75 * f[k]) & (f < 1.3 * f[k]) & (coh >= 0.3))
        kz = near[np.argmin(RdB[near])]
        modes0.append(dict(name='mode %d' % (i + 1), f=f[k], zeta=zeta, r=(f[k] / f[kz]) ** 2 - 1))
    return fit_modal(f, G, coh, Ts, Kp, Td, modes0)


# ---------------------------------------------------------------------------------------- disturbances, PES
def disturbances(dr, n, seed):
    """Disturbance sequences for n sectors. Output side (track runout, tracks): repeatable runout (RRO) at
    the spindle harmonics with phases fixed to the disk, and disk flutter (narrow-band random). Input side
    (A): windage torque on the arm (low-pass random) plus the flex-cable bias. PES noise (tracks)."""
    Ts = sample_time(dr)
    rng = np.random.default_rng(seed)
    t = np.arange(n) * Ts
    f_rot = dr['rpm'] / 60.0
    ph = np.random.default_rng(seed + 1000).uniform(0, 2 * np.pi, 60)
    n_low = len(dr['rro_track'])
    amps = list(dr['rro_track']) + [dr['rro_wi_track'] * np.sqrt((n_low + 1) / k) for k in range(n_low + 1, 61)]
    rro = sum(a * np.cos(2 * np.pi * (k + 1) * f_rot * t + ph[k]) for k, a in enumerate(amps))
    flut = np.zeros(n)
    for f0, z, rms in dr['flutter']:
        w = 2 * np.pi * f0
        F = ct.c2d(ct.tf([2 * z * w, 0], [1, 2 * z * w, w ** 2]), Ts, 'tustin', prewarp_frequency=w)
        x = signal.lfilter(np.ravel(F.num[0][0]), np.ravel(F.den[0][0]), rng.standard_normal(n))
        flut += rms * x / np.std(x)
    wf = 2 * np.pi * dr['windage_fc']
    Fw = ct.c2d(ct.tf([wf], [1, wf]), Ts, 'tustin')
    wnd = signal.lfilter(np.ravel(Fw.num[0][0]), np.ravel(Fw.den[0][0]), rng.standard_normal(n))
    wnd = dr['windage_A'] * wnd / np.std(wnd)
    noise = dr['noise_track'] * rng.standard_normal(n)
    return dict(t=t, rro=rro, flutter=flut, d_out=rro + flut, d_in=wnd + dr['bias_A'], noise=noise,
                rro_amps=np.array(amps))


def simulate_pes(dr, C_d, n_rev, seed, quant=True, tau=None, n_warm=5):
    """Track following of drive dr with the sampled controller C_d (A/track): n_rev revolutions after n_warm
    revolutions of warm-up (the drive was already following the track when the log starts). Per sector k:
    PES(k) = runout(k) - head(k) + noise(k) (quantized), u(k) = C(z) PES(k), applied tau later and held.
    Returns t (s), pes (track), u (A) and the disturbance sequences of the logged part."""
    n = (n_rev + n_warm) * dr['sectors']
    d = disturbances(dr, n, seed)
    P = plant_dt(dr, tau)
    A, B, C = (np.asarray(m) for m in (P.A, P.B, P.C))
    Pc = plant_ct(dr)
    np_ = Pc.nstates
    _, Gw = _zoh_parts(np.asarray(Pc.A), np.asarray(Pc.B), sample_time(dr))
    Bw = np.zeros(A.shape[0])                   # torque disturbance: held over the sample, no computation delay
    Bw[:np_] = Gw[:, 0]
    K = ct.ss(C_d)
    Ak, Bk, Ck, Dk = (np.asarray(m) for m in (K.A, K.B, K.C, K.D))
    x, xk = np.zeros(A.shape[0]), np.zeros(Ak.shape[0])
    pes, u_out = np.zeros(n), np.zeros(n)
    lo, hi = dr['i_lim_A']
    q = dr['pes_res_track']
    for k in range(n):
        e = d['d_out'][k] - (C @ x)[0] + d['noise'][k]
        if quant:
            e = q * np.round(e / q)
        u = min(max((Ck @ xk)[0] + Dk[0, 0] * e, lo), hi)
        xk = Ak @ xk + Bk[:, 0] * e
        x = A @ x + B[:, 0] * u + Bw * d['d_in'][k]
        pes[k], u_out[k] = e, u
    k0 = n_warm * dr['sectors']
    d = {key: (v[k0:] if isinstance(v, np.ndarray) and len(v) == n else v) for key, v in d.items()}
    t = np.arange(n - k0) * sample_time(dr)
    return dict(t=t, pes=pes[k0:], u=u_out[k0:], d=d)


def pes_stats(pes, sectors):
    """Synchronous average over revolutions: RRO(sector) = mean over the revolutions, NRRO = PES - RRO.
    Returns 3 sigma of PES, RRO and NRRO (same unit as pes) and the waveforms."""
    n_rev = len(pes) // sectors
    P = np.asarray(pes[:n_rev * sectors], float).reshape(n_rev, sectors)
    rro = P.mean(axis=0)
    nrro = P - rro
    return dict(three_sigma=3 * np.std(P), rro3=3 * np.std(rro), nrro3=3 * np.std(nrro),
                rro=rro, nrro=nrro.ravel(), n_rev=n_rev)


def rro_harmonics(rro):
    """Amplitudes of the spindle harmonics 1x, 2x, ... in one revolution of RRO (same unit)."""
    X = np.fft.rfft(rro) / len(rro)
    return 2 * np.abs(X[1:])


def predict_pes(pes, Ts, S_old, S_new):
    """PES under a new controller, predicted from a log under the old one for a linear loop with the same
    disturbances: PES_new = (S_new/S_old) PES_old on the FFT grid of the log. S_old, S_new: callables
    f (Hz) -> S(e^{j 2 pi f Ts}). The FFT treats the log as periodic, so use whole revolutions."""
    n = len(pes)
    f = np.fft.rfftfreq(n, Ts)
    X = np.fft.rfft(pes - np.mean(pes))
    X[1:] *= S_new(f[1:]) / S_old(f[1:])
    return np.fft.irfft(X, n) + np.mean(pes)
