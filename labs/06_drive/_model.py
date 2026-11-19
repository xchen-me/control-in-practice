# Author: Xu Chen, University of Washington
"""Belt-driven linear stage: servo motor + pulley, elastic timing belt, carriage on a linear guide,
linear encoder on the carriage (non-collocated). Representative values (not a specific product).

    current command i (A) -> current loop (1st order, 1 kHz) -> torque Kt i -> motor inertia J_m
    -> belt (torsional stiffness k_t = k r^2, damping c_t) -> carriage m_L (J_L = m_L r^2) -> x = r theta_L

Low frequency: rigid body x/i = Kt r / (s (J s + b)), J = J_m + J_L (-40 dB/dec above b/J).
Two-mass resonance at w_r = sqrt(k_t (1/J_m + 1/J_L)), lightly damped."""
import numpy as np
import control as ct

FS = 4000.0
TS = 1.0 / FS
N_DELAY = 2                     # samples of computation + bus delay (0.5 ms); the ZOH adds Ts/2
I_LIM = (-3.0, 3.0)             # amplifier current limit, A
ENC_RES_MM = 1e-3               # 1 um linear encoder
ENC_NOISE_MM = 0.3e-3
I_NOISE_A = 0.01                # force ripple / friction noise referred to current, A rms

# demo stage S1: 1.2 kg carriage, 1.0 m belt span, resonance near 60 Hz
STAGE_A = dict(name='S1', J_m=2.0e-5, r=0.012, m_L=1.2, k_belt=17.7e3, zeta_belt=0.04, b_L=8.0,
               Kt=0.12, f_cl=1000.0)


def derived(st):
    J_L = st['m_L'] * st['r'] ** 2
    kt = st['k_belt'] * st['r'] ** 2
    wr = np.sqrt(kt * (1 / st['J_m'] + 1 / J_L))
    Jr = 1 / (1 / st['J_m'] + 1 / J_L)                  # reduced inertia of the two-mass mode
    ct_ = 2 * st['zeta_belt'] * wr * Jr                  # belt damping giving zeta_belt
    J = st['J_m'] + J_L
    bt = st['b_L'] * st['r'] ** 2
    return dict(J_L=J_L, k_t=kt, c_t=ct_, w_r=wr, f_r=wr / 2 / np.pi, J=J, b_t=bt,
                K_rigid=st['Kt'] * st['r'] * 1e3 / J, a_rigid=bt / J)   # x(mm)/i(A) = K/(s(s + a))


def plant_ss(st, current_loop=True):
    """State space x = [th_m, w_m, th_L, w_L, (i_act)], input current command (A), output x (mm)."""
    d = derived(st)
    Jm, JL, k, c, b = st['J_m'], d['J_L'], d['k_t'], d['c_t'], d['b_t']
    A = np.array([[0, 1, 0, 0],
                  [-k / Jm, -c / Jm, k / Jm, c / Jm],
                  [0, 0, 0, 1],
                  [k / JL, c / JL, -k / JL, -(c + b) / JL]])
    B = np.array([[0], [st['Kt'] / Jm], [0], [0]])
    C = np.array([[0, 0, st['r'] * 1e3, 0]])
    P = ct.ss(A, B, C, 0)
    if current_loop:
        wc = 2 * np.pi * st['f_cl']
        P = ct.series(ct.tf([wc], [1, wc]), P)
    return P


def plant_tf(st, current_loop=True):
    return ct.tf(plant_ss(st, current_loop))


def lead(K, wc, phi_deg):
    """Lead K (s/z + 1)/(s/p + 1) with its maximum phase phi at wc: p/z = (1 + sin phi)/(1 - sin phi)."""
    s = np.sin(np.radians(phi_deg))
    alpha = (1 - s) / (1 + s)
    z, p = wc * np.sqrt(alpha), wc / np.sqrt(alpha)
    return ct.tf([K / z, K], [1 / p, 1]), z, p


def frf_margins(w, L):
    """Margins from samples of L(jw) (measured or model): crossover(s), PM, all -180 deg crossings with
    their gain margins, and the sensitivity peak Ms = max |1/(1 + L)|."""
    mag, ph = np.abs(L), np.unwrap(np.angle(L))
    ph = ph - 2 * np.pi * np.round((ph[0] + np.pi / 2) / (2 * np.pi))     # start near -90 deg (type 1)
    out = dict(wc=[], PM=[], w180=[], GM=[])
    for k in np.where(np.diff(np.sign(mag - 1)) != 0)[0]:
        f = (1 - mag[k]) / (mag[k + 1] - mag[k])
        out['wc'].append(w[k] + f * (w[k + 1] - w[k]))
        out['PM'].append(180 + np.degrees(ph[k] + f * (ph[k + 1] - ph[k])))
    im = np.imag(L)
    for k in np.where(np.diff(np.sign(im)) != 0)[0]:          # crossings of the real axis ...
        f = im[k] / (im[k] - im[k + 1])
        Lx = L[k] + f * (L[k + 1] - L[k])
        if Lx.real < 0:                                         # ... on its negative side
            out['w180'].append(w[k] + f * (w[k + 1] - w[k]))
            out['GM'].append(1 / abs(Lx))
    S = 1 / (1 + L)
    out['Ms'] = np.max(np.abs(S))
    out['w_Ms'] = w[np.argmax(np.abs(S))]
    return out


def measured_frf(df, nperseg=8192):
    """Plant FRF from a closed-loop excitation log: G = H(exc -> x) / H(exc -> i_cmd) (ctrlprac.ident.etfe).
    Returns w (rad/s), G (complex, mm/A, includes the drive delay), coherence of exc -> x."""
    from ctrlprac import ident
    t, d, u, x = (df[k].values for k in ('time_s', 'i_exc_A', 'i_cmd_A', 'x_mm'))
    w, Hdx, coh = ident.etfe(t, d, x, nperseg=nperseg)
    _, Hdu, _ = ident.etfe(t, d, u, nperseg=nperseg)
    return w, Hdx / Hdu, coh


def structured_fit(w, G, delay, band=None):
    """Fit G(s) = b0 / (s (s + a)(s^2 + 2 zeta w_r s + w_r^2)): the integrator is known physics, so fit
    s G(jw) with ctrlprac.ident.fit_tf (nb = 0, na = 3) after removing the known delay. Returns the full
    model and its parameters, K = b0 / w_r^2 being the rigid-body gain of K/(s (s + a))."""
    from ctrlprac import ident
    F = ident.fit_tf(w, 1j * w * G, 0, 3, delay=delay)
    p = ct.poles(F)
    r = p[np.argmax(abs(p.imag))]
    a = -p[np.argmin(abs(p))].real
    wr = abs(r)
    b0 = F.num[0][0][-1]
    return dict(model=F * ct.tf([1], [1, 0]), a=a, w_r=wr, f_r=wr / 2 / np.pi, zeta=-r.real / wr,
                K=b0 / wr ** 2, rigid=ct.tf([b0 / wr ** 2], [1, a, 0]))


def design_lead(Gd, wc, PM, Td):
    """Lead C = K (s/z + 1)/(s/p + 1) for loop crossover wc and phase margin PM (deg) on the design
    model Gd with delay Td: the lead supplies phi = PM - 180 - angle(Gd(jwc) e^{-jwc Td}) at wc."""
    G = Gd(1j * wc) * np.exp(-1j * wc * Td)
    phi = PM - 180 - np.degrees(np.angle(G) - 2 * np.pi * (np.angle(G) > 0))
    C1, z, p = lead(1.0, wc, phi)
    K = 1 / abs(C1(1j * wc) * G)
    C, z, p = lead(K, wc, phi)
    return C, dict(K=K, z=z, p=p, phi=phi)
