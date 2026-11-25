# Author: Xu Chen, University of Washington
# %% [markdown]
# # Homework starter, Lab 06: belt stage S2
#
# ![scene](figures/scene_hw.png)
#
# Loads the two S2 logs and plots them. Fill in the TODO steps (see homework.md).

# %%
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import control as ct
HERE = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
sys.path.insert(0, HERE)
sys.modules.pop('_model', None)
from ctrlprac import data, ident, sim
import _model as M          # M.TS, M.I_LIM, M.ENC_RES_MM, ... (drive constants)

ch, mch = data.load('06_hw_chirp')
ss, mss = data.load('06_hw_steppedsine')
Ts = mch['sample_time_s']
Td = (mch['delay_samples'] + 0.5) * Ts
print('Kt = %.3f N m/A, r = %.4f m, delay + ZOH = %.3f ms'
      % (mch['parameters']['Kt'], mch['parameters']['r'], 1e3 * Td))

fig, ax = plt.subplots(3, 1, figsize=(9, 7.5))
ax[0].plot(ch['time_s'], ch['i_cmd_A'], lw=0.3, label='i_cmd')
ax[0].plot(ch['time_s'], ch['i_exc_A'], lw=0.3, alpha=0.6, label='i_exc (chirp)')
ax[0].set(xlabel='time (s)', ylabel='current (A)', title='06_hw_chirp')
ax[1].plot(ch['time_s'], ch['x_mm'], lw=0.5, color='C2', label='x (encoder)')
ax[1].set(xlabel='time (s)', ylabel='x (mm)')
ax[2].plot(ss['time_s'], ss['x_mm'], lw=0.4, color='C2', label='x (mm)')
ax[2].plot(ss['time_s'], ss['measure'], 'k', lw=0.6, label='measure flag')
ax[2].set(xlabel='time (s)', ylabel='x (mm) / flag', title='06_hw_steppedsine')
for a in ax:
    a.grid(True, alpha=0.4)
    a.legend(fontsize=8, loc='upper right')
fig.tight_layout()
plt.show()

# %% [markdown]
# ## 1. ETFE with the excitation as reference, and coherence
# TODO: w, Hdx, coh = ident.etfe(t, i_exc, x, nperseg=...); Hdu likewise; G = Hdx / Hdu

# %% [markdown]
# ## 2. Stepped sine
# TODO: for each segment, ident.sine_response on rows with measure == 1 (exc -> x and exc -> i_cmd)

# %% [markdown]
# ## 3. Bode reading: K, J, f_r, Q, zeta
# TODO

# %% [markdown]
# ## 4-5. fit_tf (unstructured and with the integrator), dominant-pole model
# TODO: ident.fit_tf(w, G, nb, na, delay=Td)

# %% [markdown]
# ## 6-7. Lead design, margins and Ms on the MEASURED response, step simulation
# TODO: L = C(1j w) G_measured; sim.simulate_loop(ct.ss(G_fit), sim.LTIController(C, Ts, u_lim=M.I_LIM), ...,
# delay=(mch['delay_samples'] - 1) * Ts) (simulate_loop applies one sample of delay by itself)
