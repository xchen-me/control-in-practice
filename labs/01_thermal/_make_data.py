# Author: Xu Chen, University of Washington
"""Write the simulated heater-board demo log data/01_thermal_step (board A). Seeded and deterministic.
(The homework logs are generated separately and their sidecars carry no model parameters.)"""
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ctrlprac import data, sim
from _model import BOARD_A, TS, QUANT, NOISE, U_LIM, true_plant


def meta(board, desc, test, seed):
    return dict(source='simulated', description=desc,
                true_model='T = amb + K e^{-delay s} / ((tau1 s + 1)(tau2 s + 1)) Q',
                parameters={k: v for k, v in board.items() if k != 'name'},
                units=dict(K='degC/%', tau1='s', tau2='s', delay='s', amb='degC'),
                test=test, sample_time_s=TS, sensor_noise_std_degC=NOISE, sensor_resolution_degC=QUANT,
                heater_limits_pct=list(U_LIM), seed=seed)


def step_log(board, u_before, u_after, t_step, T, seed, preroll=0.0):
    """Open-loop heater step; preroll seconds at u_before are simulated and discarded (steady start)."""
    n_pre, n = int(preroll / TS), int(T / TS) + 1
    k = np.arange(n_pre + n) * TS - preroll
    u = np.where(k < t_step, u_before, u_after)
    _, ym, _ = sim.open_loop(true_plant(board), TS, u, noise_std=NOISE, quant=QUANT,
                             delay=board['delay'], y0=board['amb'], seed=seed)
    return pd.DataFrame({'time_s': k[n_pre:], 'heater_pct': u[n_pre:], 'T_C': np.round(ym[n_pre:], 2)})


if __name__ == '__main__':
    # demo: board A starts at ambient, heater 0 -> 50 % at t = 60 s
    df = step_log(BOARD_A, 0.0, 50.0, 60.0, 1200.0, seed=101)
    data.save(df, '01_thermal_step', meta(BOARD_A, 'Heater step test on board A (demo), from ambient.',
                                          dict(heater_pct=[0, 50], t_step_s=60), 101))
    print('wrote data/01_thermal_step')
