# Author: Xu Chen, University of Washington
"""load_phyphox reads every CSV format phyphox exports (fixtures in tests/fixtures/)."""
import os
import sys
import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'labs', '04_pendulum'))
sys.modules.pop('_model', None)
from _model import load_phyphox  # noqa: E402

FORMATS = ['comma_point', 'tab_point', 'tab_comma', 'semicolon_point', 'semicolon_comma']


@pytest.mark.parametrize('fmt', FORMATS)
def test_load_phyphox_formats(fmt):
    rec = load_phyphox(os.path.join(HERE, 'fixtures', 'phyphox_%s.csv' % fmt))
    # the fixtures hold 50 samples at 100 Hz; z carries the swing 1.3 sin(2 pi 0.57 t) - 0.0025 rad/s
    t_raw = np.arange(50) / 100.0 + 0.0012
    assert rec['axis'] == 'z' and set(rec['columns']) == {'x', 'y', 'z'}
    assert abs(rec['fs'] - 100.0) < 1e-6
    assert np.allclose(rec['raw'].iloc[:, 0].to_numpy(float), t_raw, atol=1e-9)
    assert np.allclose(rec['w'], 1.3 * np.sin(2 * np.pi * 0.57 * (rec['t'] + 0.0012)) - 0.0025, atol=1e-6)
