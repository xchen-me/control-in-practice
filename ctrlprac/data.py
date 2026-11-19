# Author: Xu Chen, University of Washington
"""Datasets live in data/<name>.csv with a sidecar data/<name>.json that records provenance:
"source": "measured" (with who/where/how) or "simulated" (with the model and its parameters)."""
import json
import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, 'data')


def load(name):
    """Return (DataFrame, metadata dict) for data/<name>.csv; prints the provenance line."""
    df = pd.read_csv(os.path.join(DATA, name + '.csv'), comment='#')
    meta_path = os.path.join(DATA, name + '.json')
    meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {}
    print('%s: %d samples, %s' % (name, len(df), meta.get('source', 'provenance unknown')))
    return df, meta


def save(df, name, meta):
    """Write data/<name>.csv and its provenance sidecar. meta must contain 'source' and 'description'."""
    assert meta.get('source') in ('measured', 'simulated'), "meta['source'] must be measured or simulated"
    os.makedirs(DATA, exist_ok=True)
    df.to_csv(os.path.join(DATA, name + '.csv'), index=False, float_format='%.6g')
    with open(os.path.join(DATA, name + '.json'), 'w') as fh:
        json.dump(meta, fh, indent=2)
