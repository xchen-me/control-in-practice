# Author: Xu Chen, University of Washington
"""Write notebooks/<lab>/<script>.ipynb from every demo and starter (jupytext percent format).

usage:  python tools/make_notebooks.py           (the .py files stay the source; edit those)
"""
import glob, os
import jupytext
import nbformat

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for path in sorted(glob.glob(os.path.join(HERE, 'labs', '*', '*.py'))):
    if os.path.basename(path).startswith('_'):
        continue
    lab = os.path.basename(os.path.dirname(path))
    nb = jupytext.read(path, fmt='py:percent')
    # notebooks run from notebooks/<lab>/: put the lab folder and the repo root on the path first
    script = os.path.basename(path)
    setup = ("import os, sys\n"
             "__file__ = os.path.abspath('../../labs/%s/%s')   # the demos locate data from here\n"
             "sys.path[:0] = [os.path.abspath('../..'), os.path.dirname(__file__)]" % (lab, script))
    nb.cells.insert(0, nbformat.v4.new_code_cell(setup))
    for cell in nb.cells:                 # images are referenced relative to the lab folder
        if cell.cell_type == 'markdown':
            cell.source = cell.source.replace('](figures/', '](../../labs/%s/figures/' % lab)
    out = os.path.join(HERE, 'notebooks', lab, os.path.basename(path)[:-3] + '.ipynb')
    os.makedirs(os.path.dirname(out), exist_ok=True)
    jupytext.write(nb, out)
    print(os.path.relpath(out, HERE))
