# Author: Xu Chen, University of Washington
"""Run every demo and homework starter headless; figures go to out/<lab>_<script>_<n>.png.

usage:  python run_all.py                 (everything under labs/)
        python run_all.py thermal         (paths containing 'thermal')
"""
import glob, os, runpy, sys, traceback
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out')
os.makedirs(OUT, exist_ok=True)
key = sys.argv[1] if len(sys.argv) > 1 else ''
scripts = [p for p in sorted(glob.glob(os.path.join(HERE, 'labs', '*', '*.py')))
           if key in p and not os.path.basename(p).startswith('_')]
failed = []
for path in scripts:
    name = os.path.basename(os.path.dirname(path)) + '_' + os.path.basename(path)[:-3]
    count = [0]

    def show(*args, **kwargs):          # replace plt.show(): save the figure instead
        count[0] += 1
        plt.gcf().savefig(os.path.join(OUT, '%s_%d.png' % (name, count[0])), dpi=110)
        plt.close('all')

    plt.show = show
    print('=' * 20, os.path.relpath(path, HERE))
    for mod in [m for m in sys.modules if m.startswith('_')]:   # lab helpers (_model) differ per lab
        if getattr(sys.modules[mod], '__file__', '') and HERE in sys.modules[mod].__file__:
            del sys.modules[mod]
    try:
        sys.argv = [path]
        os.chdir(os.path.dirname(path))
        runpy.run_path(path, run_name='__main__')
        print('-> ok, %d figure(s)' % count[0])
    except Exception:
        traceback.print_exc()
        failed.append(name)
    finally:
        os.chdir(HERE)
print('\nfailed: %s' % (failed or 'none'))
sys.exit(1 if failed else 0)
