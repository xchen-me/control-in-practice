# Author: Xu Chen, University of Washington
"""Write the "How to run" section at the end of every labs/<lab>/README.md (it replaces a final "## Run" or
"## How to run" section), so that all labs give the same step-by-step instructions.

usage:  python tools/write_run_section.py            (all labs)
        python tools/write_run_section.py 07_        (labs whose folder name contains '07_')
"""
import glob, os, re, sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

TEMPLATE = """## How to run

**Get the materials once.** Clone the repository or download it as a zip file (see
[Get the materials](../../README.md#get-the-materials)), and put the folder on a local disk, not inside Google
Drive, OneDrive or iCloud. Every command below is typed in a terminal (macOS: Terminal; Windows: PowerShell)
inside that folder. It is called `control-in-practice` after `git clone` and `control-in-practice-main` after
unzipping; `cd control-in-practice*` works for both.

**Option 1: JupyterLab in Docker (recommended).** Install [Docker Desktop](https://www.docker.com/products/docker-desktop/)
and start it, then run

```bash
cd control-in-practice*
docker compose up --build         # the first time builds the image (a few minutes); afterwards: docker compose up
```

When the terminal shows a line with `http://127.0.0.1:8888/lab`, open http://localhost:8888 in your browser. In
the file browser on the left, open `notebooks/{lab}/` and double-click a notebook. Run its cells in order with
Shift+Enter; the figures appear under the cells. The notebooks of this lab are

{table}

Your changes are saved in the folder on your disk. To stop, press Ctrl+C in the terminal and then run
`docker compose down`. To run all demos and the homework starter of this lab without the browser, use

```bash
docker compose run --rm lab python run_all.py {key}
```

The printed numbers appear in the terminal and the figures are saved in `out/` as `{lab}_<script>_<n>.png`.

**Option 2: Python on your own computer.** You need Python 3.11 or newer. Set up an environment once:

```bash
cd control-in-practice*
python -m venv .venv
source .venv/bin/activate         # Windows PowerShell: .venv\\Scripts\\Activate.ps1
python -m pip install -r requirements.txt
```

In every new terminal, activate the environment again and tell Python where `ctrlprac` is:

```bash
cd control-in-practice*
source .venv/bin/activate         # Windows PowerShell: .venv\\Scripts\\Activate.ps1
export PYTHONPATH=$PWD            # Windows PowerShell: $env:PYTHONPATH = $PWD
python labs/{lab}/{first}         # one demo; each figure opens in a window, close it to continue
python run_all.py {key}               # all demos and hw_starter.py of this lab, figures saved in out/
jupyter lab                       # the same notebooks as in Option 1, in your browser
```

**Homework.** Write your solution in `labs/{lab}/hw_starter.py` (run it with `python labs/{lab}/hw_starter.py`)
or in `notebooks/{lab}/hw_starter.ipynb`, and submit the file you worked in. The two start out identical but are
separate files, so pick one and stay with it.
"""


def title(path):
    m = re.search(r'^# #+ (.+)$', open(path).read(), re.M)
    return re.sub(r'^(Lab \d+, demo \d+|[A-Z][a-z]+ [a-z]+ \d+|Lab \d+ homework starter): ', '', m[1]) if m else ''


for d in sorted(glob.glob(os.path.join(HERE, 'labs', '*_*'))):
    lab = os.path.basename(d)
    if len(sys.argv) > 1 and sys.argv[1] not in lab:
        continue
    scripts = sorted(os.path.basename(p) for p in glob.glob(os.path.join(d, 'demo_*.py')))
    rows = ['| notebook | content |', '|---|---|']
    rows += ['| `%s.ipynb` | %s |' % (s[:-3], title(os.path.join(d, s))) for s in scripts]
    rows += ['| `hw_starter.ipynb` | starter for the homework |']
    text = TEMPLATE.format(lab=lab, key=lab[:3], first=scripts[0], table='\n'.join(rows))
    readme = os.path.join(d, 'README.md')
    src = open(readme).read()
    m = re.search(r'^## (How to run|Run)\s*$', src, re.M)
    src = (src[:m.start()] if m else src.rstrip() + '\n\n') + text
    open(readme, 'w').write(src)
    print(os.path.relpath(readme, HERE))
