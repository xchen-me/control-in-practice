# Lab 01: temperature control of a heater-and-sensor board

A transistor heats a small aluminium fin and a thermistor measures its temperature, in the style of the open
TCLab kit. We log a heater step, identify a first-order-plus-dead-time (FOPDT) model, see why proportional
control leaves an error, design a PI/PID controller on the root locus, run it with a saturating heater (and
show windup), and use Bode and Nyquist plots to see how the dead time caps the phase margin and bandwidth.

![scene](figures/scene.png)

## System and parameters

| quantity | demo board A | homework board B | basis |
|---|---|---|---|
| heater input | 0-100 % (PWM of the transistor) | same | TCLab-style firmware |
| static gain $K$ | 0.62 degC/% | identify from the log | representative values |
| heater/fin time constant $\tau_1$ | 150 s | identify from the log | representative (TCLab fits give a few minutes) |
| thermistor lag $\tau_2$ | 15 s | identify from the log | representative |
| transport delay | 4 s | about 6 s (board documentation) | representative |
| ambient | 21.0 degC | read from the log | |
| sampling, resolution, noise | 1 s, 0.1 degC, 0.04 degC rms | same | |

"True" plant used to simulate the logs: $T = T_a + \dfrac{K e^{-\theta_0 s}}{(\tau_1 s + 1)(\tau_2 s + 1)} Q$.
The logs in `data/01_*.csv` are **simulated** (the demo `.json` sidecar gives the exact model, noise and
seed; the homework sidecars give only the test conditions); they are not measurements from a real board.
The identification fits a simpler FOPDT model,
$G(s) = K e^{-\theta s}/(\tau s + 1)$, so there is a deliberate model mismatch to check.

## Tool map

| step | course tool | lecture |
|---|---|---|
| heater step test -> $K$, $\tau$, $\theta$ (`demo_1_identify.py`) | first-order step response, parameter reading | L6-L7 |
| P control, $e_{ss} = (r - T_a)/(1 + K K_p)$ (`demo_2_p_control.py`) | steady-state error, type-0 loop | L12 |
| Pade approximation, PI zero at $1/\tau$, gain for $\zeta = 0.7$ (`demo_3_pid_rootlocus.py`) | root locus, PI/PID design | L13-L17, L19 |
| heater saturation, windup and clamping anti-windup (`demo_3`) | PID implementation | L19 |
| controller on identified model vs true plant, fan disturbance (`demo_3`) | model mismatch, disturbance rejection | L12, L22 |
| PM, GM, delay margin, bandwidth vs delay (`demo_4_bode_nyquist.py`) | Bode, Nyquist, margins | L20-L22, L23-L26 |

Key demo results (board A): fitted $K$ = 0.620 degC/%, $\tau$ = 152 s, $\theta$ = 18.6 s (4 s transport +
1 s computation delay, since each logged heater value acts from the next sample, + 13.6 s for the sensor lag);
PI $K_p$ = 7.18 %/degC, $K_i$ = 0.047 %/(degC s); PM = 59 deg on the model, 62 deg on the true plant; without
anti-windup the 21 -> 50 degC step overshoots by about 17 %.

## Files

- `figures/scene.tex`, `figures/scene_hw.tex`: TikZ sources of the scene figures (build with `tools/build_figures.sh 01_`).
- `_make_data.py` writes the logs (run once; deterministic).
- `_model.py` lab helpers (true boards, root-locus gain search, margins with exact delay).
- `demo_1` ... `demo_4` classroom demos (jupytext percent format; open as notebooks or run as scripts).
- `homework.md`, `hw_starter.py` homework on board B (different unit and test than the demos).

## How to run

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
the file browser on the left, open `notebooks/01_thermal/` and double-click a notebook. Run its cells in order with
Shift+Enter; the figures appear under the cells. The notebooks of this lab are

| notebook | content |
|---|---|
| `demo_1_identify.ipynb` | a first-order-plus-dead-time model from a heater step test (L6-L7) |
| `demo_2_p_control.ipynb` | proportional control leaves a steady-state error (L12) |
| `demo_3_pid_rootlocus.ipynb` | PI/PID by root locus, saturation and anti-windup (L13-L19) |
| `demo_4_bode_nyquist.ipynb` | the delay limits phase margin and bandwidth (L20-L26) |
| `hw_starter.ipynb` | starter for the homework |

Your changes are saved in the folder on your disk. To stop, press Ctrl+C in the terminal and then run
`docker compose down`. To run all demos and the homework starter of this lab without the browser, use

```bash
docker compose run --rm lab python run_all.py 01_
```

The printed numbers appear in the terminal and the figures are saved in `out/` as `01_thermal_<script>_<n>.png`.

**Option 2: Python on your own computer.** You need Python 3.11 or newer. Set up an environment once:

```bash
cd control-in-practice*
python -m venv .venv
source .venv/bin/activate         # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

In every new terminal, activate the environment again and tell Python where `ctrlprac` is:

```bash
cd control-in-practice*
source .venv/bin/activate         # Windows PowerShell: .venv\Scripts\Activate.ps1
export PYTHONPATH=$PWD            # Windows PowerShell: $env:PYTHONPATH = $PWD
python labs/01_thermal/demo_1_identify.py         # one demo; each figure opens in a window, close it to continue
python run_all.py 01_               # all demos and hw_starter.py of this lab, figures saved in out/
jupyter lab                       # the same notebooks as in Option 1, in your browser
```

**Homework.** Write your solution in `labs/01_thermal/hw_starter.py` (run it with `python labs/01_thermal/hw_starter.py`)
or in `notebooks/01_thermal/hw_starter.ipynb`, and submit the file you worked in. The two start out identical but are
separate files, so pick one and stay with it.
