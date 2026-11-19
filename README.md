# Control in Practice

**Xu Chen · Introduction to Control Theory**

Practical demos and homework that take the course tools off the whiteboard. Each lab starts from a logged
experiment on a real kind of machine and carries it through the whole toolchain: a physical model, a transfer
function identified from data, a root-locus or Bode design, a PID/lead/lag implementation that has to live with
saturation, noise, quantization and delay, and a Nyquist/margin check of how robust the result is.

| Lab | System | Data you start from | Main course tools | Lectures |
|---|---|---|---|---|
| [01 thermal](labs/01_thermal/README.md) | heater-and-sensor board (TCLab style) | heater step log, P-control log | FOPDT identification, steady-state error, PI on the root locus with a Pade delay, anti-windup, delay vs phase margin | 6-7, 12, 13-19, 20-26 |
| [02 cruise](labs/02_cruise/README.md) | car speed control | throttle step, hill-climb log | physics vs data, grade as a disturbance, P vs PI vs lag, root-locus PI zero, windup on a steep hill | 2-3, 6-7, 12-19 |
| [03 motor](labs/03_motor/README.md) | DC gearmotor with encoder | voltage step on speed | dominant-pole reduction, step identification, P/PD root locus, zero effect, lead vs PD on quantized data, lag/PI | 3, 7-9, 12-18 |
| [04 pendulum](labs/04_pendulum/README.md) | **your phone** as a pendulum (phyphox app) | **your own** gyroscope recording | log decrement, second-order fit, s-plane poles, settling time, nonlinear effects | 6-8, 16 |
| [05 balance](labs/05_balance/README.md) | two-wheeled self-balancing robot | drop test, push test, tip-over log | unstable pole, Routh range for PD, root locus, Bode of an unstable plant, Nyquist with P = 1, two-sided gain margin, delay limit | 3, 10-11, 13-17, 22-26 |
| [06 drive](labs/06_drive/README.md) | belt-driven linear stage with a resonance | PRBS/chirp and stepped-sine logs | measured frequency response, coherence, resonance from the Bode peak, transfer-function fit, lead on the Bode plot, margins and sensitivity on measured data | 20-22, 23-26 |

Each lab folder has a `README.md` (the system, its parameters and a tool map), worked `demo_*.py` scripts for
class, a `homework.md` sheet (100 points, uses a different unit or test than the demos) and a `hw_starter.py`.

**About the data.** Lab 04 uses your own measurement. The other logs are simulated from physical models with
realistic sensor noise, resolution, delays and actuator limits; each `data/*.json` file says so. Instructors
with the hardware (for example a TCLab board) can drop measured logs in with the same columns.

## Get the materials

The lab materials are in the GitHub repository
[github.com/xchen-me/control-in-practice](https://github.com/xchen-me/control-in-practice). Clone it with
git, or download it as a zip file without git:

```bash
git clone https://github.com/xchen-me/control-in-practice.git
```

Zip: [control-in-practice-main.zip](https://github.com/xchen-me/control-in-practice/archive/refs/heads/main.zip).
The zip unpacks into a folder named `control-in-practice-main`; the `cd control-in-practice*` lines below work
for either folder name. Put the folder on a local disk, not inside Google Drive, OneDrive or iCloud.

## Start: Docker (recommended, same environment for everyone)

Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) and start it. Make sure the folder is on
a local disk (Docker cannot mount cloud-synced folders). Then, in a terminal (macOS: Terminal; Windows: PowerShell):

```bash
cd control-in-practice*          # control-in-practice (git clone) or control-in-practice-main (zip)
docker compose up --build        # the first time builds the image (a few minutes); afterwards: docker compose up
```

When the terminal shows a line with `http://127.0.0.1:8888/lab`, open http://localhost:8888 in your browser. In
the file browser on the left, open `notebooks/`, then the folder of your lab (for example `02_cruise/`), and
double-click a demo notebook. Run its cells in order with Shift+Enter; the figures appear under the cells. Your
changes are saved in the folder on your disk. To stop, press Ctrl+C in the terminal and run `docker compose down`.

Without the browser:

```bash
docker compose run --rm lab python run_all.py 02_    # all demos and the homework starter of Lab 02
docker compose run --rm check                        # every lab and the tests; figures in out/
```

The notebooks in `notebooks/` are made from the `.py` files in `labs/` (`python tools/make_notebooks.py`); both
contain the same code.

## Start: without Docker

You need Python 3.11 or newer. Set up an environment once:

```bash
cd control-in-practice*                # the folder from the clone or the zip
python -m venv .venv
source .venv/bin/activate              # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

In every new terminal:

```bash
cd control-in-practice*
source .venv/bin/activate              # Windows PowerShell: .venv\Scripts\Activate.ps1
export PYTHONPATH=$PWD                 # Windows PowerShell: $env:PYTHONPATH = $PWD
python labs/01_thermal/demo_1_identify.py    # one demo; each figure opens in a window, close it to continue
python run_all.py 01_                  # all demos and hw_starter.py of Lab 01, figures saved in out/
jupyter lab                            # the notebooks in your browser
```

Each lab README ends with a "How to run" section that lists its notebooks and the commands for that lab.

## Shared tools (`ctrlprac/`)

- `ident.py`: `fit_fopdt`, `fit_second_order_step`, `log_decrement`, `etfe` (frequency response from broadband
  data, with coherence), `sine_response` (stepped sine), `fit_tf` (transfer function from a frequency response).
- `sim.py`: `PID` (derivative filter, anti-windup), `LTIController` (any C(s), Tustin), `simulate_loop`
  (sampled-data loop with saturation, noise, quantization and delay), `open_loop`.
- `data.py`: `load` and `save` for the logs and their provenance files.

Tests: `pytest -q tests`. Instructors: `tools/check_in_docker.sh` runs everything in the image from any folder.
