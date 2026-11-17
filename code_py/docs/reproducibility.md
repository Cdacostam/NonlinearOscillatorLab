# Validation and reproducibility

This page lists the scripts and inputs for the bundled benchmarks. The complete map from each table and figure of the paper to its script is in `reproduce/README.md` of the archive.

## Error metrics (paper, Eq. 2)

| Symbol | Definition | Reference value |
|---|---|---|
| $\varepsilon_T$ | $\lvert T-T_{\rm ref}\rvert/T_{\rm ref}$ | closed-form elliptic period ($c=0$) or 40-digit quadrature |
| $e_x$ | $\max_j \lvert x(t_j)-x_{\rm ref}(t_j)\rvert$, 401 uniform points on $[0,t_w]$, $t_w=20$ (or $0.9\,t_{\rm esc}$ for escaping orbits) | 32/28-digit `NDSolve` (Wolfram) or DOP853 with a conditional 25-digit Taylor replacement (Python) |
| $\delta_E$ | $\max_j\lvert E_j-E_0\rvert/(1+\lvert E_0\rvert)$ | — |

## What to run

Run from the bundle root. Install `./code_py` for the Python scripts, or install
`python -m pip install './code_py[test,bench]'` for pytest, coverage and figures.
Wolfram runs require a licensed kernel.

| Result | Command | Newly generated output | Recorded evidence |
|---|---|---|---|
| Python tests (Tables 3/S3) | `python -m pytest code_py/tests -q -o pythonpath= --cov=nonlinear_oscillator_lab --cov-branch --cov-report=json:reproduce/outputs/coverage.json` | `reproduce/outputs/coverage.json` | `data/results/v0.2.2_tests/` |
| Closed-form campaign (Table 4) | `python reproduce/closed_form_campaign.py --workers 4` | `reproduce/outputs/closed_form_campaign_py*` | `data/results/v0.2.2_campaign/` |
| Period baselines (Table 5) | `python reproduce/baseline_comparison.py` | `reproduce/outputs/baseline_*` | `data/results/v0.2.2_period/` |
| Figure 3/S5 | `python reproduce/period_accuracy_figure.py` | `reproduce/outputs/period_accuracy_panel_*.csv` and `fig6_period_accuracy.*` | Python series in `data/results/v0.2.2_period/` and Wolfram series in `data/wolfram/period_accuracy.csv` |
| Nine example periods (Table S5) | `python reproduce/python_benchmarks.py` | `reproduce/outputs/python_periods.csv` | `data/results/v0.2.2_period/python_periods.csv` and `data/wolfram/benchmark_periods.csv` |
| Period near degeneracy | `python reproduce/period_boundary_check.py` | `reproduce/outputs/period_boundary_check.json` | `data/results/v0.2.2_boundary/` |
| Wolfram tests | `wolframscript -file code_wl/RunTests.wls` | `code_wl/test_report.*` | `code_wl/verification/` |

Timing depends on the platform documented in the recorded summary. The Python columns in
Table S5 use the same binary64 inputs for the double and 40-digit estimates. The Wolfram
columns use the separately archived comparison with a rational-input reference.

## Campaign design (Table 4)

* **R200** – the 200 problems of the published Wolfram campaign (`data/wolfram/closed_form_campaign.csv`,
  `SeedRandom[20261019]`, $a,b\sim U[-12,12]$, $c\sim U[0,1]$, $x_0\sim U[-2,0]$, $v_0\sim U[-2,2]$).
* **RWIDE** – NumPy `default_rng(20261103)`: $a,b\sim U[-12,12]$, $c\sim U[-1,1]$, $x_0,v_0\sim U[-2,2]$.
* **SEP** – $E=E_s\pm\lvert E_s\rvert 10^{-k}$ ($k=1..12$) for the double well, the triple well and a
  cubic–quintic double well, with absolute offsets of $10^{-k}$ when $E_s=0$.
* **REP** – $b^2=4ac$ (inflection equilibria), states at and $10^{-8}$, $10^{-4}$ away from minima,
  vanishing coefficients.
* **UNB** – $c<0$ or $b<0=c$: $E=E_b(1\pm10^{-k})$, states beyond the barrier and far above it.

Reference: SciPy DOP853, `rtol=2.5e-14`, `atol=1e-16`; every case with $e_x>10^{-10}$ is recomputed
against a 25-digit mpmath Taylor-series solution started from the exact binary inputs (column `ref = mp25`).


The fixed binary64 trajectory inputs for Table 4 are in `data/inputs/closed_form_campaign.csv`.
Default campaign execution reads those values. `--generate` exercises the distributions with a
new sample. Reported error statistics describe discrepancies against the corresponding row
references. The 30 fixed cases marked in the input CSV compare DOP853 with a 32-digit Taylor solution.
The recorded Table 4 data are in `data/results/v0.2.2_campaign/`.
