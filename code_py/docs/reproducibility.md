# Validation and reproducibility

This page lists the scripts and inputs for the archived benchmarks. The complete map from each table and figure of the paper to its script is in `reproduce/README.md` of the archive.

## Error metrics (paper, Eq. 3)

| Symbol | Definition | Reference value |
|---|---|---|
| $\varepsilon_T$ | $\lvert T-T_{\rm ref}\rvert/T_{\rm ref}$ | closed-form elliptic period ($c=0$) or 40-digit quadrature |
| $e_x$ | $\max_j \lvert x(t_j)-x_{\rm ref}(t_j)\rvert$, 401 uniform points on $[0,t_w]$, $t_w=20$ (or $0.9\,t_{\rm esc}$ for escaping orbits) | 32/28-digit `NDSolve` (Wolfram) or 25-digit mpmath Taylor solution from the exact binary inputs (Python) |
| $\delta_E$ | $\max_j\lvert E_j-E_0\rvert/(1+\lvert E_0\rvert)$ | — |

## What to run

From the root of the archive, after `pip install ./code_py` (Python) or with `wolframscript` (Wolfram):

| Result in the paper | Command | Output |
|---|---|---|
| Wolfram tests (37) | `wolframscript -file code_wl/RunTests.wls` | `code_wl/test_report.*` |
| Wolfram figures and tables | `wolframscript -file code_wl/RunAll.wls` | `figures/`, `tables/` at the archive root (archived: `data/figures/`, `data/wolfram/`) |
| Python tests (Table 3) | `cd code_py && python -m pytest -q --cov=nonlinear_oscillator_lab --cov-branch` | terminal |
| Closed-form campaigns (Table 4) | `python reproduce/closed_form_campaign.py --n 1000` (about 15 min on 4 cores; resumable) | `reproduce/outputs/closed_form_campaign_py*` |
| Period baselines (Table 5) | `python reproduce/baseline_comparison.py` | `reproduce/outputs/baseline_*` |

Archived outputs of every command are in `data/` and are never overwritten. Timings depend on the machine; the reference platform of the paper is listed in
`data/results/baseline_summary.json` (`_meta`).

## Campaign design (Table 4)

* **R200** – the 200 problems of the published Wolfram campaign (`data/wolfram/closed_form_campaign.csv`,
  `SeedRandom[20261019]`, $a,b\sim U[-12,12]$, $c\sim U[0,1]$, $x_0\sim U[-2,0]$, $v_0\sim U[-2,2]$).
* **RWIDE** – NumPy `default_rng(20261103)`: $a,b\sim U[-12,12]$, $c\sim U[-1,1]$, $x_0,v_0\sim U[-2,2]$.
* **SEP** – $E=E_s\pm\lvert E_s\rvert 10^{-k}$ ($k=1..12$) for the double well, the triple well and a
  cubic–quintic double well, initial states in every well.
* **REP** – $b^2=4ac$ (inflection equilibria), states at and $10^{-8}$, $10^{-4}$ away from minima,
  vanishing coefficients.
* **UNB** – $c<0$ or $b<0=c$: $E=E_b(1\pm10^{-k})$, states beyond the barrier and far above it.

Reference: SciPy DOP853, `rtol=2.5e-14`, `atol=1e-16`; every case with $e_x>10^{-10}$ is recomputed
against a 25-digit mpmath Taylor-series solution started from the exact binary inputs (column `ref = mp25`;
`mp25bin` in the re-referenced archive file).
