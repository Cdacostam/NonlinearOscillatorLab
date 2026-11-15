# Reproduce the tables and figures

Run from the archive root after `python -m pip install ./code_py` (Python) or with `wolframscript`
(Wolfram). The figure script `period_accuracy_figure.py` also needs `matplotlib`. Python outputs go to
`reproduce/outputs/` (or to the folder in the environment variable `NOL_OUTPUTS`); the Wolfram script
`code_wl/MakeFigures.wls` writes to `figures/` and `tables/` at the archive root. Nothing in `data/` is
overwritten: compare the new outputs with the archived files listed below.

## Map: table or figure → command → archived file

"Fig." and "Table" refer to the paper, "S" to the supplementary material.

| Paper | Supplement | Command | Archived file |
|---|---|---|---|
| Fig. 1 (potentials) | Fig. S1 | `wolframscript -file code_wl/MakeFigures.wls` | `data/figures/fig1_potentials.pdf` |
| — | Fig. S2 (phase portraits) | `wolframscript -file code_wl/MakeFigures.wls` | `data/figures/fig2_phase_portraits.pdf` |
| Fig. 2 (closed form vs numerics) | Fig. S3 | `wolframscript -file code_wl/MakeFigures.wls` | `data/figures/fig3_closed_form_vs_numerics.pdf` |
| Fig. 3 (period accuracy) | Fig. S5 | `python reproduce/period_accuracy_figure.py` (panel (a) from `data/wolfram/period_accuracy.csv`; panel (b) computed with `period_mp` at 40 digits) | `data/figures/fig6_period_accuracy.pdf`, `data/results/period_accuracy_panel_b.csv` |
| Fig. 4 (Wolfram campaign, 200 problems) | Fig. S6 | `wolframscript -file code_wl/MakeFigures.wls` (`SeedRandom[20261019]`) | `data/figures/fig4_closed_form_campaign.pdf`, `data/wolfram/closed_form_campaign.csv`, `data/wolfram/closed_form_campaign_summary.txt` |
| — | Fig. S4 (period–energy curves) | `python reproduce/python_period_energy.py`, then `wolframscript -file code_wl/MakeFigures.wls` | `benchmarks/python_period_energy.csv`, `data/figures/fig5_period_energy.pdf` |
| — | precision and tolerance study (S3.2) | `wolframscript -file code_wl/MakeFigures.wls` | `data/figures/fig7_precision_and_tolerance.pdf` |
| Table 3 (tests) | Table S3 | `cd code_py && python -m pytest -q --cov=nonlinear_oscillator_lab --cov-branch`; `wolframscript -file code_wl/RunTests.wls` | Wolfram: `code_wl/verification_v0.2.1/round2/test_report.txt` |
| Table 4 (closed-form campaigns, 1436 problems) | Table S4 | `python reproduce/closed_form_campaign.py --workers 4` (about 15 min on 4 cores; resumable) | `data/results/closed_form_campaign_py_exactref.csv` and `_summary.json` |
| Table 4 caption (accuracy of the DOP853 reference) | S3.2 | `python reproduce/dop853_reference_check.py` | `data/results/dop853_reference_check.json` |
| — | S3.2 (25-digit references from decimal strings vs exact binary inputs) | `python reproduce/rereference_campaign.py` | `data/results/closed_form_campaign_py_exactref.csv` |
| — | Table S5 (nine benchmark periods) | `wolframscript -file code_wl/MakeFigures.wls`; Python column: `python reproduce/python_benchmarks.py` | `data/wolfram/benchmark_periods.csv`, `benchmarks/python_periods.csv` |
| Table 5 (period baselines) and trajectory timings (Section 3.3) | Table S6 | `python reproduce/baseline_comparison.py` | `data/results/v0.2.1_period/baseline_summary.json`, `baseline_periods.csv`, `baseline_trajectories.csv` |
| Section 3.2, Wolfram numbers (medians, campaign, timings) | S3.2 | `wolframscript -file code_wl/MakeFigures.wls` | `data/wolfram/numbers.tex` (`\pythonTests` there is a fixed placeholder; Python test counts come from pytest) |
| — | S7.2 (Algorithm 1 on the 1436 problems) | `python -I reproduce/independent/campaign_weierstrass.py --workers 4` | `data/results/campaign_weierstrass.csv` |
| — | S8 (independent checks of the formulas) | `python -I reproduce/independent/<script>.py` (see below) | printed output |

`code_wl/MakeFigures.wls` also writes a `figures/fig6_period_accuracy.pdf` of its own, with an
additional reference series; the figure of the paper and supplement is the one drawn by
`reproduce/period_accuracy_figure.py` from the same `period_accuracy.csv`.

`reproduce/independent/` contains checks written from the formulas of the supplement only (NumPy,
SciPy and mpmath; they do not import the package): `verify_sextic.py`, `verify_weierstrass.py`,
`algorithm_weierstrass.py`, `campaign_weierstrass.py --workers 4` and `verify_appendix_algorithm.py`.
`verify_sextic.py --tex <dir>` compares its polynomials with the LaTeX sources of the supplement,
which are not part of this archive.

The other files in `data/results/` (`closed_form_campaign_py.csv` and its summary, `baseline_*` at
the top level, `test_summary.json`, `python_period_energy_v0.2.0.csv`) are runs of version 0.2.0,
kept for reference; the paper and supplement do not quote them.

Timings depend on the machine. Errors against DOP853 and the Algorithm 1 statistics may differ in
the last digits between platforms; counts and representation decisions are reproducible.
