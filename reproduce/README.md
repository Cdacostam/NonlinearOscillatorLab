# Reproduce the tables and figures

Run from the archive root after `python -m pip install ./code_py` (Python) or with `wolframscript`
(Wolfram). Outputs go to `reproduce/outputs/` (or to the folder in the environment variable
`NOL_OUTPUTS`) and can be compared with the archived files in `data/`.

| Paper | Command | Archived output |
|---|---|---|
| Figures 1–4, Supplementary Figures and Table S6 | `wolframscript -file code_wl/RunAll.wls` (writes `figures/`, `tables/`) | `data/figures/`, `data/wolfram/` |
| Table 3 (tests) | `cd code_py && python -m pytest -q --cov=nonlinear_oscillator_lab --cov-branch`; `wolframscript -file code_wl/RunTests.wls` | `data/results/test_summary.json`, `code_wl/test_report.*` |
| Table 4 (closed-form campaigns, 1436 problems) | `python reproduce/closed_form_campaign.py --workers 4` (about 15 min on 4 cores; resumable) | `data/results/closed_form_campaign_py_exactref.csv` |
| Table 4, reference check of v0.2.0 rows | `python reproduce/rereference_campaign.py` | `data/results/closed_form_campaign_py_exactref.csv` |
| Table 4 caption, DOP853 reference accuracy | `python reproduce/dop853_reference_check.py` | `data/results/dop853_reference_check.json` |
| Table 5 (period baselines) and trajectory timing | `python reproduce/baseline_comparison.py` | `data/results/baseline_*` |
| LaTeX tables from the result files | `python reproduce/make_revision_tables.py` | — |
| Python periods read by the Wolfram tests | `python reproduce/python_benchmarks.py` | `benchmarks/python_periods.csv` |
| Period–energy curves (Supplementary figure) | `python reproduce/python_period_energy.py` | `benchmarks/python_period_energy.csv` |
| Algorithm 1 and appendix formulas (Supplement S8–S9) | `python -I reproduce/independent/<script>.py` (see below) | `data/results/campaign_weierstrass.csv` |

`reproduce/independent/` contains checks written from the formulas of the supplement only (NumPy
and SciPy; they do not import the package): `verify_sextic.py`, `verify_weierstrass.py`,
`algorithm_weierstrass.py`, `campaign_weierstrass.py --workers 4` and
`verify_appendix_algorithm.py`. `verify_sextic.py --tex <dir>` compares its polynomials with the
LaTeX sources of the supplement, which are not part of this archive.

Timings depend on the machine. Errors against DOP853 and the Algorithm 1 statistics may differ in
the last digits between platforms; counts and representation decisions are reproducible.
