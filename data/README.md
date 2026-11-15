# Data: inputs and archived outputs

Nothing in this folder is overwritten by the Python scripts, which write to reproduce/outputs/. The Wolfram
script code_wl/MakeFigures.wls writes to figures/ and tables/ at the archive root. The Python outputs that
the Wolfram scripts read are in ../benchmarks/ and ../examples/.

| Folder | Contents | Produced by |
|---|---|---|
| `wolfram/` | Wolfram campaign of 200 problems (closed_form_campaign.csv, SeedRandom[20261019]), benchmark periods, period-accuracy data, environment, LaTeX numbers | code_wl/MakeFigures.wls (Mathematica 15.0, 15 Nov 2026) |
| `results/` | closed_form_campaign_py.csv and summary (v0.2.0 run), closed_form_campaign_py_exactref.csv and summary (25-digit references from exact inputs, as in Table 4 of the paper), dop853_reference_check.json, campaign_weierstrass.csv (Algorithm 1), baseline_periods.csv, baseline_trajectories.csv, baseline_summary.json, test_summary.json (v0.2.0 runs); v0.2.1_period/ (baseline rerun with the v0.2.1 `period`; errors as in Table 5, timings machine-dependent); python_period_energy_v0.2.0.csv (replaced in ../benchmarks/ by the v0.2.1 values) | reproduce/ (see its README) |
| `figures/` | Figures of the paper and supplement (PDF and PNG) | code_wl/MakeFigures.wls |
