# Data: inputs and archived outputs

Nothing in this folder is overwritten by the Python scripts, which write to reproduce/outputs/; reproduce/README.md maps every table and figure to its file here. The Wolfram
script code_wl/MakeFigures.wls writes to figures/ and tables/ at the archive root. The Python outputs that
the Wolfram scripts read are in ../benchmarks/ and ../examples/.

| Folder | Contents | Produced by |
|---|---|---|
| `wolfram/` | Wolfram campaign of 200 problems (closed_form_campaign.csv, SeedRandom[20261019]), benchmark periods, period-accuracy data, environment, LaTeX numbers | code_wl/MakeFigures.wls (Mathematica 15.0, 15 Nov 2026) |
| `results/` | **Quoted in the paper:** closed_form_campaign_py_exactref.csv and summary (Table 4 / S4; 25-digit references from exact inputs), dop853_reference_check.json (Table 4 caption), campaign_weierstrass.csv (Algorithm 1, S7.2), v0.2.1_period/ (baseline run of version 0.2.1: Table 5 / S6 and the trajectory timings), period_accuracy_panel_b.csv (panel (b) of Fig. 3 / S5, 40-digit `period_mp`). **Kept for reference (version 0.2.0 runs, not quoted):** closed_form_campaign_py.csv and summary, baseline_periods.csv, baseline_trajectories.csv, baseline_summary.json, test_summary.json, python_period_energy_v0.2.0.csv | reproduce/ (see its README) |
| `figures/` | Figures of the paper and supplement (PDF and PNG) | code_wl/MakeFigures.wls; fig6_period_accuracy (Fig. 3 / S5) by reproduce/period_accuracy_figure.py |
