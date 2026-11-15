# Changelog

All notable changes are listed here. Versions follow [semantic versioning](https://semver.org).

## [0.2.1] - 2026-11 (second revision of the SoftwareX article)

### Changed
- Closed-form failure message (Python `ClosedFormNotAvailable`, Wolfram `Failure["NotCovered"]`):
  "No candidate representation passed the acceptance checks". The v0.2.0 text claimed that no
  real representation exists, which the finite acceptance checks cannot establish. A parenthesis
  "(unbounded motion in all tested cases)", used in an intermediate 0.2.1 build, was dropped: the
  failure also occurs at separatrix levels (e.g. (1, -4, 3, 0, sqrt(4/27))) and for bounded
  orbits outside the sampled ranges (e.g. (1, 0, 1, 100, 1)); use `motion_type` for the motion.
- Archive layout: `code_py/`, `code_wl/`, `data/`, `reproduce/`; scripts write to
  `reproduce/outputs/` and never overwrite the archived data.
- The 15 legacy v0.1.0 tests (`tests/test_core.py`) are no longer distributed; the suite is the
  158-test suite of v0.2.0 plus 8 tests for the corrections below (166 in total).

### Fixed
- Python `period`: E - V(x) = (x - x-)(x+ - x) Q(x) is deflated backwards (from the constant term)
  when the forward quotient does not reproduce E - V to 1e-12 relative to its terms. For
  (a, b, c, x0, v0) = (1, -1, 1e-8, 0, 1) v0.2.0 raised "Degenerate turning point" (T = 11.41683767063561);
  for c = 1e-7 it was 4e-3 off. On the 320 baseline orbits only 12 outer double-well periods next to
  the separatrix change; their maximum error falls from 1.1e-11 to 5e-14
  (`data/results/v0.2.1_period/`). `tests/test_v021_corrections.py`.
- Wolfram `NOLPeriod`: the same backward deflation (`deflatedQuotient`). At machine precision the
  period of (1, -1, 1e-8, 0, 1) was 8.50929202969664, 25% off, without any message; it is now
  11.41683767063562 (relative error 9e-16). On the 320 baseline orbits the maximum Wolfram error
  falls from 2.4e-11 to 1.3e-11 and the outer double-well errors to 5.0e-14.
- Wolfram: the degenerate-equilibria correction (`NOLClassifyEnergy` key `"Degenerate"`, critical
  levels in `"OnSeparatrixEnergy"`) is merged into the package; three tests added. Three more
  cover the failure message and the period above (37 in total). Verified with Mathematica 15.0.0
  (code_wl/verification_v0.2.1/: first round; round2/: final files, 37/37).
- Campaign script: the 25-digit mpmath reference starts from the exact binary inputs
  (`mpf(float(v))`). The decimal strings used before moved the reference by up to 8.1e-7 near
  barrier energies; the closed-form solutions themselves are unchanged.

### Reproduction
- `reproduce/period_accuracy_figure.py` draws the period-accuracy figure (Figure 3 of the paper,
  Figure S5 of the supplement): panel (a) from `data/wolfram/period_accuracy.csv`, panel (b) computed
  with `period_mp` at 40 digits (`data/results/period_accuracy_panel_b.csv`, equal to the Wolfram
  values in double precision and to the elliptic-K periods to 3e-40). Its output replaces
  `data/figures/fig6_period_accuracy.{pdf,png}`, which `code_wl/MakeFigures.wls` drew with an
  additional v0.1.0 series.
- `reproduce/README.md` maps every table and figure of the paper and supplement to its command and
  archived file; Table 5 comes from `data/results/v0.2.1_period/`.
- `reproduce/make_revision_tables.py` removed: it wrote LaTeX tables from the v0.2.0 runs.
- Documentation: 166 tests in `code_py/docs/installation.md`; paths `code_wl/` in script headers.

## [0.2.0] - 2026-11 (first revision of the SoftwareX article)

### Added
- `closed_form.py`: Python port of the closed-form solver (`exact_solution`,
  `duffing_solution`), previously available only in the Wolfram Language. It selects the same
  representation as the Wolfram solver in all 200 cases of the published campaign.
- `period_mp`: arbitrary-precision period (mpmath), the analogue of
  `NOLPeriod[..., WorkingPrecision -> d]`; reproduces the published 40-digit benchmark periods.
- `motion_type`: equilibrium / separatrix / bounded / unbounded classification.
- `forced.py`: damped, harmonically forced extension (`simulate_forced`, `poincare_section`,
  `energy_balance_residual`) in autonomous form with the phase variable θ = ωt + φ.
- Test suite reorganised by type (unit, integration, regression, failure): 158 tests, 92 %
  branch coverage (was 15 tests).
- Campaigns (`benchmarks/closed_form_campaign.py`) and baseline benchmarks against SciPy
  (`benchmarks/baseline_comparison.py`), with results in `benchmarks/results/`.
- Documentation in `docs/` (installation, tutorial, API references, troubleshooting,
  reproducibility, scope/roadmap), `CONTRIBUTING.md`, `RELEASING.md`, CI workflows,
  archive metadata, `mathematica/PacletInfo.wl`, `examples/quickstart.py`.

### Fixed
- Degenerate equilibria (V'' = 0, e.g. b² = 4ac) were not recognised as critical levels:
  `period(Model(1, -2, 1), 0, 3**-0.5)` returned a finite, meaningless value (5684.85) instead
  of reporting an infinite period. `classify_energy` now lists them under `degenerate` and
  includes their energies in `on_separatrix_energy`; `period` also rejects multiple turning
  points. A corresponding patch for the Wolfram package is in
  `mathematica/patches/degenerate_equilibria.wl`.
- Oscillations about a minimum whose energy lies less than `tol` above the minimum (amplitude
  about 1e-8) were classified as unbounded, and `period` refused them, because the test
  E - V > 0 inside the interval used the absolute tolerance `tol`. The test now uses a
  rounding-scaled threshold; oscillations too small to resolve in double precision return
  the harmonic limit 2 pi / sqrt(V''(x_min)) (relative error O(tol)).

### Improved (Python closed form)
- The construction stage (roots of the cubic, parameters A, lambda, mu, p, q, modulus,
  frequency, phase) is computed in 40-digit arithmetic, the denominator
  1 + lambda y^2 + mu y^4 is evaluated in factored form, and Jacobi functions whose
  (transformed) parameter lies within 1e-8 of 1 are evaluated with mpmath. Near separatrices
  these steps previously limited the accuracy (e_x up to 1e-1 within 1e-7 of a barrier
  energy); the maximum is now 4.5e-8 in the near-separatrix campaign and 8.1e-7 over all
  campaigns (within 1e-10 of an escape barrier). The selected representation is unchanged
  (200/200 agreement with the Wolfram solver).
- Exact equilibria return the constant solution instead of raising.

### Changed
- `classify_energy` returns the additional keys `minima` and `degenerate`.
- `mpmath` is now a dependency.

## [0.1.0] - 2026-10-19
- Initial release accompanying the SoftwareX submission.
