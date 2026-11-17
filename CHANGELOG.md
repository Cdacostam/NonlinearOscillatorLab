# Changelog

## 0.2.2

- Support scalar and array Fraction coordinates in double-precision potential and period calculations.
- Detect unresolved minima in the factored binary64 period integrand and require convergence from
  QUADPACK. The error message directs users to the arbitrary-precision period routine.
- Recover representable potential values when intermediate powers overflow and reject nonfinite
  stationary or initial energies explicitly.
- Resolve repeated stationary roots using the quadratic discriminant and evaluate their curvature
  analytically in the energy classifier and period routines.
- Evaluate the analytic endpoint limits in QAWS and provide a divided-difference QAWS comparison.
  Compute period discrepancies from the native numerical outputs.
- Supply the fixed binary64 inputs for the 1436-case trajectory campaign, the reference comparison,
  the period and trajectory benchmarks, and the scripts that generate the manuscript tables.
- Document error metrics, failure handling, installation and reproduction commands.
- Supply the Wolfram component and its paclet as version 0.2.1 within software bundle 0.2.2.
  Include scripts for paclet construction and installation checks in a clean kernel.
- Distribute the built paclet and the recorded Mathematica 15.0.0 Windows 64-bit verification.
  Installation in a new paclet base, loading with Needs and the 40-digit harmonic-period check passed.
  The source-package suite separately passed 37 tests.
- Include the installed-wheel Python test and coverage reports for 171 passing tests, with
  92.43% combined coverage and 175 of 204 branches covered.
