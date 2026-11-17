# Contributing to NonlinearOscillatorLab

Thank you for helping. Bug reports, test cases (especially orbits where a result looks wrong),
documentation fixes and new features within the [scope](docs/scope.md) are all welcome.

## Reporting a bug

Open an issue with: the model coefficients `a, b, c`, the initial state `x0, v0`, the exact call,
the full error message or the wrong value, the expected value if known, and your versions
(`python -c "import nonlinear_oscillator_lab, numpy, scipy, mpmath, sys; print(nonlinear_oscillator_lab.__version__, numpy.__version__, scipy.__version__, mpmath.__version__, sys.version)"`,
or `$Version` and `NOLVersion` in the Wolfram Language).

## Development set-up

```bash
python -m pip install -e "./code_py[test,bench]"
python -m pytest code_py/tests -q --cov=nonlinear_oscillator_lab --cov-branch
wolframscript -file code_wl/RunTests.wls       # if a Wolfram Engine is available
```

## Rules for pull requests

1. One topic per pull request, with a test that fails before and passes after the change.
   Mark tests with `unit`, `integration`, `regression` or `failure` (see `pyproject.toml`).
2. New numerical routines need an independent reference in the tests (closed form,
   arbitrary precision, or the other implementation), not only a self-consistency check.
3. A change of a published number (anything in `data/`) must be
   explained in `CHANGELOG.md`.
4. Keep the Python and Wolfram Language APIs in correspondence (`docs/api_*.md`); if only one
   side is changed, say so in the pull request.
5. Run the local test command above and review any changes to branch coverage.
6. By contributing you agree that your contribution is released under the MIT licence.

## Releases

Versions follow semantic versioning; every release is archived with a DOI and recorded in `CHANGELOG.md`.

## Conduct

Be respectful and constructive; harassment of any kind is not tolerated. Concerns can be sent
privately to the maintainer (ahsalass@unal.edu.co).
