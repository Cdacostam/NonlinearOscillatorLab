# Installation

## Python (3.10 or newer)

Dependencies: NumPy ≥ 1.23, SciPy ≥ 1.10, mpmath ≥ 1.2 (installed automatically).

```bash
# from the archived release / a clone of the repository
python -m pip install .            # or: python -m pip install -e ".[test]"  (development)
python -m pytest -q                # 158 tests, under a minute
python examples/quickstart.py
```

Once the package is published on PyPI (see `RELEASING.md`):

```bash
python -m pip install nonlinear-oscillator-lab
```

Verified environments: Linux x86-64 with Python 3.13, NumPy 2.5, SciPy 1.18 (this release);
the continuous-integration workflow of the development repository tests Python 3.10–3.13 on
Linux, macOS and Windows.

## Wolfram Language (Mathematica or Wolfram Engine 14.0 or newer)

No installation is required:

```mathematica
Get["/path/to/NonlinearOscillatorLab-0.2.1/code_wl/NonlinearOscillatorLab.wl"]
NOLVersion
```

or install it as a paclet (the folder `code_wl/` contains `PacletInfo.wl`):

```mathematica
PacletInstall["/path/to/NonlinearOscillatorLab-0.2.1.paclet"]   (* built with CreatePacletArchive *)
Needs["NonlinearOscillatorLab`"]
```

Verification: `wolframscript -file code_wl/RunTests.wls` (37 tests).

### Licensing of the Wolfram Engine

The Wolfram Language package is MIT-licensed, but it runs on proprietary software:

* **Mathematica / Wolfram|One**, under the user's own licence (many universities hold site
  licences that also cover the Wolfram Engine);
* the **free Wolfram Engine for Developers / Community Edition**, whose terms allow
  pre-production development, personal non-commercial projects and exploring the
  language, but *not* running it to produce output for commercial or organisational use or
  in production. Check the current terms at <https://www.wolfram.com/legal/terms/wolfram-engine.html>
  before using it for published work.

The Python package needs no proprietary software. Since v0.2.0 it contains every capability of
the Wolfram package except symbolic (exact) equilibria and the phase-portrait plot; arbitrary
precision is provided through mpmath (`period_mp`) and closed forms through `exact_solution`.
