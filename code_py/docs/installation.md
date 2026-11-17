# Installation

## Python (3.10 or newer)

Dependencies: NumPy ≥ 1.23, SciPy ≥ 1.10, mpmath ≥ 1.2 (installed automatically).

From the software bundle root:

```bash
python -m pip install ./code_py
python code_py/examples/quickstart.py
```

Run `python -m pip install './code_py[test]'` before `python -m pytest code_py/tests -q`.
The v0.2.2 Python wheel and suite were verified on macOS arm64 with Python 3.15.1.
Install `mkdocs-material` with pip, then run `mkdocs build --strict` from `code_py/`.

## Wolfram Language (Mathematica or Wolfram Engine 14.0 or newer)

Load the bundled Wolfram source directly:

```mathematica
Get["/path/to/NonlinearOscillatorLab-0.2.2/code_wl/NonlinearOscillatorLab.wl"]
NOLVersion
```

The Wolfram component is version 0.2.1. Install the included paclet and load its context,
using the absolute path to the extracted bundle:

```mathematica
PacletInstall["/path/to/NonlinearOscillatorLab-0.2.2/code_wl/paclets/NonlinearOscillatorLab-0.2.1.paclet"]
Needs["NonlinearOscillatorLab`"]
NOLPeriod[NOLModel[4, 0, 0], 2/5, 0, WorkingPrecision -> 40]
```

The last command gives Pi to 40 digits. The paclet was built, installed into a new paclet base
and loaded with `Needs` in a fresh Mathematica 15.0.0 kernel on Windows 64-bit.
The source suite separately passed 37/37 tests. Logs, reports and the environment record
are in `code_wl/verification/`.

To repeat construction and the installation check with a licensed kernel, run the following
from the bundle root:

```bash
test_root=$(cd "$(mktemp -d)" && pwd -P)
mkdir "$test_root/paclets"
NOL_PACLET_OUTPUT="$test_root/archives" wolframscript -file code_wl/BuildPaclet.wls
NOL_PACLET_PATH="$test_root/archives/NonlinearOscillatorLab-0.2.1.paclet" \
NOL_PACLET_BASE="$test_root/paclets" \
  WolframKernel -noinit -pacletbase "$test_root/paclets" -script code_wl/VerifyPaclet.wls
```

The build stages `PacletInfo.wl`, the MIT licence and `Kernel/NonlinearOscillatorLab.wl`.
`CreatePacletArchive` returns a string path. The verification script passes that archive path
to `PacletInstall`, then loads the installed context with `Needs` in the fresh kernel and
checks the harmonic period. [Wolfram build
reference](https://reference.wolfram.com/language/ref/CreatePacletArchive.html)

The `Get` run and 37 passing tests are in `code_wl/verification/test_report.txt`.
Use `wolframscript -file code_wl/RunTests.wls` to rerun them with a licensed engine.

### Licensing of the Wolfram Engine

The Wolfram Language package is MIT-licensed, but it runs on proprietary software:

* **Mathematica / Wolfram|One**, under the user's own licence (many universities hold site
  licences that also cover the Wolfram Engine);
* the **free Wolfram Engine for Developers / Community Edition**, whose terms allow
  pre-production development, personal non-commercial projects and exploring the
  language, but *not* running it to produce output for commercial or organisational use or
  in production. Check the current terms at <https://www.wolfram.com/engine/faq/>
  before using it for published work.

The Python package needs no proprietary software. Since v0.2.0 it contains every capability of
the Wolfram package except symbolic (exact) equilibria and the phase-portrait plot; arbitrary
precision is provided through mpmath (`period_mp`) and closed forms through `exact_solution`.
