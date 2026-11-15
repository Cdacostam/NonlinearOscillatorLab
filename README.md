# NonlinearOscillatorLab 0.2.1

Python and Wolfram Language toolkit for the conservative Duffing and cubic–quintic oscillator
x'' + a x + b x^3 + c x^5 = 0: energy-based motion classification, singularity-free and
arbitrary-precision periods, closed-form solutions x = A y / sqrt(1 + lambda y^2 + mu y^4), and a
numerical extension to damped, harmonically forced motion. MIT License.

Archive of the version described in the SoftwareX article
"NonlinearOscillatorLab: A Python and Wolfram Language toolkit for energy-based classification, accurate period computation and closed-form solution of Duffing and cubic–quintic oscillators".
Companion demonstrator (no installation): https://fizmako.com/softwarex/. It is an independent
server-side implementation for single initial-value problems, not the code archived here; the
reproduction scripts in `reproduce/` are the reproducible capsule.

## Contents

| Folder | Contents |
|---|---|
| `code_py/` | Python package `nonlinear_oscillator_lab` (src/, 166 tests, examples, documentation) |
| `code_wl/` | Wolfram Language package, paclet metadata, 37 tests, figure scripts, v0.2.1 verification evidence |
| `data/` | Inputs and archived outputs: Wolfram campaign and tables, campaign and baseline results, figures |
| `benchmarks/`, `examples/` | Python outputs read by the Wolfram scripts (python_periods.csv, period–energy curves, orbit.csv), kept at the paths the tested scripts expect |
| `reproduce/` | One script per table or figure of the paper; new outputs go to `reproduce/outputs/` |

`MANIFEST.sha256` lists the SHA-256 of every file (`sha256sum -c MANIFEST.sha256`).

## Quick start

```bash
python -m pip install ./code_py            # Python >= 3.10, NumPy, SciPy, mpmath
cd code_py && python -m pytest -q          # 166 tests
python examples/quickstart.py
```

```
wolframscript -file code_wl/RunTests.wls   # Mathematica or Wolfram Engine 14+; 37 tests
```

The Wolfram suite of v0.2.1 (37/37, Mathematica 15.0.0, Windows 64-bit) was run on these files by the
Wolfram co-worker; the evidence is in `code_wl/verification_v0.2.1/round2/`.

`reproduce/README.md` maps every table and figure of the paper to its command.

## Version 0.2.1

Python and Wolfram: the closed-form failure message states only that no candidate passed the
acceptance checks (it no longer claims that no real representation exists, nor anything about the
motion). Python: `period` uses backward deflation where forward deflation of E - V is inaccurate
(wide wells such as (a, b, c) = (1, -1, 1e-8), and the outer double-well orbits next to the
separatrix). Wolfram: `NOLPeriod` uses the same backward deflation, and critical
levels of degenerate equilibria are detected (formerly a separate patch). Campaign script: the
25-digit reference starts from the exact binary inputs. See `CHANGELOG.md`.

## Citation

See `CITATION.cff`. Author of the software: Alvaro H. Salas (Universidad Nacional de Colombia, Sede
Manizales).
