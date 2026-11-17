# NonlinearOscillatorLab 0.2.2

Python and Wolfram Language toolkit for the conservative Duffing and cubic–quintic oscillator
x'' + a x + b x^3 + c x^5 = 0: energy-based motion classification, singularity-free and
arbitrary-precision periods, closed-form solutions x = A y / sqrt(1 + lambda y^2 + mu y^4), and a
numerical extension to damped, harmonically forced motion. MIT License.

Archive of the version described in the SoftwareX article
"NonlinearOscillatorLab: A Python and Wolfram Language toolkit for energy-based classification, accurate period computation and closed-form solution of Duffing and cubic–quintic oscillators".
The companion demonstrator at https://fizmako.com/softwarex/ illustrates single initial-value
problems with its own server-side implementation. The commands and fixed inputs for
reproduction are supplied in `reproduce/` and `data/`. The accompanying
`supplementary_material.pdf` provides the formulas, algorithms and verification methods.

## Contents

| Folder | Contents |
|---|---|
| `code_py/` | Python package `nonlinear_oscillator_lab` (src/, 171 tests, examples, documentation) |
| `code_wl/` | Wolfram Language source and built paclet, 37 tests, figure scripts and verification logs |
| `data/` | Inputs and archived outputs: Wolfram campaign and tables, campaign and baseline results, figures |
| `benchmarks/`, `examples/` | Python outputs read by the Wolfram scripts (python_periods.csv, period–energy curves, orbit.csv), kept at the paths the tested scripts expect |
| `reproduce/` | Scripts for the tables and figures of the paper; new outputs go to `reproduce/outputs/`. `reproduce/README.md` maps every table and figure to its command and archived file |

`MANIFEST.sha256` lists the SHA-256 of every file (`sha256sum -c MANIFEST.sha256`).

## Quick start

```bash
python -m pip install ./code_py            # Python >= 3.10, NumPy, SciPy, mpmath
python code_py/examples/quickstart.py
python -m pip install './code_py[test]'
python -m pytest code_py/tests -q -o pythonpath=
```

```
wolframscript -file code_wl/RunTests.wls   # Mathematica or Wolfram Engine 14+; 37 tests
```

The Wolfram component and its paclet are version 0.2.1. The paclet is in
`code_wl/paclets/NonlinearOscillatorLab-0.2.1.paclet`. It was built, installed into a new
paclet base and loaded with `Needs` in a fresh Mathematica 15.0.0 kernel on Windows 64-bit.
The harmonic period at 40 digits agreed with Pi within 1e-30. The source suite separately
passed 37/37 tests. The logs, reports and environment record are in `code_wl/verification/`.
See `code_py/docs/installation.md` for installation and verification commands.

`reproduce/README.md` maps every table and figure of the paper and the supplement to its command and
archived file. The period-accuracy figure (`python reproduce/period_accuracy_figure.py`) also needs
`matplotlib`.

## Availability

The software archive for the article is identified by
[DOI 10.5281/zenodo.24248864](https://doi.org/10.5281/zenodo.24248864).
The [development repository](https://github.com/Cdacostam/NonlinearOscillatorLab)
provides the issue tracker and contribution workflow.

## Citation

See `CITATION.cff`. Authors of the software: Alvaro H. Salas (Universidad Nacional de Colombia, Sede
Manizales; ORCID 0000-0001-9343-6062) and Carlos Daniel Acosta-Medina (Universidad Nacional de
Colombia; ORCID 0000-0002-6477-8984).
