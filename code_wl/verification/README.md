# Wolfram package and paclet verification

The files record the Mathematica 15.0.0 run on Windows 64-bit on 16 November 2026.
The paclet is supplied in `../paclets/NonlinearOscillatorLab-0.2.1.paclet`.

| File | Recorded result |
|---|---|
| `01-build.log` | Path of the archive returned by the build |
| `02-paclet.log` | PASS after installation in a new paclet base, loading with Needs and a harmonic-period check |
| `03-tests.log` | Source-package suite, 37 successful tests |
| `test_report.txt`, `test_report.csv` | Individual outcomes and times for the 37 source-package tests |
| `ENTORNO.txt` | Kernel and platform, input ZIP hash, invocation options and exit codes |

The installation check used a fresh kernel with `-noinit` and a new paclet base
specified by `-pacletbase`. It called `PacletInstall` on the local archive, loaded the
context with `Needs` and computed the harmonic period at 40 digits. The tolerance
against Pi was 1e-30. The 37-test suite used source loading with `Get` in a separate process.

The paclet SHA-256 is
`e5d185a2de5203b1d6684d9cdad1f3a6b759a406d59fd5f2c7941ad3c34585c6`.
Its kernel source, PacletInfo.wl and licence match the corresponding source files in the
software bundle byte for byte. Build and verification commands are in
`../../code_py/docs/installation.md`.
