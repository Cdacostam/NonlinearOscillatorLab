# NonlinearOscillatorLab

NonlinearOscillatorLab is an MIT-licensed toolkit for the conservative oscillator

$$\ddot x + a x + b x^3 + c x^5 = 0,\qquad V(x)=\tfrac a2x^2+\tfrac b4x^4+\tfrac c6x^6,$$

with two implementations whose common numerical results are compared: a **Python package** (`nonlinear_oscillator_lab`)
and a **Wolfram Language package** (`code_wl/NonlinearOscillatorLab.wl`).

| Page | Content |
|---|---|
| [Installation](installation.md) | pip, editable install, Wolfram Engine / Mathematica, licences |
| [Tutorial](tutorial.md) | a complete workflow in both languages |
| [Python API](api_python.md) | every public function, arguments, returns, exceptions |
| [Wolfram Language API](api_wolfram.md) | every public symbol, options, `Failure` objects |
| [Troubleshooting and failure handling](troubleshooting.md) | what each error means and what to do |
| [Validation and reproducibility](reproducibility.md) | tests, campaigns, benchmarks, how to rerun the archived benchmarks |
| [Scope and roadmap](scope.md) | target users, what the package does not do, planned work |

## Who is it for?

* **Researchers developing analytic approximations** (harmonic balance, perturbation series,
  Jacobi-elliptic ansätze) for Duffing-type equations, who need trustworthy reference periods
  and trajectories, including near separatrices and in multi-well potentials.
* **Numerical analysts** testing integrators or quadrature rules on a problem family with
  closed-form and arbitrary-precision references.
* **Teachers** of nonlinear dynamics (phase planes, separatrices, period-energy curves,
  multi-well potentials, forced response).

## What does it add to SciPy / Mathematica?

General-purpose tools integrate any ODE but leave the *energy geometry* to the user: which
interval of the energy surface contains the initial state, whether the orbit is closed, whether
the energy is critical, and how to integrate the period without the endpoint singularity.
NonlinearOscillatorLab packages these decisions, returns explicit failures instead of silently
wrong numbers, and supplies closed forms and arbitrary-precision references for verification.
See Table 1 of the paper and [reproducibility](reproducibility.md) for quantitative comparisons.
