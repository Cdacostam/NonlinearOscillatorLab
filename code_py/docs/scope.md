# Scope, target users and roadmap

## In scope (v0.2.2)

* Autonomous conservative oscillators $\ddot x+ax+bx^3+cx^5=0$ with arbitrary real $a,b,c$:
  energy geometry, motion classification, periods (double and arbitrary precision),
  closed forms, integration and energy diagnostics.
* Numerical extension to linear damping and harmonic forcing,
  $\ddot x+\delta\dot x+ax+bx^3+cx^5=F\cos(\omega t+\varphi)$: integration in autonomous form,
  stroboscopic Poincaré map and energy-balance check.

## Out of scope

* Coupled multi-degree-of-freedom systems, unilateral contact/impacts, non-polynomial
  restoring forces, parametric excitation.
* Continuation and bifurcation analysis (use AUTO-07p, MatCont, BifurcationKit.jl, PyDSTool/PyCont).
* Interval-certified error bounds.
* Energy classification, periods and closed forms for the damped/forced system (no
  conserved energy there).

## Roadmap

| Version | Planned |
|---|---|
| 0.3 | General odd and even polynomial potentials $V(x)=\sum_k c_k x^k$ for the period and classification routines; event handling (impacts) in `simulate`; Wolfram counterpart of the forced extension. |
| 0.4 | Frequency-response curves of the forced system by numerical continuation of the Poincaré-map fixed points; interval enclosures of periods (python-flint / Arb). |

## Maintenance

The package is maintained by the author (Universidad Nacional de Colombia). Policy:
semantic versioning; every release is archived with a DOI and listed in
`CHANGELOG.md`; bug reports are answered by the author; contributors run the local Python tests before submitting changes; supported versions follow the NumPy/SciPy support
window. See `CONTRIBUTING.md`.
