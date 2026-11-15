# Wolfram Language API reference (`NonlinearOscillatorLab.wl`)

In v0.2.1 the Wolfram Language package treats degenerate equilibria as critical levels and words the
closed-form failure message as in Python; in v0.2.0 it was unchanged except for the paclet metadata
(`PacletInfo.wl`); its full reference with examples is Appendix A of the paper. Every symbol
has a `::usage` message (`?NOLPeriod`).

| Call | Result |
|---|---|
| `NOLModel[a, b, c]` (`NOLModel[a, b]` sets $c=0$) | inert model expression |
| `NOLPotential[m, x]`, `NOLForce[m, x]`, `NOLEnergy[m, x, v]`, `NOLCurvature[m, x]` | $V$, $-V'$, $v^2/2+V$, $V''$ |
| `NOLStationaryPoints[m]` | sorted real roots of $V'$ (exact for exact input) |
| `NOLClassifyEnergy[m, x0, v0]` | `<|"Energy","Equilibria","Minima","Maxima","OnSeparatrixEnergy"|>` |
| `NOLTurningPoints[m, E]` | sorted real roots of $V=E$ |
| `NOLBoundedComponent[m, x0, v0]` | `{x-, x+}` or `Failure["Unbounded", ...]` |
| `NOLPeriod[m, x0, v0, opts]` | period; options `WorkingPrecision`, `PrecisionGoal`, `"Tolerance"`; `Infinity` at an equilibrium; `Failure` on separatrix / unbounded motion |
| `NOLPeriodFromEnergy[m, E, x0, opts]` | period of the energy-$E$ orbit through the well containing $x_0$ |
| `NOLDuffingExactPeriod[a, b, A]`, `NOLDuffingExactSolution[a, b, A, t]` | closed-form $c=0$ references |
| `NOLSeparatrixPeriodAsymptotics[m, xs, dE, k]` | $(k/\lambda)\ln(1/\lvert dE\rvert)$ |
| `NOLDuffingSolution[p, q, y0, y1]` | closed-form Duffing solution (Association) |
| `NOLExactSolution[m, x0, v0]` | closed form, key `"x"` is a pure function; or `Failure["NotCovered", ...]` |
| `NOLExactSolutionError[m, x0, v0, sol, tmax]` | max deviation from a 32-digit `NDSolve` reference |
| `NOLSimulate[m, x0, v0, tEnd, opts]` | `<|"x","v","tEnd",...|>` (explicit RK8, dense output) |
| `NOLEnergyDrift[m, sol]` | $\delta_E$ on 2001 samples |
| `NOLPhasePortrait[m, {x1,x2}, {v1,v2}]` | level sets, separatrices in red |

Correspondence with Python: `NOLPeriod[..., WorkingPrecision -> d]` ↔ `period_mp(..., dps=d)`;
`NOLExactSolution` ↔ `exact_solution` (same representation in all 200 published campaign
cases); `Failure[...]` ↔ exceptions.

Known differences between the Wolfram code and the Python package (v0.2.1):

* The Python closed-form solver computes its parameters in 40-digit arithmetic and
  evaluates the denominator in factored form, which keeps it accurate close to separatrices
  (see `CHANGELOG.md`); the Wolfram solver works in machine precision.
* Very small oscillations about a minimum (energy within `"Tolerance"` of the minimum) are
  classified as bounded in Python (harmonic-limit period).
* Since v0.2.1 both packages treat the levels of degenerate equilibria ($V''=0$, e.g.
$b^2=4ac$) as critical: `NOLClassifyEnergy` lists them under `"Degenerate"` and `NOLPeriod` fails
there. The former separate patch, `code_wl/patches/degenerate_equilibria.wl`, is merged and kept for
the record (see `CHANGELOG.md`).
