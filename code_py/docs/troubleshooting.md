# Troubleshooting and failure handling

| Symptom | Meaning | What to do |
|---|---|---|
| `ValueError: Separatrix energy: infinite period` / `NOLPeriod::sep` | $E$ coincides (within `tol`, default $10^{-9}$ relative) with a barrier energy or with the energy of a degenerate equilibrium. The period is infinite. | If you really mean an energy just off the separatrix, lower `tol` (`period(..., tol=1e-14)`, `period_mp(..., tol=1e-30)`, `"Tolerance" -> 10^-30`) and use arbitrary precision; the error near a separatrix grows like $\epsilon_{\rm mach}/\lvert E-E_s\rvert$ in double precision. |
| `ValueError: Initial state is not in a verified bounded component` / `NOLPeriod::unb` | No unique bounded admissible interval containing the initial state could be verified. | Check `motion_type(m, x0, v0)`. Use `simulate` or try `exact_solution` on its admissible interval instead of `period`. |
| `ValueError: Degenerate (multiple) turning point` | $V'(x_\pm)=0$: the orbit ends at an equilibrium, period infinite. | As for separatrices. |
| `period` returns `inf` | $(x_0,v_0)$ is an equilibrium. | Expected. |
| `ClosedFormNotAvailable` / `Failure["NotCovered", ...]` | No candidate representation passed the acceptance checks. Use the energy classifier to determine the motion type. | Use `simulate` / `NOLSimulate`. The coverage observed in the campaigns is listed in the paper (Table 4). |
| `ValueError: Period unresolved in double precision` | The factored integrand is unresolved at the computed precision or QUADPACK did not meet its tolerance. | Use `period_mp` with the same inputs and compare increasing precision. |
| `ValueError: Initial energy outside ...` or `Stationary energy outside ...` | Energy outside the finite double range. | Rescale the model or use arbitrary precision for a period. |
| `RuntimeError` from `simulate` | The ODE integrator could not reach `t_end`, usually finite-time blow-up. | Shorten `t_end`; check `motion_type`. |
| Period near a separatrix differs from Wolfram in the last digits | The period becomes sensitive to small energy and turning-point errors. | Use `period_mp` / `WorkingPrecision`. |
| `period_mp` with `0.4` gives a different 30th digit than with `"0.4"` | `0.4` is the binary float $0.40000000000000002220\ldots$ | Use strings or `Fraction` for exact decimals. |
| Plotting the Wolfram closed form shows spurious gaps | `Plot` treats `Abs`/`Sign` as discontinuities. | `Plot[..., Exclusions -> None]`. |
| `wolframscript: command not found` | No Wolfram Engine. | Install Mathematica or the Wolfram Engine (see [installation](installation.md), including licence terms), or use the Python package. |

Reporting problems: write to the author (see `CITATION.cff`) with the model coefficients,
the initial state, the call and the full error message (see `CONTRIBUTING.md`).
