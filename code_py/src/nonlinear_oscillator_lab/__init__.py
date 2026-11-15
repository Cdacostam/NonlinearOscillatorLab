"""NonlinearOscillatorLab: conservative polynomial-oscillator diagnostics.

Python API (v0.2.0)
-------------------
Energy geometry : Model, stationary_points, classify_energy, turning_points, motion_type
Periods         : period (double precision), period_mp (arbitrary precision, mpmath)
Closed forms    : exact_solution, duffing_solution (Python port of the Wolfram solver)
Simulation      : simulate, energy_drift
Forced/damped   : simulate_forced, poincare_section, energy_balance_residual
"""
from .core import (Model, simulate, stationary_points, classify_energy, turning_points,
                   period, period_mp, motion_type, energy_drift)
from .closed_form import (exact_solution, duffing_solution, ClosedFormSolution,
                          DuffingSolution, ClosedFormNotAvailable)
from .forced import simulate_forced, poincare_section, energy_balance_residual

__version__ = "0.2.1"
__all__ = ["Model", "simulate", "stationary_points", "classify_energy", "turning_points",
           "period", "period_mp", "motion_type", "energy_drift",
           "exact_solution", "duffing_solution", "ClosedFormSolution", "DuffingSolution",
           "ClosedFormNotAvailable",
           "simulate_forced", "poincare_section", "energy_balance_residual"]
