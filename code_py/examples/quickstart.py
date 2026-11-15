"""Quick start: a complete NonlinearOscillatorLab workflow in Python (v0.2.0).

Triple-well cubic-quintic oscillator x'' + x - 4 x^3 + 3 x^5 = 0.
Run:  python examples/quickstart.py
"""
from fractions import Fraction

import numpy as np

from nonlinear_oscillator_lab import (Model, classify_energy, motion_type, period, period_mp,
                                      exact_solution, simulate, energy_drift, ClosedFormNotAvailable,
                                      simulate_forced, poincare_section)

m = Model(a=1, b=-4, c=3)

# 1. Energy geometry: equilibria, barriers, separatrix energy
info = classify_energy(m, x0=1.1, v0=0.0)
print("equilibria      :", np.round(info["equilibria"], 6))
print("barriers (x, E) :", [(round(x, 6), round(e, 6)) for x, e in info["maxima"]])

# 2. Which kind of motion?  (equilibrium / bounded / separatrix / unbounded)
for x0, v0 in [(1.1, 0.0), (0.0, np.sqrt(4 / 27)), (1.2, 0.3)]:
    print(f"motion through ({x0}, {v0:.4f}):", motion_type(m, x0, v0))

# 3. Period: double precision and 40 digits (exact rational input)
print("T (double)      :", period(m, 1.1, 0.0))
print("T (40 digits)   :", period_mp(Model(1, -4, 3), Fraction(11, 10), 0, dps=40))

# 4. Closed form x = A y / sqrt(1 + lam y^2 + mu y^4) and comparison with DOP853
cf = exact_solution(m, 1.1, 0.0)
sim = simulate(m, 1.1, 0.0, t_end=20.0, n=2001, rtol=1e-12, atol=1e-14)
print("representation  :", cf.representation)
print("max |x_cf - x_num| on [0,20]: %.1e" % np.max(np.abs(cf.x(sim.t) - sim.y[0])))
print("energy drift (DOP853)       : %.1e" % energy_drift(m, sim))

# 5. Failure handling: separatrix and unbounded motion are reported, never approximated
try:
    period(m, 0.0, np.sqrt(4 / 27))
except ValueError as exc:
    print("period on separatrix ->", exc)
try:
    exact_solution(1, 1, -0.5, 0.0, 3.0)          # c < 0, energy far above the barrier
except ClosedFormNotAvailable as exc:
    print("closed form, escaping orbit ->", exc)

# 6. Beyond the conservative case: damped, harmonically forced response (Poincare map)
P = poincare_section(m, 1.0, 0.0, delta=0.1, F=0.05, omega=1.4, n_periods=50, transient=200)
print("Poincare points after transient: %d distinct" % len(np.unique(np.round(P, 6), axis=0)))
