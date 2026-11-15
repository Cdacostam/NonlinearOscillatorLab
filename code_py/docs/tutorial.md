# Tutorial: a complete workflow

The example uses the triple-well potential $V(x)=x^2/2-x^4+x^6/2$ ($a=1,b=-4,c=3$), whose
barriers at $x=\pm1/\sqrt3$ have the energy $2/27$. The full script is
`examples/quickstart.py`.

## 1. Energy geometry

```python
from nonlinear_oscillator_lab import Model, classify_energy, motion_type
m = Model(a=1, b=-4, c=3)
info = classify_energy(m, x0=1.1, v0=0.0)
info["equilibria"]      # [-1, -0.57735, 0, 0.57735, 1]
info["maxima"]          # [(-0.57735, 0.074074), (0.57735, 0.074074)]   barriers (x, V)
motion_type(m, 1.1, 0.0)                  # 'bounded'
motion_type(m, 0.0, (4/27)**0.5)          # 'separatrix'
```

```mathematica
m = NOLModel[1, -4, 3];
NOLStationaryPoints[m]                   (* {-1, -1/Sqrt[3], 0, 1/Sqrt[3], 1}  exact *)
NOLClassifyEnergy[m, 11/10, 0]["Maxima"]
```

## 2. Period

```python
from fractions import Fraction
from nonlinear_oscillator_lab import period, period_mp
period(m, 1.1, 0.0)                                   # 3.460792103704556
period_mp(Model(1, -4, 3), Fraction(11, 10), 0, dps=40)
# 3.460792103704556623922823899094473019154
```

Use exact input (`Fraction`, `int`, or decimal strings such as `"0.4"`) for high precision:
a binary float such as `0.4` is taken at its exact binary value.

```mathematica
NOLPeriod[NOLModel[1, -4, 3], 11/10, 0, WorkingPrecision -> 40]
```

## 3. Closed-form solution

```python
from nonlinear_oscillator_lab import exact_solution, simulate
cf = exact_solution(m, 1.1, 0.0)
cf.representation           # 'I / y: Outside barrier: blows up (time-shifted)'
cf.A, cf.lam, cf.mu         # x = A y / sqrt(1 + lam y^2 + mu y^4)
cf.x(3.0)                   # scalar or numpy array of times
sim = simulate(m, 1.1, 0.0, t_end=20, rtol=1e-12, atol=1e-14)
abs(cf.x(sim.t) - sim.y[0]).max()      # ~ 2e-12
```

The representation label tells which family (I: $\mu=0$, II: $\mu\ne0$) and which Jacobi
branch of the auxiliary Duffing solution $y$ were used. A $y$ that "blows up" is not an error:
$x$ stays bounded because the denominator grows with $y$.

## 4. Failure handling

```python
period(m, 0.0, (4/27)**0.5)       # ValueError: Separatrix energy: infinite period
exact_solution(1, 1, -0.5, 0, 3)  # ClosedFormNotAvailable: no real representation
```

Every non-periodic or unsupported request raises an exception (Python) or returns a
`Failure` object with a message (Wolfram Language). Nothing is silently approximated. See
[troubleshooting](troubleshooting.md).

## 5. Damped and forced response (numerical extension)

```python
from nonlinear_oscillator_lab import simulate_forced, poincare_section, energy_balance_residual
s = simulate_forced(m, 1.0, 0.0, delta=0.1, F=0.05, omega=1.4, t_end=300)
energy_balance_residual(m, s)            # integration check: dE/dt = -delta v^2 + F v cos(theta)
P = poincare_section(m, 1.0, 0.0, delta=0.1, F=0.05, omega=1.4, n_periods=200, transient=200)
```

The forcing is made autonomous with the phase variable $\theta=\omega t+\varphi$; `P` holds the
stroboscopic samples at $t_k=2\pi k/\omega$ (one point for a period-1 response).
