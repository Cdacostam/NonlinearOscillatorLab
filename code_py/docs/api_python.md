# Python API reference (`nonlinear_oscillator_lab`, v0.2.1)

All functions accept plain Python/NumPy numbers. Docstrings contain the same information
(`help(nonlinear_oscillator_lab.period)`).

## Model

| Call | Returns |
|---|---|
| `Model(a=1.0, b=0.0, c=0.0)` | frozen dataclass for $\ddot x+ax+bx^3+cx^5=0$ |
| `m.potential(x)` | $V(x)=ax^2/2+bx^4/4+cx^6/6$ (vectorised) |
| `m.force(x)` | $-V'(x)$ |
| `m.energy(x, v)` | $v^2/2+V(x)$ |

## Energy geometry

| Call | Returns | Notes |
|---|---|---|
| `stationary_points(m, tol=1e-10)` | sorted array of real roots of $V'(x)$ | roots of $c y^2+by+a$, $y=x^2$ |
| `turning_points(m, E, tol=1e-9)` | sorted array of real roots of $V(x)=E$ | Newton-polished |
| `classify_energy(m, x0, v0, tol=1e-9)` | `dict`: `energy`, `equilibria`, `maxima`, `minima`, `degenerate`, `on_separatrix_energy` | `degenerate` lists equilibria with $V''=0$ (e.g. $b^2=4ac$); their level is treated as critical |
| `motion_type(m, x0, v0, tol=1e-9)` | `'equilibrium'`, `'separatrix'`, `'bounded'`, `'unbounded'` | |

## Periods

| Call | Returns | Raises |
|---|---|---|
| `period(m, x0, v0, *, tol=1e-9)` | `float` period (double precision); `inf` at an equilibrium | `ValueError` on separatrix energies, unbounded motion, degenerate turning points |
| `period_mp(m, x0, v0, *, dps=50, tol=1e-9)` | `mpmath.mpf` with about `dps` correct digits | same as `period` |

`period_mp` accepts `str` (exact decimal), `fractions.Fraction`, `int` and `float` (exact binary
value) for the coefficients and initial data. Both use the factored integrand
$T=2\int_{-\pi/2}^{\pi/2}d\theta/\sqrt{2Q(\bar x+h\sin\theta)}$ with
$E-V=(x-x_-)(x_+-x)Q(x)$. Lower `tol` (e.g. `1e-25`) when working closer than $10^{-9}$ to a
separatrix at high precision.

## Closed forms

| Call | Returns |
|---|---|
| `exact_solution(m, x0, v0)` or `exact_solution(a, b, c, x0, v0)` | `ClosedFormSolution`, or raises `ClosedFormNotAvailable` |
| `duffing_solution(p, q, y0, y1)` | `DuffingSolution` for $\ddot y+py+qy^3=0$ |

`ClosedFormSolution` attributes: `representation` (str), `family` (`'I'`/`'II'`), `A`, `lam`,
`mu`, `p`, `q`, `y0`, `y1`, `duffing` (the `DuffingSolution`), `time_shift`, `bounded_y`;
method `x(t)` (also `sol(t)`) evaluates $x=Ay/\sqrt{1+\lambda y^2+\mu y^4}$ for scalar or array `t`.

`DuffingSolution` attributes: `case` (`'cn: orbit around 0'`, `'dn: oscillation inside one well'`,
`'Separatrix (sech)'`, `'Crosses 0 and blows up'`, `'Outside barrier: blows up'`, `'Linear'`,
`'Equilibrium'`), `y(t)`, `m`, `omega`, `c1`, `period`, `blow_up_time`, `amplitude`, `bounded`.

Applicability: see the paper (Section 4.3) and `data/results/closed_form_campaign_py_exactref_summary.json`.
A representation exists for all bounded orbits tested and for finite-time escape beyond a
barrier; it does not exist for orbits exactly on a separatrix level or for most escaping orbits
with $c<0$ above the barrier, where `ClosedFormNotAvailable` is raised.

## Simulation

| Call | Returns |
|---|---|
| `simulate(m, x0, v0, t_end=40., n=2001, rtol=1e-10, atol=1e-12)` | SciPy `OdeResult` (DOP853); `RuntimeError` if the integrator fails (e.g. blow-up) |
| `energy_drift(m, sol)` | $\delta_E=\max_j\lvert E_j-E_0\rvert/(1+\lvert E_0\rvert)$ |

## Damped / forced extension

$\ddot x+\delta\dot x+ax+bx^3+cx^5=F\cos(\omega t+\varphi)$, integrated as the autonomous
system $(x,v,\theta)$ with $\dot\theta=\omega$.

| Call | Returns |
|---|---|
| `simulate_forced(m, x0, v0, *, delta=0, F=0, omega=1, phi=0, t_end=100, n=5001, rtol=1e-10, atol=1e-12)` | `OdeResult` with rows x, v, θ and dense output |
| `poincare_section(m, x0, v0, *, delta, F, omega, phi=0, n_periods=200, transient=100)` | array `(n_periods, 2)` of $(x,v)$ at $t_k=2\pi k/\omega$ |
| `energy_balance_residual(m, sol, samples=20001)` | normalised residual of $E(t)-E(0)-\int_0^t(-\delta v^2+Fv\cos\theta)\,ds$ |

## Exceptions

| Exception | When |
|---|---|
| `ValueError('Separatrix energy: infinite period')` | energy equals a barrier or degenerate-equilibrium energy within `tol` |
| `ValueError('Initial state is not in a verified bounded component')` | unbounded motion |
| `ValueError('Degenerate (multiple) turning point ...')` | $V'(x_\pm)=0$ |
| `ClosedFormNotAvailable` (subclass of `ValueError`) | no real representation passes the checks |
| `RuntimeError` | ODE integrator failure |
