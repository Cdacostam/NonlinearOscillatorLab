"""Damped and harmonically forced oscillator (numerical extension, v0.2.0).

    x'' + delta x' + a x + b x^3 + c x^5 = F cos(omega t + phi)

Following the standard device of making a periodically forced system autonomous
(e.g. Hartman, *Ordinary Differential Equations*, 1964), the phase
``theta = omega t + phi`` is added as a third state variable,

    x' = v,   v' = -delta v - (a x + b x^3 + c x^5) + F cos(theta),   theta' = omega,

so that the extended system is autonomous on the cylinder R^2 x S^1 and the
stroboscopic (Poincare) map is the return map to the section theta = phi (mod 2 pi).

Energy-based classification, periods and closed forms apply only to the
conservative case delta = F = 0.  For the forced/damped case the package offers
integration, the Poincare map and an energy-balance diagnostic.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

__all__ = ["simulate_forced", "poincare_section", "energy_balance_residual"]


def _rhs(model, delta, F, omega):
    a, b, c = model.a, model.b, model.c

    def f(t, z):
        x, v, th = z
        return [v, -delta * v - (a * x + b * x ** 3 + c * x ** 5) + F * np.cos(th), omega]
    return f


def simulate_forced(model, x0, v0, *, delta=0.0, F=0.0, omega=1.0, phi=0.0,
                    t_end=100.0, n=5001, rtol=1e-10, atol=1e-12, method="DOP853"):
    """Integrate the damped/forced IVP in the autonomous form (x, v, theta).

    Returns a SciPy ``OdeResult`` whose ``y`` has rows x, v, theta and which carries
    the attributes ``delta``, ``F``, ``omega``, ``phi`` for later diagnostics.
    With ``delta = F = 0`` the result coincides with :func:`simulate`.
    """
    if t_end <= 0 or n < 2:
        raise ValueError("t_end must be positive and n >= 2")
    if omega <= 0 and F != 0:
        raise ValueError("omega must be positive when F != 0")
    t = np.linspace(0.0, t_end, n)
    res = solve_ivp(_rhs(model, delta, F, omega), (0.0, t_end), [x0, v0, phi],
                    method=method, t_eval=t, rtol=rtol, atol=atol, dense_output=True)
    if not res.success:
        raise RuntimeError(res.message)
    res.delta, res.F, res.omega, res.phi = delta, F, omega, phi
    return res


def poincare_section(model, x0, v0, *, delta=0.0, F=0.0, omega=1.0, phi=0.0,
                     n_periods=200, transient=100, rtol=1e-10, atol=1e-12):
    """Stroboscopic samples (x_k, v_k) at t_k = 2 pi k / omega, k > transient.

    Returns an array of shape (n_periods, 2).  A period-1 response gives one point,
    period-n responses give n points, chaotic responses give a fractal set.
    """
    if omega <= 0:
        raise ValueError("omega must be positive")
    T = 2 * np.pi / omega
    k = np.arange(transient, transient + n_periods + 1)
    sol = solve_ivp(_rhs(model, delta, F, omega), (0.0, k[-1] * T), [x0, v0, phi],
                    method="DOP853", t_eval=k * T, rtol=rtol, atol=atol)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:2, 1:].T


def energy_balance_residual(model, sol, samples=20001):
    """Normalised residual of the energy balance along a forced/damped trajectory.

    For E = v^2/2 + V(x) one has dE/dt = -delta v^2 + F v cos(theta), hence
    E(t) - E(0) = int_0^t (-delta v^2 + F v cos theta) ds.  The function returns
    max_t |E(t) - E(0) - int_0^t (...) ds| / (1 + max_t |E(t)|), evaluated with the
    dense output and the trapezoidal rule on ``samples`` points.  It reduces to the
    energy drift of the conservative case when delta = F = 0.
    """
    t = np.linspace(sol.t[0], sol.t[-1], samples)
    x, v, th = sol.sol(t)
    E = model.energy(x, v)
    p = -sol.delta * v * v + sol.F * v * np.cos(th)
    work = np.concatenate([[0.0], np.cumsum(0.5 * (p[1:] + p[:-1]) * np.diff(t))])
    return float(np.max(np.abs(E - E[0] - work)) / (1 + np.max(np.abs(E))))
