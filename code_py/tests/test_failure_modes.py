"""Documented failure modes: every invalid or non-periodic request fails explicitly."""
import numpy as np
import pytest

from nonlinear_oscillator_lab import (Model, period, period_mp, simulate, exact_solution,
                                      ClosedFormNotAvailable, simulate_forced, poincare_section, motion_type)

pytestmark = pytest.mark.failure


@pytest.mark.parametrize("fn", [period, lambda m, x, v: period_mp(m, x, v, dps=20)])
def test_separatrix_raises(fn):
    with pytest.raises(ValueError, match="Separatrix"):
        fn(Model(-1, 1), 1.0, 2 ** -0.5)


@pytest.mark.parametrize("fn", [period, lambda m, x, v: period_mp(m, x, v, dps=20)])
def test_unbounded_raises(fn):
    with pytest.raises(ValueError, match="bounded component"):
        fn(Model(1, -1), 2.0, 0.0)


def test_triple_well_separatrix_raises():
    m = Model(1, -4, 3)
    v0 = np.sqrt(2 * (2 / 27))
    with pytest.raises(ValueError):
        period(m, 0.0, v0)


def test_inflection_equilibrium_energy_is_critical():
    # b^2 = 4ac: V'(x) = x (x^2 - 1)^2; the level E = V(1) = 1/6 has infinite period
    with pytest.raises(ValueError):
        period(Model(1, -2, 1), 0.0, np.sqrt(1 / 3))


def test_unresolved_period_at_split_inflection_level():
    model = Model(1, -2, np.nextafter(1., np.inf))
    v = np.sqrt(1/3)
    assert motion_type(model, 0, v) == "bounded"
    with pytest.raises(ValueError, match="use period_mp"):
        period(model, 0, v)
    assert float(period_mp(model, 0, v, dps=50)) == pytest.approx(3187.7658725993868, rel=1e-12)


def test_period_does_not_return_unconverged_quadpack_value(monkeypatch):
    from nonlinear_oscillator_lab import core

    def quad_with_roundoff(*args, full_output, **kwargs):
        assert full_output == 1
        return 1.0, 0.1, {}, "Roundoff prevented convergence"

    monkeypatch.setattr(core, "quad", quad_with_roundoff)
    with pytest.raises(ValueError, match="use period_mp"):
        period(Model(4), 0.5, 0)


@pytest.mark.parametrize("fn", [period, lambda m, x, v: period_mp(m, x, v, dps=20)])
def test_equilibrium_gives_infinite_period(fn):
    assert np.isinf(float(fn(Model(1, -4, 3), 1.0, 0.0)))


def test_free_particle_unbounded():
    with pytest.raises(ValueError):
        period(Model(0, 0, 0), 0.0, 1.0)


def test_closed_form_unbounded_escape_raises():
    # c < 0, energy far above the barrier: no real representation passes the checks
    with pytest.raises(ClosedFormNotAvailable):
        exact_solution(1, 1, -0.5, 0.0, 3.0)


def test_closed_form_finite_time_escape_is_represented():
    # softening Duffing beyond the barrier escapes in finite time; the closed form
    # follows the orbit up to the escape (|x| -> inf)
    from scipy.integrate import solve_ivp
    s = exact_solution(1, -1, 0, 2.0, 0.0)
    f = lambda t, z: [z[1], -(z[0] - z[0] ** 3)]
    ts = np.linspace(0, 0.9, 200)                     # escape time ~ 0.99966
    ref = solve_ivp(f, (0, 0.9), [2.0, 0.0], method="DOP853", t_eval=ts, rtol=2.5e-14, atol=1e-16).y[0]
    assert np.max(np.abs(s.x(ts) - ref)) < 1e-10


@pytest.mark.parametrize("kw", [dict(t_end=0), dict(t_end=-1), dict(n=1)])
def test_simulate_invalid_arguments(kw):
    with pytest.raises(ValueError):
        simulate(Model(1), 0.1, 0.0, **kw)


def test_forced_invalid_frequency():
    with pytest.raises(ValueError):
        simulate_forced(Model(1), 0.1, 0.0, F=1.0, omega=0.0)
    with pytest.raises(ValueError):
        poincare_section(Model(1), 0.1, 0.0, omega=-1.0)


def test_simulate_reports_solver_failure():
    # finite-time blow-up of x'' = x^5 (c < 0): the integrator cannot reach t_end
    with pytest.raises(RuntimeError):
        simulate(Model(0, 0, -1), 2.0, 2.0, t_end=10.0)


@pytest.mark.parametrize("model,x0,v0", [
    (Model(1, 0.4, 0.2), 0.0, 1e-8), (Model(-1, 1), -1.0, 1e-8), (Model(-1, 1), 1.0, 1e-8),
    (Model(1, -4, 3), 0.0, 1e-8), (Model(1, 0.4, 0.2), 0.0, 1e-4)])
def test_tiny_oscillation_is_bounded_with_harmonic_period(model, x0, v0):
    """v0.1.0 classified oscillations of amplitude ~1e-8 about a minimum as unbounded
    (absolute tolerance on E - V); v0.2.0 reports them as bounded with the correct period."""
    from nonlinear_oscillator_lab import motion_type, stationary_points
    assert motion_type(model, x0, v0) == "bounded"
    k = model.a + 3 * model.b * x0 ** 2 + 5 * model.c * x0 ** 4
    assert period(model, x0, v0) == pytest.approx(2 * np.pi / np.sqrt(k), rel=1e-7)


def test_start_on_barrier_just_above_separatrix_is_bounded():
    from nonlinear_oscillator_lab import motion_type
    assert motion_type(Model(-1, 1), 0.0, np.sqrt(2e-9)) == "bounded"
