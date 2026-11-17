"""Unit tests added in v0.2.0: energy geometry, arbitrary precision, closed forms, forcing."""
import numpy as np
import mpmath as mp
import pytest
from fractions import Fraction as Fr
from scipy.integrate import solve_ivp

from nonlinear_oscillator_lab import (Model, stationary_points, classify_energy, turning_points,
                                      motion_type, period, period_mp, simulate, energy_drift,
                                      exact_solution, duffing_solution, simulate_forced,
                                      poincare_section, energy_balance_residual)
from nonlinear_oscillator_lab.closed_form import _ellipj_real

pytestmark = pytest.mark.unit


def _ref(a, b, c, x0, v0, ts):
    f = lambda t, z: [z[1], -(a * z[0] + b * z[0] ** 3 + c * z[0] ** 5)]
    return solve_ivp(f, (0, ts[-1]), [x0, v0], method="DOP853", t_eval=ts, rtol=2.5e-14, atol=1e-16).y[0]


# ----------------------------------------------------------------- model / geometry
@pytest.mark.parametrize("a,b,c", [(1, 0.4, 0.2), (-1, 1, 0), (1, -4, 3), (2, -1, -0.3)])
def test_force_is_minus_gradient(a, b, c):
    m = Model(a, b, c); x = np.linspace(-1.3, 1.3, 27); h = 1e-6
    dV = (m.potential(x + h) - m.potential(x - h)) / (2 * h)
    assert np.allclose(m.force(x), -dV, atol=1e-7)


def test_energy_definition():
    m = Model(1, 0.4, 0.2)
    assert m.energy(0.8, 0.1) == pytest.approx(0.005 + 0.32 + 0.4 * 0.8 ** 4 / 4 + 0.2 * 0.8 ** 6 / 6)


@pytest.mark.parametrize("a,b,c,expected", [
    (1, 0, 0, [0]), (1, 1, 0, [0]), (-1, 1, 0, [-1, 0, 1]),
    (1, -4, 3, [-1, -3 ** -0.5, 0, 3 ** -0.5, 1]), (1, -2, 1, [-1, 0, 1]),
    (1, 0, -1, [-1, 0, 1])])
def test_stationary_points(a, b, c, expected):
    assert np.allclose(stationary_points(Model(a, b, c)), expected)


def test_stationary_discriminant_and_critical_motion():
    model = Model(1, -2, 1)
    v = np.sqrt(1 / 3)
    geometry = classify_energy(model, 0, v)
    assert geometry["equilibria"] == [-1, 0, 1]
    assert geometry["maxima"] == []
    assert [x for x, _ in geometry["minima"]] == [0]
    assert [x for x, _ in geometry["degenerate"]] == [-1, 1]
    assert motion_type(model, 0, v) == "separatrix"
    assert motion_type(model, 1, 0) == "equilibrium"
    for measure in (period, period_mp):
        with pytest.raises(ValueError, match="Separatrix energy"):
            measure(model, 0, v)
    assert len(stationary_points(Model(np.nextafter(1., 0.), -2, 1))) == 5
    assert len(stationary_points(Model(np.nextafter(1., 2.), -2, 1))) == 1

    assert np.allclose(stationary_points(Model(-1e100, 1e100, 1e-300)), [-1, 0, 1])
    wide = stationary_points(Model(1e100, 0, -1e-300))
    assert np.allclose(wide, [-1e100, 0, 1e100])

    # Non-binary coefficients and an exactly repeated stationary point at x=5/3.
    exact = Model("1.25", "-0.9", Fr(81, 500))
    with pytest.raises(ValueError, match="Separatrix energy"):
        period_mp(exact, 0, "1.075828707279838", dps=40)


def test_representable_wide_barrier_and_explicit_energy_overflow():
    model = Model(1, 0, -1e-220)
    barrier_energy = 1e110/3
    maxima = classify_energy(model, 0, 1)["maxima"]
    assert [x for x, _ in maxima] == pytest.approx([-1e55, 1e55])
    assert [e for _, e in maxima] == pytest.approx([barrier_energy]*2, rel=1e-14)
    assert np.allclose(model.potential([[0, 1e55], [-1e55, 0]]),
                       [[0, barrier_energy], [barrier_energy, 0]], rtol=1e-14)
    outside_range = Model(1e130, 0, -1e-310)
    with pytest.raises(ValueError, match="Stationary energy outside"):
        classify_energy(outside_range, 0, 1)
    with pytest.raises(ValueError, match="Initial energy outside"):
        classify_energy(Model(1, 0, 1), 1e100, 1)


def test_fraction_coordinates_in_potential_and_period():
    model = Model(4, 0, 0)
    x = Fr(2, 5)
    energy = model.potential(0.4)
    assert model.potential(x) == energy
    np.testing.assert_array_equal(model.potential([x, -x]), [energy, energy])
    coords = np.array([[x, -x], [Fr(1, 5), 0]], dtype=object)
    np.testing.assert_array_equal(model.potential(coords), model.potential(coords.astype(float)))
    rational_model = Model(4, 0, Fr(1, 3))
    np.testing.assert_array_equal(rational_model.potential(coords), Model(4, 0, 1/3).potential(coords.astype(float)))
    assert model.energy(x, 0) == energy
    assert classify_energy(model, x, 0) == classify_energy(model, 0.4, 0)
    assert period(model, x, 0) == pytest.approx(np.pi, rel=1e-14)
    with mp.workdps(55):
        assert abs(period_mp(model, x, 0, dps=50) - mp.pi) < mp.mpf("1e-49")


def test_turning_points_are_roots():
    m = Model(1, -4, 3)
    for E in (0.02, 0.07, 0.1, 1.0):
        x = turning_points(m, E)
        assert len(x) in (2, 4, 6)
        assert np.allclose(m.potential(x), E, atol=1e-12)


def test_classify_energy_keys_and_barriers():
    d = classify_energy(Model(1, -4, 3), 0.5, 0.0)
    assert {"energy", "equilibria", "maxima", "minima", "degenerate", "on_separatrix_energy"} <= set(d)
    assert [round(e, 12) for _, e in d["maxima"]] == [round(2 / 27, 12)] * 2
    assert len(d["minima"]) == 3


@pytest.mark.parametrize("model,x0,v0,expected", [
    (Model(1, 0.4, 0.2), 0.8, 0.1, "bounded"),
    (Model(-1, 1), 1.0, 0.0, "equilibrium"),
    (Model(-1, 1), 0.0, 0.0, "equilibrium"),
    (Model(-1, 1), 1.0, 2 ** -0.5, "separatrix"),
    (Model(1, -1), 2.0, 0.0, "unbounded"),
    (Model(1, 1, -0.5), 0.0, 3.0, "unbounded"),
    (Model(0, 0, 0), 0.0, 1.0, "unbounded"),
    (Model(1, -2, 1), 0.0, np.sqrt(1 / 3), "separatrix"),   # inflection equilibria at +-1
])
def test_motion_type(model, x0, v0, expected):
    assert motion_type(model, x0, v0) == expected


# ----------------------------------------------------------------- periods
def test_period_mp_harmonic_40_digits():
    T = period_mp(Model(4), Fr(2, 5), 0, dps=40)
    with mp.workdps(45):
        assert abs(T - mp.pi) < mp.mpf(10) ** -38


@pytest.mark.parametrize("a,b,A", [(1, 1, 1), (-1, 1, Fr(6, 5)), (-1, 1, Fr(8, 5)), (1, Fr(-1, 2), Fr(9, 10))])
def test_period_mp_matches_elliptic_50_digits(a, b, A):
    mp.mp.dps = 60
    a_, b_, A_ = (mp.mpf(Fr(v).numerator) / Fr(v).denominator for v in (a, b, A))
    E = a_ * A_ ** 2 / 2 + b_ * A_ ** 4 / 4
    if a_ < 0 < b_ and E < 0:
        Tex = 2 * mp.ellipk(2 + 2 * a_ / (b_ * A_ ** 2)) / (mp.sqrt(b_ / 2) * A_)
    else:
        w2 = a_ + b_ * A_ ** 2; Tex = 4 * mp.ellipk(b_ * A_ ** 2 / (2 * w2)) / mp.sqrt(w2)
    T = period_mp(Model(a, b), A, 0, dps=50)
    mp.mp.dps = 15
    assert abs(T / Tex - 1) < mp.mpf(10) ** -48


def test_period_mp_string_input_is_exact_decimal():
    T1 = period_mp(Model("1", "0.4", "0.2"), "0.8", "0.1", dps=30)
    T2 = period_mp(Model(1, Fr(2, 5), Fr(1, 5)), Fr(4, 5), Fr(1, 10), dps=30)
    assert abs(T1 - T2) < mp.mpf(10) ** -28


def test_period_double_vs_mp():
    for args in [(Model(1, 0.4, 0.2), 0.8, 0.1), (Model(1, -4, 3), 1.1, 0.0), (Model(1, -4, 3), 1.2, 0.3)]:
        assert abs(period(*args) / float(period_mp(*args, dps=30)) - 1) < 1e-14


def test_period_closes_orbit():
    m = Model(1, -4, 3); T = period(m, 1.2, 0.3)
    s = simulate(m, 1.2, 0.3, t_end=T, n=2, rtol=1e-13, atol=1e-15)
    assert abs(s.y[0, -1] - 1.2) < 1e-9 and abs(s.y[1, -1] - 0.3) < 1e-9


# ----------------------------------------------------------------- closed forms
@pytest.mark.parametrize("m", [-3.7, -0.4, 0.0, 0.25, 0.9, 1.0, 1.3, 4.2])
def test_ellipj_real_parameter_transformations(m):
    u = np.linspace(-2.3, 3.1, 11)
    sn, cn, dn = _ellipj_real(u, m)
    for k, arr in zip(("sn", "cn", "dn"), (sn, cn, dn)):
        ref = np.array([complex(mp.ellipfun(k, float(x), m=m)).real for x in u])
        assert np.allclose(arr, ref, atol=1e-12)


@pytest.mark.parametrize("p,q,y0,y1,case", [
    (1.0, 1.0, 0.5, 0.2, "cn"), (-1.0, 1.0, 1.2, 0.0, "dn"), (-1.0, 1.0, 1.6, 0.0, "cn"),
    (1.0, -1.0, 0.5, 0.0, "cn"), (1.0, -1.0, 0.0, 2.0, "Crosses"), (1.0, -1.0, 1.5, 0.0, "Outside"),
    (2.0, 0.0, 0.3, 0.4, "Linear")])
def test_duffing_solution_branches_satisfy_ode(p, q, y0, y1, case):
    s = duffing_solution(p, q, y0, y1)
    assert s.case.startswith(case)
    tmax = 0.8 * s.blow_up_time if s.blow_up_time else 6.0
    ts = np.linspace(0, tmax, 200)
    f = lambda t, z: [z[1], -(p * z[0] + q * z[0] ** 3)]
    ref = solve_ivp(f, (0, tmax), [y0, y1], method="DOP853", t_eval=ts, rtol=1e-13, atol=1e-15).y[0]
    assert np.max(np.abs(s.y(ts) - ref)) < 1e-8 * (1 + np.max(np.abs(ref)))


@pytest.mark.parametrize("a,b,c,x0,v0", [
    (1, 0.4, 0.2, 0.8, 0.1), (1, -4, 3, 0.3, 0.0), (1, -4, 3, 1.1, 0.0), (1, -4, 3, 1.2, 0.3),
    (-1, 1, 0, 1.2, 0.0), (-1, 1, 0, 1.6, 0.0), (1, 1, 0, 1.0, 0.0), (1, -0.5, 0, 0.9, 0.0)])
def test_exact_solution_benchmark_orbits(a, b, c, x0, v0):
    s = exact_solution(a, b, c, x0, v0)
    ts = np.linspace(0, 20, 401)
    assert np.max(np.abs(s.x(ts) - _ref(a, b, c, x0, v0, ts))) < 1e-10
    assert s.x(0.0) == pytest.approx(x0, abs=1e-12)


def test_exact_solution_model_call_and_attributes():
    s = exact_solution(Model(1, 0.4, 0.2), 0.8, 0.1)
    assert s.representation == "II / y: cn: orbit around 0"
    assert s.family == "II" and s.bounded_y
    assert np.isfinite([s.A, s.lam, s.mu, s.p, s.q]).all()
    assert np.shape(s(np.zeros((2, 3)))) == (2, 3)


# ----------------------------------------------------------------- forced / damped
def test_forced_reduces_to_conservative():
    m = Model(1, 0.4, 0.2)
    a = simulate(m, 0.8, 0.1, t_end=20, n=201, rtol=1e-12, atol=1e-14)
    b = simulate_forced(m, 0.8, 0.1, t_end=20, n=201, rtol=1e-12, atol=1e-14)
    assert np.max(np.abs(a.y[0] - b.y[0])) < 1e-10


def test_linear_forced_damped_steady_state_amplitude():
    w0, d, F, w = 1.0, 0.2, 0.5, 1.3
    P = poincare_section(Model(w0 ** 2), 0.0, 0.0, delta=d, F=F, omega=w, n_periods=5, transient=150)
    amp = np.hypot(P[-1, 0], P[-1, 1] / w)
    assert amp == pytest.approx(F / np.hypot(w0 ** 2 - w ** 2, d * w), rel=1e-8)
    assert np.ptp(P, axis=0).max() < 1e-9            # period-1 response: a single point


def test_energy_balance_residual_small():
    m = Model(-1, 1)
    s = simulate_forced(m, 1.0, 0.0, delta=0.25, F=0.3, omega=1.0, t_end=60, rtol=1e-11, atol=1e-13)
    assert energy_balance_residual(m, s) < 1e-6
    s0 = simulate_forced(m, 1.2, 0.0, t_end=60, rtol=1e-11, atol=1e-13)
    assert energy_balance_residual(m, s0) < 1e-9


@pytest.mark.slow
@pytest.mark.parametrize("a,b,c,Es,x0,dE", [(1, -4, 3, 2 / 27, 1.0, 7e-10), (-1, 1, 0, 0.0, 0.0, 7e-12),
                                            (-1, -0.5, 0.5, 0.0, 0.0, 7e-10)])
def test_exact_solution_near_separatrix(a, b, c, Es, x0, dE):
    """40-digit construction and factored denominator keep e_x small close to separatrices
    (in v0.1.0-style double-precision construction e_x reached 1e-1)."""
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "reproduce"))
    from closed_form_campaign import mp_reference
    m = Model(a, b, c)
    v0 = float(np.sqrt(2 * (Es + dE - float(m.potential(x0)))))
    ts = np.linspace(0, 20, 101)
    s = exact_solution(a, b, c, x0, v0)
    assert np.max(np.abs(s.x(ts) - mp_reference(a, b, c, x0, v0, ts, dps=25))) < 1e-8
