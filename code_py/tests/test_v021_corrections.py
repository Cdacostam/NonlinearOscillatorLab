"""v0.2.1 corrections found by the audit of the online demonstrator (November 2026)."""
import pytest

from nonlinear_oscillator_lab import (Model, period, period_mp, motion_type, exact_solution,
                                      ClosedFormNotAvailable)

MSG = "No candidate representation passed the acceptance checks"


@pytest.mark.regression
@pytest.mark.parametrize("c", [1e-2, 1e-5, 1e-7, 1e-8, 1e-10])
def test_period_wide_well_small_quintic(c):
    # a=1, b=-1, x0=0, v0=1: regular orbit whose turning points grow like (1.5/c)**0.5
    # for small c; v0.2.0 raised "Degenerate turning point" for c <= 1e-8 and was inaccurate for c = 1e-7
    m = Model(1.0, -1.0, c)
    T, ref = period(m, 0.0, 1.0), float(period_mp(m, 0.0, 1.0, dps=30))
    assert abs(T - ref) / ref < 1e-13
    if c == 1e-8:
        assert abs(T - 11.416837670635609) < 1e-12


@pytest.mark.failure
@pytest.mark.parametrize("args, kind", [((1.0, -4.0, 3.0, 0.0, (4 / 27) ** 0.5), "separatrix"),
                                        ((1.0, 0.0, 1.0, 100.0, 1.0), "bounded"),
                                        ((1.0, 0.0, -1.0, 0.0, 1.0), "unbounded")])
def test_failure_message_makes_no_claim_on_the_motion(args, kind):
    assert motion_type(Model(*args[:3]), *args[3:]) == kind
    with pytest.raises(ClosedFormNotAvailable) as exc:
        exact_solution(*args)
    assert str(exc.value) == MSG
