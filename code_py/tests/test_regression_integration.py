"""Regression tests against the published benchmark values and the Wolfram
implementation, and integration tests of complete workflows."""
import csv
import subprocess
import sys
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from nonlinear_oscillator_lab import (Model, period, period_mp, exact_solution, simulate,
                                      energy_drift, classify_energy, motion_type)

ROOT = Path(__file__).resolve().parents[1]          # code_py/
REPO = ROOT.parent                                    # package root (data/, reproduce/)
BENCH = list(csv.DictReader(open(REPO / "data" / "wolfram" / "benchmark_periods.csv")))
CAMPAIGN = list(csv.DictReader(open(REPO / "data" / "wolfram" / "closed_form_campaign.csv")))


# --------------------------------------------------------------- regression
@pytest.mark.regression
@pytest.mark.parametrize("row", BENCH, ids=[r["case"] for r in BENCH])
def test_benchmark_period_double_precision(row):
    a, b, c, x0, v0 = (float(Fr(row[k])) for k in ("a", "b", "c", "x0", "v0"))
    assert abs(period(Model(a, b, c), x0, v0) / float(row["T_40digit"]) - 1) < 3e-15


@pytest.mark.regression
@pytest.mark.parametrize("row", BENCH, ids=[r["case"] for r in BENCH])
def test_benchmark_period_arbitrary_precision_matches_wolfram(row):
    """period_mp with exact rational input reproduces the 16 printed digits of the
    40-digit Wolfram Language value."""
    T = period_mp(Model(*(Fr(row[k]) for k in "abc")), Fr(row["x0"]), Fr(row["v0"]), dps=40)
    assert mp.nstr(T, 16, strip_zeros=False) == mp.nstr(mp.mpf(row["T_40digit"]), 16, strip_zeros=False)


@pytest.mark.regression
@pytest.mark.parametrize("row", CAMPAIGN[::5], ids=lambda r: "case")
def test_closed_form_representation_matches_wolfram(row):
    """The Python port selects the same representation as the Wolfram solver
    (40 of the 200 published campaign cases; all 200 in benchmarks/)."""
    s = exact_solution(*(float(row[k]) for k in ("a", "b", "c", "x0", "v0")))
    assert s.representation == row["rep"]


# --------------------------------------------------------------- integration
@pytest.mark.integration
def test_full_workflow_triple_well():
    """classify -> period -> closed form -> simulation -> cross-check, as in docs/tutorial.md."""
    m = Model(1, -4, 3)
    info = classify_energy(m, 1.1, 0.0)
    assert not info["on_separatrix_energy"] and motion_type(m, 1.1, 0.0) == "bounded"
    T = period(m, 1.1, 0.0)
    cf = exact_solution(m, 1.1, 0.0)
    sim = simulate(m, 1.1, 0.0, t_end=2 * T, n=801, rtol=1e-12, atol=1e-14)
    assert np.max(np.abs(cf.x(sim.t) - sim.y[0])) < 1e-9
    assert abs(cf.x(T) - 1.1) < 1e-9                     # closed form is T-periodic
    assert energy_drift(m, sim) < 1e-11


@pytest.mark.integration
@pytest.mark.parametrize("script", ["examples/demo.py", "examples/quickstart.py"])
def test_example_scripts_run(script, tmp_path):
    r = subprocess.run([sys.executable, str(ROOT / script)], capture_output=True, text=True,
                       cwd=tmp_path, timeout=300)
    assert r.returncode == 0, r.stderr


@pytest.mark.integration
@pytest.mark.slow
def test_benchmark_script_runs(tmp_path):
    r = subprocess.run([sys.executable, str(REPO / "reproduce" / "python_benchmarks.py")],
                       capture_output=True, text=True, timeout=300, cwd=tmp_path,
                       env={**__import__("os").environ, "NOL_OUTPUTS": str(tmp_path)})
    assert r.returncode == 0, r.stderr
