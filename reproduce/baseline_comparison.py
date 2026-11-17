"""Accuracy, run-time and memory of NonlinearOscillatorLab versus SciPy baselines.

Usage:  python reproduce/baseline_comparison.py

Part A - periods (8 orbit families x 40 energies = 320 periods, the grid of the paper,
         approaching the separatrices to 1e-6 relative distance):
  NOL period      factored, singularity-free quadrature (double precision)
  NOL period_mp   same integrand, mpmath tanh-sinh, 30 digits
  SciPy quad      direct quadrature of 1/sqrt(2(E-V)) between the turning points
                  (scipy.integrate.quad, epsabs=1e-14, epsrel=1e-13, limit=400)
  SciPy QAWS      quad with weight='alg' (x-x_-)^(-1/2)(x_+-x)^(-1/2), remaining factor
                  evaluated directly inside the interval, with analytic values at the endpoints;
                  same tolerances
  SciPy QAWS DD   same weight and limits, smooth factor via divided differences of V
  (all quadrature baselines use the package's turning points, isolating the quadrature)
  SciPy events    solve_ivp(DOP853, rtol=1e-12, atol=1e-14) with v = 0 events;
                  T = time between the 1st and 3rd turning point
  SciPy events*   the same with SciPy's default tolerances (rtol=1e-3, atol=1e-6)
  Reference: period_mp at 40 digits.  Metric: eps_T = |T - T_ref| / T_ref,
  evaluated at 50 digits from the native period outputs.

Part B - trajectories on [0, 20] for the nine benchmark orbits of Table 3:
  NOL exact_solution (closed form, Python port) versus solve_ivp DOP853 at
  rtol = 1e-10 and 1e-13.  Reference: 25-digit mpmath Taylor solution.
  Metric: e_x = max_t |x(t) - x_ref(t)| on 401 uniform points.

Timings: median of 5 repetitions (one for period_mp at 30 digits), serial.
Memory: peak Python heap allocation of one call (tracemalloc); allocations inside compiled
        (C/Fortran) code are not recorded, so memory values are indicative only.
Outputs: reproduce/outputs/baseline_periods.csv, baseline_trajectories.csv,
         baseline_summary.json
"""
from __future__ import annotations

import csv
import json
import platform
import sys
import time
import tracemalloc
import warnings
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp, IntegrationWarning

ROOT = Path(__file__).resolve().parents[1]                       # package root
OUTPUTS = Path(__import__("os").environ.get("NOL_OUTPUTS", ROOT / "reproduce" / "outputs"))
sys.path.insert(0, str(ROOT / "code_py" / "src"))
sys.path.insert(0, str(ROOT / "reproduce"))
from nonlinear_oscillator_lab import Model, period, period_mp, turning_points, exact_solution  # noqa
from closed_form_campaign import mp_reference  # noqa

OUT = OUTPUTS
REPS = 5


def approach(es, n=40, lo=0.02, hi=6.0):
    return es * (1.0 - 10.0 ** -np.linspace(lo, hi, n))


FAMILIES = {
    "hardening": (1.0, 1.0, 0.0, 0.0, np.linspace(0.05, 5.0, 40)),
    "softening": (1.0, -0.5, 0.0, 0.0, approach(0.5)),
    "dw_inner": (-1.0, 1.0, 0.0, 1.0, -0.25 * 10.0 ** -np.linspace(0.005, 6.0, 40)),
    "dw_outer": (-1.0, 1.0, 0.0, 0.0, 10.0 ** np.linspace(-6.0, 0.7, 40)),
    "cq_demo": (1.0, 0.4, 0.2, 0.0, np.linspace(0.05, 5.0, 40)),
    "tw_center": (1.0, -4.0, 3.0, 0.0, approach(2.0 / 27.0)),
    "tw_side": (1.0, -4.0, 3.0, 1.0, approach(2.0 / 27.0)),
    "tw_outer": (1.0, -4.0, 3.0, 0.0, 2.0 / 27.0 + 10.0 ** np.linspace(-6.0, 0.0, 40)),
}


# ------------------------------------------------------------------ baselines
def scipy_quad_period(m, x0, v0):
    E = float(m.energy(x0, v0))
    r = turning_points(m, E)
    l, rr = [(p, q) for p, q in zip(r[:-1], r[1:])
             if E - float(m.potential((p + q) / 2)) > 0 and p - 1e-12 <= x0 <= q + 1e-12][0]
    f = lambda x: 1.0 / np.sqrt(max(2 * (E - float(m.potential(x))), 1e-300))
    return 2 * quad(f, l, rr, epsabs=1e-14, epsrel=1e-13, limit=400)[0]


def scipy_qaws_period(m, x0, v0, *, divided_difference=False):
    """quad with the algebraic endpoint weight (x-x_-)^(-1/2) (x_+-x)^(-1/2) (QUADPACK QAWS);
    the smooth factor uses analytic endpoint values, with either direct evaluation
    in the interior or the polynomial second divided difference of V."""
    E = float(m.energy(x0, v0))
    r = turning_points(m, E)
    l, rr = [(p, q) for p, q in zip(r[:-1], r[1:])
             if E - float(m.potential((p + q) / 2)) > 0 and p - 1e-12 <= x0 <= q + 1e-12][0]

    left_force, right_force = float(m.force(l)), float(m.force(rr))
    if left_force <= 0 or right_force >= 0:
        raise ValueError("QAWS requires two simple turning points")
    g_left = np.sqrt((rr - l) / (2 * left_force))
    g_right = np.sqrt((rr - l) / (-2 * right_force))

    def g(x):
        if x == l:
            return g_left
        if x == rr:
            return g_right
        if divided_difference:
            # V[l,x,r] for even monomials. The h_j are the complete homogeneous
            # polynomials in the three nodes, given by their generating recurrence.
            e1 = l + x + rr
            e2 = l*x + l*rr + x*rr
            e3 = l*x*rr
            h2 = e1*e1 - e2
            h3 = e1*h2 - e2*e1 + e3
            h4 = e1*h3 - e2*h2 + e3*e1
            q = m.a/2 + m.b*h2/4 + m.c*h4/6
            if q <= 0:
                raise ValueError("Nonpositive interior divided difference in QAWS")
            return 1 / np.sqrt(2 * q)
        d = 2 * (E - float(m.potential(x)))
        if d <= 0:
            raise ValueError("Nonpositive interior energy gap in QAWS")
        return np.sqrt((x - l) * (rr - x) / d)

    return 2 * quad(g, l, rr, weight="alg", wvar=(-0.5, -0.5), epsabs=1e-14, epsrel=1e-13, limit=400)[0]


def scipy_event_period(m, x0, v0, rtol=1e-12, atol=1e-14):
    f = lambda t, z: [z[1], float(m.force(z[0]))]
    # turning points are the zeros of v; an event at t ~ 0 (start at a turning
    # point) is discarded, after which the start state recurs at the 2nd event
    ev = lambda t, z: z[1]
    ev.terminal = 4
    sol = solve_ivp(f, (0, 1e4), [x0, v0], method="DOP853", rtol=rtol, atol=atol, events=ev)
    te = sol.t_events[0]
    te = te[te > 1e-8]
    if v0 == 0:
        return float(te[1]) if len(te) >= 2 else np.nan
    return float(te[2] - te[0]) if len(te) >= 3 else np.nan


def timed(fn, *args, reps=REPS):
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter(); out = fn(*args); ts.append(time.perf_counter() - t0)
    return out, float(np.median(ts))


def peak_mem(fn, *args):
    tracemalloc.start(); fn(*args); peak = tracemalloc.get_traced_memory()[1]; tracemalloc.stop()
    return peak


METHODS = {
    "NOL period": lambda m, x, v: period(m, x, v),
    "NOL period_mp(30)": lambda m, x, v: period_mp(m, x, v, dps=30, tol=1e-25),
    "SciPy quad": scipy_quad_period,
    "SciPy QAWS, direct": scipy_qaws_period,
    "SciPy QAWS, divided differences": lambda m, x, v: scipy_qaws_period(m, x, v, divided_difference=True),
    "SciPy events": scipy_event_period,
    "SciPy events (default tol)": lambda m, x, v: scipy_event_period(m, x, v, rtol=1e-3, atol=1e-6),
}


def part_a():
    rows = []
    error_ctx = mp.mp.clone()
    error_ctx.dps = 50
    for fam, (a, b, c, x0, Es) in FAMILIES.items():
        m = Model(a, b, c)
        for E in Es:
            v0 = float(np.sqrt(2.0 * (E - float(m.potential(x0)))))
            Tref = error_ctx.mpf(period_mp(m, x0, v0, dps=40, tol=1e-25))
            for name, fn in METHODS.items():
                reps = 1 if "mp" in name else REPS
                with warnings.catch_warnings(record=True) as recorded:
                    warnings.simplefilter("always", IntegrationWarning)
                    try:
                        T, dt = timed(fn, m, x0, v0, reps=reps)
                        err = float(abs(error_ctx.mpf(T) - Tref) / Tref)
                        failure = "" if np.isfinite(err) else "Nonfinite period estimate"
                    except (ArithmeticError, ValueError, IndexError) as exc:
                        T, dt, err = np.nan, np.nan, np.nan
                        failure = str(exc)
                mem = ""
                if fam == "cq_demo" and E == Es[20] and not failure:
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", IntegrationWarning)
                        mem = peak_mem(fn, m, x0, v0)
                rows.append(dict(family=fam, E=repr(float(E)), method=name,
                                 T=error_ctx.nstr(T, 32) if "mp" in name else T,
                                 T_ref=error_ctx.nstr(Tref, 42), rel_err=err, time_s=dt,
                                 peak_mem_bytes=mem, failure=failure,
                                 quadrature_warning=next((str(w.message) for w in recorded
                                                          if w.category is IntegrationWarning), "")))
        print("done", fam, flush=True)
    return rows


# ------------------------------------------------------------------ trajectories
BENCH = [("harmonic", 4, 0, 0, 0.4, 0), ("hardening", 1, 1, 0, 1, 0), ("softening", 1, -0.5, 0, 0.9, 0),
         ("dw_inner", -1, 1, 0, 1.2, 0), ("dw_outer", -1, 1, 0, 1.6, 0), ("cq_demo", 1, 0.4, 0.2, 0.8, 0.1),
         ("tw_center", 1, -4, 3, 0.3, 0), ("tw_side", 1, -4, 3, 1.1, 0), ("tw_outer", 1, -4, 3, 1.2, 0.3)]


def part_b():
    rows = []
    ts = np.linspace(0, 20, 401)
    for name, a, b, c, x0, v0 in BENCH:
        a, b, c, x0, v0 = map(float, (a, b, c, x0, v0))
        ref = mp_reference(a, b, c, x0, v0, ts, dps=25)

        def cf():
            return exact_solution(a, b, c, x0, v0).x(ts)

        def ode(rtol):
            f = lambda t, z: [z[1], -(a * z[0] + b * z[0] ** 3 + c * z[0] ** 5)]
            return solve_ivp(f, (0, 20), [x0, v0], method="DOP853", t_eval=ts, rtol=rtol,
                             atol=rtol * 1e-2).y[0]
        for meth, fn in [("NOL closed form", cf), ("DOP853 rtol=1e-10", lambda: ode(1e-10)),
                         ("DOP853 rtol=1e-13", lambda: ode(1e-13))]:
            x, dt = timed(fn)
            rows.append(dict(case=name, method=meth, max_abs_err=float(np.max(np.abs(x - ref))),
                             time_s=dt, peak_mem_bytes=peak_mem(fn)))
        print("done", name, flush=True)
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    A = part_a()
    with open(OUT / "baseline_periods.csv", "w", newline="") as f:
        w = csv.DictWriter(f, list(A[0])); w.writeheader(); w.writerows(A)
    B = part_b()
    with open(OUT / "baseline_trajectories.csv", "w", newline="") as f:
        w = csv.DictWriter(f, list(B[0])); w.writeheader(); w.writerows(B)
    summ = {"periods": {}, "trajectories": {}}
    for meth in METHODS:
        R = [r for r in A if r["method"] == meth]
        e = np.array([r["rel_err"] for r in R], float); t = np.array([r["time_s"] for r in R], float)
        mem = [r["peak_mem_bytes"] for r in R if r["peak_mem_bytes"] != ""]
        summ["periods"][meth] = dict(
            n=len(R), failed=int(np.sum(~np.isfinite(e))),
            quadrature_warning=sum(bool(r["quadrature_warning"]) for r in R), err_median=float(np.nanmedian(e)),
            err_max=float(np.nanmax(e)), frac_below_1e10=float(np.mean(np.nan_to_num(e, nan=1) < 1e-10)),
            time_median_ms=1e3 * float(np.nanmedian(t)), peak_mem_kib=(mem[0] / 1024 if mem else None),
            err_max_by_family={fam: float(np.nanmax([r["rel_err"] for r in R if r["family"] == fam]))
                               for fam in FAMILIES})
    for meth in ["NOL closed form", "DOP853 rtol=1e-10", "DOP853 rtol=1e-13"]:
        R = [r for r in B if r["method"] == meth]
        summ["trajectories"][meth] = dict(
            err_median=float(np.median([r["max_abs_err"] for r in R])),
            err_max=float(np.max([r["max_abs_err"] for r in R])),
            time_median_ms=1e3 * float(np.median([r["time_s"] for r in R])),
            peak_mem_kib=float(np.median([r["peak_mem_bytes"] for r in R])) / 1024)
    summ["_meta"] = dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                         mpmath=mp.__version__, platform=platform.platform(), processor=platform.processor(),
                         repeats={name: 1 if "mp" in name else REPS for name in METHODS},
                         trajectory_repeats=REPS, error_dps=50, reference_dps=40)
    json.dump(summ, open(OUT / "baseline_summary.json", "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
