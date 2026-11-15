"""Check of the universal Weierstrass formula (eq:wp) with the invariants (eq:g2g3), Supplementary
Appendix B, for the cubic Duffing problem

    y'' + p y + q y^3 = 0,   y(0) = y0,  y'(0) = y1.

Classes (all with uniform grids t_j = j t_w / 400):
  bounded     random p, q in [-3, 3], y0, y1 in [-1.5, 1.5], orbits that stay bounded; DOP853 reference
              (rtol 1e-12, atol 1e-14) on t_w = 20
  separatrix  p < 0 < q, h_y = 0; exact reference y = sgn(y0) alpha sech(sqrt(-p) t + c1): (a) formula
              with h_y = 0 exactly on t_w = 20; (b) h_y from the rounded data, while y > 1e-6 alpha
  linear      q = 0, p of either sign; exact reference cos/cosh on t_w = 20
  unbounded   orbits that reach |y| = 1e3; DOP853 reference on t_w = 0.9 t_esc; the first genuine pole
              of (eq:wp) is compared with the DOP853 escape time plus the
              analytic tail sqrt(2/(-q))/1e3 beyond |y| = 1e3
              Zeros of the denominator at which the numerator also vanishes (removable, y finite) are
              counted and reported separately.
  mutation    every sign in (eq:wp) flipped in turn: the error must then exceed the printed formula's
              largest error by orders of magnitude
Metric: hat e_y = max_j |y - y_ref| / max_j |y_ref|, the paper's amplitude-normalised metric applied to y.

wp is evaluated from its classical Jacobi representations with the roots of 4z^3 - g2 z - g3 taken in
closed form (z = p/6, z = -p/12 +- sqrt(-q h_y/2)/2, see appendix_formulas.wp_roots), not from a library
WeierstrassP. Usage:  python -I verify_weierstrass.py [n_draws] [seed]      NumPy and SciPy only.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from appendix_formulas import wp_y, wp_invariants, _wp_parts, metrics   # noqa: E402


def duffing_ref(p, q, y0, y1, t_max=20.0, y_esc=1e3):
    f = lambda t, u: [u[1], -p * u[0] - q * u[0] ** 3]

    def esc(t, u):
        return abs(u[0]) - y_esc
    esc.terminal = True
    pre = solve_ivp(f, (0, t_max), [y0, y1], method="DOP853", rtol=1e-12, atol=1e-14, events=esc)
    t_esc = float(pre.t_events[0][0]) if len(pre.t_events[0]) else np.inf
    tw = min(t_max, 0.9 * t_esc)
    ts = np.linspace(0, tw, 401)
    s = solve_ivp(f, (0, tw), [y0, y1], method="DOP853", t_eval=ts, rtol=1e-12, atol=1e-14)
    return ts, s.y[0], t_esc


def y_formula(ts, p, q, y0, y1, hy=None):
    Ny, Dy = wp_y(ts, p, q, y0, y1, hy)
    return Ny / Dy


def y_mutated(ts, p, q, y0, y1, k):
    """(eq:wp) with the sign of its k-th term flipped (k = 0..5)."""
    hy, _, _ = wp_invariants(p, q, y0, y1)
    u, w = _wp_parts(ts, p, q, hy)
    s = [1.0] * 6
    s[k] = -1.0
    kk = (p + 3 * q * y0 * y0) / 12
    num = s[0] * y1 * w + s[1] * y0 * (p + q * y0 * y0) * (u + s[2] * kk * u * u) + s[3] * 0.5 * q * y0 * y1 * y1 * u * u
    den = 2 * (1 + s[4] * kk * u) ** 2 + s[5] * 0.25 * q * y1 * y1 * u * u
    return y0 - num / den


def blowup_time(p, q, y0, y1, t_max):
    """First pole of y on (0, t_max]: a sign change of the denominator Dy of (eq:wp) at which y itself
    diverges. Returns (t_pole, n_removable), n_removable = sign changes of Dy passed on the way at which the
    numerator vanishes too (removable 0/0 points of (eq:wp): y stays finite there)."""
    ts = np.linspace(1e-9, t_max, 20001)
    _, Dy = wp_y(ts, p, q, y0, y1)
    idx = np.nonzero(np.sign(Dy[1:]) != np.sign(Dy[:-1]))[0]
    n_rem = 0
    for i in idx:
        lo, hi = ts[i], ts[i + 1]
        flo = wp_y(np.array([lo]), p, q, y0, y1)[1][0]
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            fm = wp_y(np.array([mid]), p, q, y0, y1)[1][0]
            if np.sign(fm) == np.sign(flo):
                lo, flo = mid, fm
            else:
                hi = mid
        tm = 0.5 * (lo + hi)
        yy = y_formula(np.array([tm - 1e-7, tm + 1e-7]), p, q, y0, y1)
        if np.min(np.abs(yy)) > 1e4:              # y ~ 1/(t - t_pole) on both sides: genuine pole
            return tm, n_rem
        n_rem += 1
    return np.inf, n_rem


def main(n=600, seed=20261109):
    rng = np.random.default_rng(seed)
    res = {"bounded": [], "unbounded": [], "separatrix": [], "linear": []}
    tstar = []
    mut = []
    removable = []
    example = None
    for _ in range(n):
        p, q = rng.uniform(-3, 3, 2)
        y0, y1 = rng.uniform(-1.5, 1.5, 2)
        ts, yr, t_esc = duffing_ref(p, q, y0, y1)
        y = y_formula(ts, p, q, y0, y1)
        if not np.all(np.isfinite(y)):
            continue
        cls = "bounded" if np.isinf(t_esc) else "unbounded"
        res[cls].append(metrics(y, yr)[1])
        if cls == "unbounded":
            tb, nr = blowup_time(p, q, y0, y1, 2 * t_esc + 1)
            removable.append(nr)
            if nr and example is None:
                example = (p, q, y0, y1)
            # DOP853 stops at |y| = 1e3; the remaining time to y = infinity is int_{1e3}^inf dy/(sqrt(-q/2) y^2)
            # to leading order, i.e. sqrt(2/(-q))/1e3 (relative correction O(1e-6))
            t_ref = t_esc + np.sqrt(2.0 / (-q)) / 1e3
            tstar.append(abs(tb - t_ref) / t_ref)
        elif len(mut) < 40:
            errs = []
            for k in range(6):
                ym = y_mutated(ts, p, q, y0, y1, k)
                errs.append(metrics(ym, yr)[1] if np.all(np.isfinite(ym)) else np.inf)
            mut.append(min(errs))
    # separatrix: p < 0 < q, h_y = 0.  (a) the formula's degenerate limit h_y = 0 exactly, t_w = 20;
    # (b) h_y recomputed from the rounded data (|h_y| ~ 1e-16): the orbit is then a nearby periodic or
    # escaping one, so it is compared only while the separatrix solution is above 1e-6 alpha.
    res["separatrix_float"] = []
    for _ in range(100):
        p = -rng.uniform(0.2, 3)
        q = rng.uniform(0.2, 3)
        al = np.sqrt(-2 * p / q)
        y0 = rng.choice([-1, 1]) * rng.uniform(0.05, 0.99) * al
        y1 = rng.choice([-1, 1]) * np.sqrt(max(-p * y0 * y0 - 0.5 * q * y0 ** 4, 0.0))
        ts = np.linspace(0, 20, 401)
        eps2 = -1.0 if y0 * y1 > 0 else 1.0
        c1 = eps2 * np.arccosh(al / abs(y0))
        yr = np.sign(y0) * al / np.cosh(np.sqrt(-p) * ts + c1)
        res["separatrix"].append(metrics(y_formula(ts, p, q, y0, y1, hy=0.0), yr)[1])
        tw = min(20.0, (np.log(2e6) - c1) / np.sqrt(-p))
        ts2 = np.linspace(0, tw, 401)
        yr2 = np.sign(y0) * al / np.cosh(np.sqrt(-p) * ts2 + c1)
        res["separatrix_float"].append(metrics(y_formula(ts2, p, q, y0, y1), yr2)[1])
    # linear: q = 0
    for _ in range(100):
        p = rng.uniform(-3, 3)
        y0, y1 = rng.uniform(-1.5, 1.5, 2)
        tw = 20.0 if p > 0 else min(20.0, 5.0 / np.sqrt(-p))
        ts = np.linspace(0, tw, 401)
        if p > 0:
            w = np.sqrt(p)
            yr = y0 * np.cos(w * ts) + y1 / w * np.sin(w * ts)
        else:
            w = np.sqrt(-p)
            yr = y0 * np.cosh(w * ts) + y1 / w * np.sinh(w * ts)
        res["linear"].append(metrics(y_formula(ts, p, q := 0.0, y0, y1), yr)[1])

    print(f"(eq:wp) against references, seed {seed}:")
    for k, v in res.items():
        v = np.array(v)
        print(f"  {k:10s} {len(v):4d} cases   hat e_y median {np.median(v):.1e}   max {v.max():.1e}")
    tstar = np.array(tstar)
    print(f"  blow-up time: first genuine pole of (eq:wp) vs DOP853 (+ analytic tail beyond |y|=1e3):"
          f" relative difference median {np.median(tstar):.1e}, max {tstar.max():.1e}")
    removable = np.array(removable)
    print(f"  removable zeros of the denominator of (eq:wp) before the pole: in {np.count_nonzero(removable)} of "
          f"{len(removable)} unbounded orbits (numerator and denominator vanish together, y finite)")
    if example is not None:
        pe, qe, y0e, y1e = example
        print(f"    example: p={pe:.6g}, q={qe:.6g}, y0={y0e:.6g}, y1={y1e:.6g}")
    mut = np.array(mut)
    worst_b = max(res["bounded"])
    print(f"  mutation test: with any one of the 6 signs flipped the error is >= {mut.min():.1e} on all "
          f"{len(mut)} bounded test orbits (median of the smallest flip {np.median(mut):.1e}), i.e. >= "
          f"{mut.min() / worst_b:.0e} times the largest error of the printed formula")
    worst = max(max(v) for k, v in res.items())
    ok = worst < 1e-8 and mut.min() > 100 * worst_b and tstar.max() < 1e-6
    print("PASS" if ok else "CHECK: see the maxima above")
    return ok


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    s = int(sys.argv[2]) if len(sys.argv) > 2 else 20261109
    raise SystemExit(0 if main(n, s) else 1)
