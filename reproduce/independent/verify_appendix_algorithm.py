"""Algorithm 2 of Supplementary Appendix B (Jacobi form, the formulation of NOLExactSolution/exact_solution),
implemented independently from the printed formulas, without importing nonlinear_oscillator_lab.

Part A (identities). For random problems, every family-I and family-II candidate of step 2 satisfies
  - Psi(lambda, 0) = 0 (family I) and Phi = Psi = 0 (family II), with Phi, Psi parsed from the printed LaTeX;
  - S = 1 + lambda x0^2 + mu x0^4 = s, and the general relations (eq:pq)
        p = a - 6 C0 lambda / S,  q = (Phi - 4 C0 mu)/S,  y'(0) = S v0/(1 - mu x0^4)
    reduce to the explicit expressions (eq:famI)-(eq:famII).
Part B (algorithm). Steps 1-8 (roots, candidates, Jacobi cases (a)-(e), solution, acceptance tests
(i)-(iv), selection, time shift, failure) on random problems a,b in [-2,2], c in [-1.5,1.5],
|x0|,|v0| <= 1.2 and on fixed cases, compared with DOP853 (rtol 2.5e-14) on t_j = j t_w/400,
t_w = min(20, 0.9 t_esc); hat e_x as in the paper. Coverage is compared with Algorithm 1.

Jacobi functions with m < 0 or m > 1 and F(phi|m) for m > 1 are reduced to 0 <= m <= 1 by the standard
transformations (DLMF 22.17, 19.7.4-19.7.5), so no library convention outside [0, 1] enters.
Double precision throughout (the packages use 40 digits in steps 1-3; Supplement, Numerical remarks).

Usage:  python -I verify_appendix_algorithm.py [n_random] [seed]
NumPy and SciPy only.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.special import ellipkinc, ellipk

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from appendix_formulas import (POLY, ev, ev_abs, C0_of, jacobi_sncndn, dop853_reference,   # noqa: E402
                               metrics)


# ------------------------------------------------------------------------------ Jacobi, general real m
def sncndn(u, m):
    u = np.asarray(u, dtype=float)
    if 0.0 <= m <= 1.0:
        return jacobi_sncndn(u, m, 1.0 - m)
    if m < 0.0:                                   # DLMF 22.17.6: parameter mu = -m/(1-m) in (0,1)
        k = -m
        mu = k / (1 + k)
        v = u * np.sqrt(1 + k)
        s, c, d = jacobi_sncndn(v, mu, 1 - mu)
        return s / (d * np.sqrt(1 + k)), c / d, 1 / d
    s, c, d = jacobi_sncndn(u * np.sqrt(m), 1 / m, 1 - 1 / m)   # DLMF 22.17.4: m > 1
    return s / np.sqrt(m), d, c


def F(phi, m):
    """Incomplete elliptic integral of the first kind F(phi | m), real m (m sin^2 phi <= 1)."""
    if m <= 1.0:
        return float(ellipkinc(phi, m))
    th = np.arcsin(min(1.0, np.sqrt(m) * abs(np.sin(phi))))      # DLMF 19.7.4
    return float(np.sign(phi) * ellipkinc(th, 1 / m) / np.sqrt(m))


def K(m):
    return float(ellipk(m)) if m < 1 else float(ellipk(1 / m) / np.sqrt(m))


def cn_inv(v, m):
    return F(np.arccos(np.clip(v, -1, 1)), m)


def sn_inv(v, m):
    return F(np.arcsin(np.clip(v, -1, 1)), m)


def dn_inv(v, m):
    return F(np.arcsin(np.sqrt(np.clip((1 - v * v) / m, 0, 1))), m)


# ------------------------------------------------------------------------------ step 1-2
def roots_ordered(a, b, c, h):
    if c == 0:
        r = list(np.roots([h, -a, -b / 2])) + [0.0] if h != 0 else list(np.roots([-a, -b / 2])) + [0.0, 0.0]
    else:
        co = [h, -a, -b / 2, -c / 3]
        while abs(co[0]) < 1e-300:
            co.pop(0)
        r = list(np.roots(co))
        while len(r) < 3:
            r.append(np.inf)                      # root at infinity when h = 0 (excluded below)
    r = [complex(z) for z in r]
    real = sorted([z.real for z in r if abs(z.imag) <= 1e-12 * max(1, abs(z)) and np.isfinite(z.real)],
                  reverse=True)
    cplx = [z for z in r if abs(z.imag) > 1e-12 * max(1, abs(z))]
    others = [z for z in r if not np.isfinite(z.real)]
    return [complex(x) for x in real] + cplx + others, len(real)


def candidates(a, b, c, x0, v0):
    C0 = C0_of(a, b, c, x0, v0)
    h = 2 * C0
    X = x0 * x0
    P, nreal = roots_ordered(a, b, c, h)
    out = []
    for j in range(nreal):
        Pj = P[j].real
        kl = [P[i] for i in range(3) if i != j]
        if not all(np.isfinite(z.real) for z in kl):
            sig = rho = None
        else:
            sig, rho = kl[0] + kl[1], kl[0] * kl[1]
        if abs(1 - Pj * X) > 1e-14:                                       # family I
            s = 1 / (1 - Pj * X)
            out.append(dict(fam="I", P=Pj, s=s, lam=s * Pj, mu=0.0, p=a - 3 * h * Pj,
                            q=s * (b + 4 * a * Pj - 6 * h * Pj * Pj), yd0=s * v0))
        if sig is not None and abs(sig.imag) < 1e-12 * max(1, abs(sig)) and abs(rho.imag) < 1e-12 * max(1, abs(rho)):
            sig, rho = sig.real, rho.real
            Dl = 1 - sig * X + rho * X * X
            if Dl >= 0:
                s = 4 / (2 - sig * X + 2 * np.sqrt(Dl))
                lam, mu = s * sig / 2, s * s * (sig * sig - 4 * rho) / 16
                if abs(1 - mu * X * X) > 1e-14 and mu != 0:
                    out.append(dict(fam="II", P=Pj, s=s, lam=lam, mu=mu, p=h * (Pj - sig / 2),
                                    q=-h * s * (sig * sig - 4 * rho) / 8, yd0=s * v0 / (1 - mu * X * X),
                                    sig=sig, rho=rho))
    return out, C0, h


# ------------------------------------------------------------------------------ step 3: Jacobi cases
def duffing_jacobi(p, q, y0, yd0):
    """Returns (case, y(t) callable, t_star, bounded)."""
    hy = yd0 * yd0 + p * y0 * y0 + 0.5 * q * y0 ** 4
    Y = y0 * y0
    e1 = 1.0 if yd0 >= 0 else -1.0
    sg0 = np.sign(y0) if y0 != 0 else 1.0
    if q == 0:                                                           # (e)
        if p > 0:
            w = np.sqrt(p)
            return "e", (lambda t: y0 * np.cos(w * t) + yd0 / w * np.sin(w * t)), np.inf, True
        if p < 0:
            w = np.sqrt(-p)
            return "e", (lambda t: y0 * np.cosh(w * t) + yd0 / w * np.sinh(w * t)), np.inf, False
        return "e", (lambda t: y0 + yd0 * t), np.inf, False
    disc = p * p + 2 * q * hy
    z = []
    if disc >= 0:
        sq = np.sqrt(disc)
        z = [(-p + sq) / q, (-p - sq) / q]
        if yd0 == 0:
            i = int(np.argmin([abs(zz - Y) for zz in z]))
            z[i] = Y
    pos = sorted(zz for zz in z if zz > 0)
    if hy > 0 and pos and Y <= pos[0] * (1 + 1e-14):                     # (a)
        Z = pos[0]
        om = p + q * Z
        m = q * Z / (2 * om)
        c1 = -e1 * cn_inv(y0 / np.sqrt(Z), m)
        so = np.sqrt(om)
        return "a", (lambda t: np.sqrt(Z) * sncndn(so * t + c1, m)[1]), np.inf, True
    if q > 0 and hy != 0 and len(z) == 2:                                 # (b)
        B2, b2 = max(z), min(z)
        B = np.sqrt(B2)
        m = 1 - b2 / B2
        om = q * B2 / 2
        e2 = -1.0 if y0 * yd0 > 0 else 1.0
        c1 = e2 * dn_inv(min(abs(y0) / B, 1.0), m)
        so = np.sqrt(om)
        return "b", (lambda t: sg0 * B * sncndn(so * t + c1, m)[2]), np.inf, True
    if q > 0 and hy == 0 and p < 0:                                       # (b) separatrix
        al = np.sqrt(-2 * p / q)
        e2 = -1.0 if y0 * yd0 > 0 else 1.0
        c1 = e2 * np.arccosh(al / abs(y0))
        return "b", (lambda t: sg0 * al / np.cosh(np.sqrt(-p) * t + c1)), np.inf, True
    k = -q
    if q < 0 and hy > 0 and not pos:                                      # (c)
        R = np.sqrt(2 * hy / k)
        m = 0.5 + p / (2 * k * R)
        om = 2 * k * R
        c1 = sg0 * cn_inv((R - Y) / (R + Y), m)
        so = np.sqrt(om)

        def yc(t):
            s, c, _ = sncndn(c1 + e1 * so * t, m)
            return np.sqrt(R) * s / (1 + c)
        return "c", yc, (2 * K(m) - e1 * c1) / so, False
    if q < 0:                                                             # (d)
        ua_c = [zz for zz in pos if zz <= Y * (1 + 1e-14)]
        if not ua_c:
            return None
        ua = max(ua_c)
        below = sorted([zz for zz in (z + [0.0]) if zz < ua * (1 - 1e-15)], reverse=True)
        if len(below) < 2:
            return None
        r, s_ = below[0], below[1]
        m = (r - s_) / (ua - s_)
        om = k * (ua - s_) / 2
        e3 = 1.0 if y0 * yd0 >= 0 else -1.0
        arg = np.sqrt(max((Y - ua) / (Y - r), 0.0)) if Y != r else 0.0
        c1 = e3 * sn_inv(arg, m)
        so = np.sqrt(om)

        def yd(t):
            sn, cn, _ = sncndn(c1 + so * t, m)
            return sg0 * np.sqrt(np.maximum(ua - r * sn * sn, 0.0)) / cn
        return "d", yd, (K(m) - c1) / so, False
    return None


# ------------------------------------------------------------------------------ steps 4-8
def x_of(cand, case, yfun, x0):
    A = np.sqrt(cand["s"])
    lam, mu = cand["lam"], cand["mu"]
    signed = case in ("a", "e") or (lam == 0 and mu == 0)
    sg = np.sign(x0) if x0 != 0 else 1.0

    def x(t):
        y = yfun(np.asarray(t, dtype=float))
        yt = y if signed else sg * np.abs(y)
        with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
            return A * yt / np.sqrt(1 + lam * y * y + mu * y ** 4), 1 + lam * y * y + mu * y ** 4, y
    return x


def accept(x, a, b, c, x0, v0, t_star, lam, mu):
    tw = 6.0
    if lam == 0 and mu == 0 and np.isfinite(t_star):
        tw = min(6.0, 0.9 * t_star)
    ti = np.arange(1, 401) * tw / 400
    xv, D, y = x(ti)
    fin = np.isfinite(y) & (np.abs(y) < 1e6)
    if not np.any(fin) or np.any(D[fin] <= 0):                                    # (i)
        return False
    x0c = x(np.array([0.0]))[0][0]
    if not abs(x0c - x0) <= 1e-8:                                                 # (ii)
        return False
    dl = 1e-3
    xdd = -(a * x0 + b * x0 ** 3 + c * x0 ** 5)
    pred = v0 * dl + xdd * dl * dl / 2
    xdl = x(np.array([dl]))[0][0]
    if not abs(xdl - x0 - pred) <= 0.05 * abs(pred):                              # (iii)
        return False
    hh = 1e-5                                                                     # (iv)
    xp = x(ti + hh)[0]
    xm = x(ti - hh)[0]
    xd = (xp - xm) / (2 * hh)
    C0 = C0_of(a, b, c, x0, v0)
    X2 = xv * xv
    Pc = 2 * C0 - a * X2 - b / 2 * X2 ** 2 - c / 3 * X2 ** 3
    ok = np.abs(xd * xd - Pc) <= 1e-3 * (1 + np.abs(Pc))
    ok = ok | ~fin
    return bool(np.all(np.isfinite(xv[fin])) and np.all(ok))


def algorithm2(a, b, c, x0, v0, _shift=True):
    cands, C0, h = candidates(a, b, c, x0, v0)
    accepted = []
    for cd in cands:
        if not np.isfinite(cd["s"]) or cd["s"] <= 0:
            continue
        r = duffing_jacobi(cd["p"], cd["q"], x0, cd["yd0"])
        if r is None:
            continue
        case, yfun, t_star, bnd = r
        try:
            x = x_of(cd, case, yfun, x0)
            if accept(x, a, b, c, x0, v0, t_star, cd["lam"], cd["mu"]):
                accepted.append((cd, case, x, t_star, bnd))
        except (ValueError, FloatingPointError, ZeroDivisionError):
            continue
    for item in accepted:                                                         # step 6
        if item[4]:
            return item, 0.0
    if accepted:
        return accepted[0], 0.0
    if _shift and v0 == 0:                                                        # step 7
        X = x0 * x0
        co = [-c / 3, -b / 2, -a, 2 * C0]
        while len(co) > 1 and abs(co[0]) < 1e-300:
            co.pop(0)
        rr = [z.real for z in np.roots(co) if abs(z.imag) < 1e-12 and z.real > 0 and abs(z.real - X) > 1e-12 * max(1, X)]
        Pc = lambda u: 2 * C0 - a * u - b / 2 * u * u - c / 3 * u ** 3
        rr = [u for u in rr if Pc((u + X) / 2) > 0]
        if rr:
            uo = min(rr, key=lambda u: abs(u - X))
            x0s = (np.sign(x0) or 1.0) * np.sqrt(uo)
            cands2, _, _ = candidates(a, b, c, x0s, 0.0)
            for cd in cands2:
                if not np.isfinite(cd["s"]) or cd["s"] <= 0:
                    continue
                r = duffing_jacobi(cd["p"], cd["q"], x0s, cd["yd0"])
                if r is None or not np.isfinite(r[2]):
                    continue
                case, yfun, t_star, bnd = r
                x = x_of(cd, case, yfun, x0s)
                if accept(x, a, b, c, x0s, 0.0, t_star, cd["lam"], cd["mu"]):
                    return (cd, case, x, t_star, bnd), t_star
    return None, 0.0                                                              # step 8


# ------------------------------------------------------------------------------ part A
def part_A(n=300, seed=7):
    rng = np.random.default_rng(seed)
    w = dict(phi=0.0, psi=0.0, psiI=0.0, S=0.0, p=0.0, q=0.0, yd=0.0)
    for _ in range(n):
        a, b, c = rng.uniform(-2, 2, 3)
        x0, v0 = rng.uniform(-1.2, 1.2, 2)
        cands, C0, h = candidates(a, b, c, x0, v0)
        X = x0 * x0
        for cd in cands:
            lam, mu = cd["lam"], cd["mu"]
            env = dict(a=a, b=b, c=c, x=x0, C=C0, l=lam, m=mu)
            aenv = {k: abs(v) for k, v in env.items()}
            phi = ev(POLY["Phi"], env)
            if cd["fam"] == "II":
                w["phi"] = max(w["phi"], abs(phi) / ev_abs(POLY["Phi"], aenv))
                w["psi"] = max(w["psi"], abs(ev(POLY["Psi"], env)) / ev_abs(POLY["Psi"], aenv))
            else:
                w["psiI"] = max(w["psiI"], abs(ev(POLY["Psi"], env)) / ev_abs(POLY["Psi"], aenv))
            S = 1 + lam * X + mu * X * X
            w["S"] = max(w["S"], abs(S - cd["s"]) / abs(cd["s"]))
            p_pq = a - 6 * C0 * lam / S
            q_pq = (phi - 4 * C0 * mu) / S
            yd_pq = S * v0 / (1 - mu * X * X)
            sc = 1 + abs(a) + abs(6 * C0 * lam / S)
            w["p"] = max(w["p"], abs(p_pq - cd["p"]) / sc)
            scq = (ev_abs(POLY["Phi"], aenv) + abs(4 * C0 * mu)) / abs(S)
            w["q"] = max(w["q"], abs(q_pq - cd["q"]) / max(scq, 1e-300))
            w["yd"] = max(w["yd"], abs(yd_pq - cd["yd0"]) / max(abs(cd["yd0"]), 1e-300))
    print(f"A  identities, {n} random problems (seed {seed}), relative deviations:")
    print(f"   family II: Phi {w['phi']:.1e}, Psi {w['psi']:.1e};  family I: Psi(lambda,0) {w['psiI']:.1e}")
    print(f"   S = s: {w['S']:.1e};  (eq:pq) vs (eq:famI)/(eq:famII): p {w['p']:.1e}, q {w['q']:.1e}, "
          f"y'(0) {w['yd']:.1e}")
    return max(w.values()) < 1e-12


# ------------------------------------------------------------------------------ part B
def part_B(n=200, seed=20261109, extra=()):
    rng = np.random.default_rng(seed)
    probs = list(extra)
    for _ in range(n):
        a, b = rng.uniform(-2, 2, 2)
        c = rng.uniform(-1.5, 1.5)
        x0, v0 = rng.uniform(-1.2, 1.2, 2)
        probs.append((a, b, c, x0, v0))
    try:
        from algorithm_weierstrass import solve as alg1
    except ImportError:
        alg1 = None
    stats = {"bounded": [0, 0, 0, []], "unbounded": [0, 0, 0, []]}
    cases = {}
    shifted = 0
    for a, b, c, x0, v0 in probs:
        ts, xr, t_esc = dop853_reference(a, b, c, x0, v0)
        cls = "bounded" if np.isinf(t_esc) else "unbounded"
        st = stats[cls]
        st[0] += 1
        (item, tsh) = algorithm2(a, b, c, x0, v0)
        if alg1 is not None and alg1(a, b, c, x0, v0) is not None:
            st[2] += 1
        if item is None:
            continue
        cd, case, x, t_star, bnd = item
        xv = x(ts + tsh)[0]
        if not np.all(np.isfinite(xv)):
            continue
        st[1] += 1
        shifted += tsh > 0
        cases[f"{cd['fam']}/{case}"] = cases.get(f"{cd['fam']}/{case}", 0) + 1
        st[3].append(metrics(xv, xr)[1])
    print(f"B  Algorithm 2 on {len(probs)} problems (seed {seed}) against DOP853:")
    for cls, (nc, nrep, n1, errs) in stats.items():
        summary = (
            f"hat e_x median {np.median(errs):.1e}, max {np.max(errs):.1e}"
            if errs else "no trajectory errors"
        )
        print(f"   {cls:9s} {nc:4d} cases: Algorithm 2 represents {nrep:4d}, Algorithm 1 {n1:4d};"
              f"  {summary}")
    print(f"   returned candidates by family/Jacobi case: {dict(sorted(cases.items()))};  time-shifted: {shifted}")
    nb = stats["bounded"]
    return nb[1] == nb[0] and (np.max(nb[3]) if nb[3] else 0) < 1e-6


FIXED = [(1, 0.4, 0.2, 0.8, 0.1), (1, -4, 3, 0.3, 0.05), (-1, 1, 0.1, 0.5, 0.3), (1, 1, -0.05, 0.5, 0.4),
         (2, -1, 0.3, -0.7, 0.9), (1, -4, 3, 1.1, 0.05), (-1, 0.5, 0.2, 1.0, 0.0), (1, 2, -1, 0.9, 0.2)]


def main():
    args = sys.argv[1:]
    n = int(args[0]) if len(args) > 0 else 200
    s = int(args[1]) if len(args) > 1 else 20261109
    okA = part_A()
    okB = part_B(n, s, FIXED)
    print("PASS" if okA and okB else "CHECK: see above")
    return okA and okB


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
