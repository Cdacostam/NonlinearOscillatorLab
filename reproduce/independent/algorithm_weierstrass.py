"""Algorithm 1 of Supplementary Appendix B (sextic + universal Weierstrass formula), implemented from the
printed formulas alone, without importing nonlinear_oscillator_lab.

    x'' + a x + b x^3 + c x^5 = 0,   x(0) = x0,  x'(0) = v0.

Steps (numbering as in the Algorithm-1 box):
  1  trivial cases (equilibrium; c = 0 -> x = y from (eq:wp));
  2  real roots lambda of the sextic (eq:sextic), coefficients e_0..e_6 parsed from the printed LaTeX;
     roots by numpy, then Newton-polished in 30-digit decimal arithmetic;
  3  mu from the quadratic Phi(lambda, mu) = 0 (numerically stable formula), kept if Psi = 0 and mu != 0;
     (lambda, mu) is then Newton-polished on Phi = Psi = 0 in 30-digit arithmetic (roots of Psi in mu are
     used as extra starting values: Phi alone is ill-conditioned for small C0 and x0);
     plus the family-I pairs lambda = P_j/(1 - P_j x0^2), mu = 0, for the real roots of (eq:resolvent);
     S = 1 + lambda x0^2 + mu x0^4 > 0, mu x0^4 != 1;
  4  p, q, y'(0) from (eq:pq); g2, g3 from (eq:g2g3);
  5  range V = [v-, v+] of y^2;
  6  selection: a pair is global if 1 + lambda v + mu v^2 > 0 on V; prefer mu != 0 and, among those,
     S closest to 1 (the "+" branch of Section B.1; the "-" branch has S ~ x0^-4 and is ill-conditioned);
  7  x(t) = sqrt(S) y / sqrt(1 + lambda y^2 + mu y^4) (sgn(x0)|y| for family I with v- > 0).

Usage
  python -I algorithm_weierstrass.py [n_random] [seed]
Runs 12 fixed problems, n_random (default 200) random problems with a,b in [-2,2], c in [-1.5,1.5],
|x0|,|v0| <= 1.2, and 14 special cases, and compares with DOP853 (rtol 1e-12) up to t = 40 or |x| = 10.
Error: hat e_x = max_j |x - x_ref| / max_j |x_ref| on 401 points (paper, Eq. metrics). Cases with
hat e_x > 1e-9 against DOP853 are re-checked against a 30-digit Taylor-series reference.
NumPy and SciPy only.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from decimal import Decimal, getcontext

import numpy as np

sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from appendix_formulas import (POLY, ev, ev_abs, coeffs_in, C0_of, wp_invariants, wp_y, x_from_y,
                               dop853_reference, taylor_reference, metrics)

DPS = 30


@dataclass
class Pair:
    lam: float
    mu: float
    family: str            # "I" or "II"
    S: float = 0.0
    p: float = 0.0
    q: float = 0.0
    y1: float = 0.0
    vm: float = 0.0
    vp: float = np.inf
    glob: bool = False
    margin: float = -np.inf


@dataclass
class Representation:
    kind: str              # "equilibrium" | "cubic" | "pair"
    a: float
    b: float
    c: float
    x0: float
    v0: float
    pair: Pair | None = None
    bounded: bool = True
    candidates: int = 0

    def __call__(self, t):
        t = np.asarray(t, dtype=float)
        if self.kind == "equilibrium":
            return np.full_like(t, self.x0)
        if self.kind == "cubic":
            Ny, Dy = wp_y(t, self.a, self.b, self.x0, self.v0)
            with np.errstate(divide="ignore", invalid="ignore"):
                return Ny / Dy
        P = self.pair
        Ny, Dy = wp_y(t, P.p, P.q, self.x0, P.y1)
        signed = not (P.family == "I" and P.vm > 0)
        x, Q = x_from_y(Ny, Dy, P.S, P.lam, P.mu, signed=signed, sgn_x0=np.sign(self.x0) or 1.0)
        x = np.where(Q > 0, x, np.nan)
        return x


# --------------------------------------------------------------------------------------------- step 2
def _polish(coeffs_hi_first, z, steps=6):
    """Newton polishing of a real root in DPS-digit decimal arithmetic."""
    getcontext().prec = DPS + 5
    co = [Decimal(repr(float(c))) for c in coeffs_hi_first]
    x = Decimal(repr(float(z)))
    for _ in range(steps):
        f = Decimal(0)
        df = Decimal(0)
        for c in co:
            df = df * x + f
            f = f * x + c
        if df == 0:
            break
        x -= f / df
    return float(x)


def sextic_real_roots(a, b, c, x0, C0):
    env = dict(a=a, b=b, c=c, x=x0, C=C0)
    e = [ev(POLY[f"e{k}"], env) for k in range(7)]
    sc = [ev_abs(POLY[f"e{k}"], env) for k in range(7)]
    # drop leading coefficients that vanish (x0 = 0 -> cubic)
    while len(e) > 1 and abs(e[-1]) <= 1e-13 * max(sc[len(e) - 1], 1e-300) + 1e-300:
        e.pop()
        sc.pop()
    if len(e) <= 1:
        return []
    hi = e[::-1]
    r = np.roots(hi)
    out = []
    for z in r:
        if abs(z.imag) <= 1e-7 * max(1.0, abs(z)):
            out.append(_polish(hi, z.real))
    return out


# --------------------------------------------------------------------------------------------- step 3
def mu_from_phi(lam, a, b, c, x0, C0):
    env = dict(a=a, b=b, c=c, x=x0, C=C0, l=lam)
    k0, k1, k2 = (coeffs_in(POLY["Phi"], "m", env) + [0.0, 0.0])[:3]
    if abs(k2) < 1e-300:
        return [] if abs(k1) < 1e-300 else [-k0 / k1]
    disc = k1 * k1 - 4 * k2 * k0
    if disc < 0:
        if disc > -1e-12 * k1 * k1:
            disc = 0.0
        else:
            return []
    sq = np.sqrt(disc)
    qq = -0.5 * (k1 + np.copysign(sq, k1))      # stable quadratic formula
    roots = [qq / k2]
    if qq != 0:
        roots.append(k0 / qq)
    return roots


def _dec_eval(terms, env):
    tot = Decimal(0)
    for c, d in terms:
        t = Decimal(c.numerator) / Decimal(c.denominator)
        for v, p in d.items():
            t *= env[v] ** p
        tot += t
    return tot


def _dec_deriv(terms, var):
    out = []
    for c, d in terms:
        p = d.get(var, 0)
        if p:
            dd = dict(d)
            dd[var] = p - 1
            out.append((c * p, dd))
    return out


_DPHI = {v: _dec_deriv(POLY["Phi"], v) for v in "lm"}
_DPSI = {v: _dec_deriv(POLY["Psi"], v) for v in "lm"}


def polish_pair(lam, mu, a, b, c, x0, C0, steps=8):
    """Newton iteration on Phi = Psi = 0 in DPS-digit decimal arithmetic (both unknowns at once).
    Returns (lam, mu, relative residual)."""
    getcontext().prec = DPS + 5
    env = {k: Decimal(repr(float(v))) for k, v in dict(a=a, b=b, c=c, x=x0, C=C0).items()}
    L, M = Decimal(repr(float(lam))), Decimal(repr(float(mu)))
    for _ in range(steps):
        env["l"], env["m"] = L, M
        f1, f2 = _dec_eval(POLY["Phi"], env), _dec_eval(POLY["Psi"], env)
        j11, j12 = _dec_eval(_DPHI["l"], env), _dec_eval(_DPHI["m"], env)
        j21, j22 = _dec_eval(_DPSI["l"], env), _dec_eval(_DPSI["m"], env)
        det = j11 * j22 - j12 * j21
        if det == 0:
            break
        L -= (f1 * j22 - f2 * j12) / det
        M -= (j11 * f2 - j21 * f1) / det
    envf = dict(a=a, b=b, c=c, x=x0, C=C0, l=float(L), m=float(M))
    env["l"], env["m"] = L, M
    r1 = abs(float(_dec_eval(POLY["Phi"], env))) / max(ev_abs(POLY["Phi"], envf), 1e-300)
    r2 = abs(float(_dec_eval(POLY["Psi"], env))) / max(ev_abs(POLY["Psi"], envf), 1e-300)
    return float(L), float(M), max(r1, r2)


def mu_candidates(lam, a, b, c, x0, C0):
    """Step 3: mu from the quadratic Phi(lam, mu) = 0; the real roots in mu of Psi(lam, mu) = 0 are added
    as starting values because Phi alone is ill-conditioned when C0 and x0 are both small."""
    out = list(mu_from_phi(lam, a, b, c, x0, C0))
    env = dict(a=a, b=b, c=c, x=x0, C=C0, l=lam)
    co = coeffs_in(POLY["Psi"], "m", env)[::-1]
    while len(co) > 1 and abs(co[0]) <= 1e-14 * max(map(abs, co)):
        co.pop(0)
    if len(co) > 1:
        out += [z.real for z in np.roots(co) if abs(z.imag) <= 1e-8 * max(1.0, abs(z))]
    return out


def psi_ok(lam, mu, a, b, c, x0, C0, tol=1e-8):
    env = dict(a=a, b=b, c=c, x=x0, C=C0, l=lam, m=mu)
    return abs(ev(POLY["Psi"], env)) <= tol * max(ev_abs(POLY["Psi"], env), 1e-300)


def resolvent_roots(a, b, c, h):
    co = [h, -a, -b / 2, -c / 3]
    while len(co) > 1 and abs(co[0]) < 1e-300:
        co.pop(0)
    if len(co) <= 1:
        return []
    return [z.real for z in np.roots(co) if abs(z.imag) <= 1e-10 * max(1.0, abs(z))]


# --------------------------------------------------------------------------------------------- steps 4-6
def _finish_pair(P: Pair, a, b, c, x0, v0, C0):
    X = x0 * x0
    P.S = 1 + P.lam * X + P.mu * X * X
    # S > 0 (numerically: not a cancellation residue) and mu x0^4 != 1
    if P.S <= 1e-12 * (1 + abs(P.lam) * X + abs(P.mu) * X * X) or abs(1 - P.mu * X * X) < 1e-14:
        return False
    env = dict(a=a, b=b, c=c, x=x0, C=C0, l=P.lam, m=P.mu)
    P.p = a - 6 * C0 * P.lam / P.S
    phi = 0.0 if P.family == "II" else ev(POLY["Phi"], env)   # Phi = 0 for family II by (eq:system)
    P.q = (phi - 4 * C0 * P.mu) / P.S
    P.y1 = P.S * v0 / (1 - P.mu * X * X)
    hy, _, _ = wp_invariants(P.p, P.q, x0, P.y1)
    # step 5: range of v = y^2, interval of {v >= 0 : g(v) >= 0} containing X, g(v) = hy - p v - q v^2/2
    g = lambda v: hy - P.p * v - 0.5 * P.q * v * v
    co = [-0.5 * P.q, -P.p, hy]
    while len(co) > 1 and abs(co[0]) < 1e-300:
        co.pop(0)
    rr = sorted(z.real for z in np.roots(co) if abs(z.imag) < 1e-12 * max(1, abs(z))) if len(co) > 1 else []
    tol = 1e-12 * max(1.0, X)
    if abs(P.y1) < 1e-300 or abs(g(X)) <= 1e-14 * max(1.0, abs(hy)):
        eps = 1e-7 * max(X, 1e-7)
        right_neg, left_neg = g(X + eps) < 0, (X - eps < 0 or g(X - eps) < 0)
        if right_neg and left_neg:
            P.vm = P.vp = X                                    # y stays at an equilibrium
        elif right_neg:
            P.vp = X
            below = [r for r in rr if 0 < r < X - tol]
            P.vm = max(below) if below else 0.0
        else:
            P.vm = X
            above = [r for r in rr if r > X + tol]
            P.vp = min(above) if above else np.inf
    else:
        below = [r for r in rr if 0 < r < X]
        above = [r for r in rr if r > X]
        P.vm = max(below) if below else 0.0
        P.vp = min(above) if above else np.inf
    if P.family == "I" and P.vm == 0.0 and np.isinf(P.vp):
        return False
    # step 6: 1 + lam v + mu v^2 > 0 on V ?
    f = lambda v: 1 + P.lam * v + P.mu * v * v
    pts = [P.vm] + ([P.vp] if np.isfinite(P.vp) else [])
    vals = [f(v) for v in pts]
    ok = all(val > 0 for val in vals)
    if np.isinf(P.vp):
        ok = ok and (P.mu > 0 or (P.mu == 0 and P.lam >= 0))
    if P.mu > 0:
        vv = -P.lam / (2 * P.mu)
        if P.vm < vv < P.vp:
            vals.append(f(vv))
            ok = ok and f(vv) > 0
    P.glob = ok
    P.margin = min(vals)
    return True


def solve(a, b, c, x0, v0) -> Representation | None:
    a, b, c, x0, v0 = map(float, (a, b, c, x0, v0))
    # step 1
    if v0 == 0 and a * x0 + b * x0 ** 3 + c * x0 ** 5 == 0:
        return Representation("equilibrium", a, b, c, x0, v0)
    if c == 0:
        return Representation("cubic", a, b, c, x0, v0)
    C0 = C0_of(a, b, c, x0, v0)
    h = 2 * C0
    X = x0 * x0
    pairs: list[Pair] = []
    # steps 2-3, family II
    seen = []
    for lam0 in sextic_real_roots(a, b, c, x0, C0):
        for mu0 in mu_candidates(lam0, a, b, c, x0, C0):
            if not np.isfinite(mu0):
                continue
            try:
                lam, mu, res = polish_pair(lam0, mu0, a, b, c, x0, C0)
            except (ArithmeticError, ValueError):        # Newton diverged / singular Jacobian
                continue
            if not (np.isfinite(lam) and np.isfinite(mu)) or res > 1e-20 or mu == 0:
                continue
            if any(abs(lam - l2) <= 1e-12 * max(1, abs(lam)) and abs(mu - m2) <= 1e-12 * max(1, abs(mu))
                   for l2, m2 in seen):
                continue
            seen.append((lam, mu))
            pairs.append(Pair(lam, mu, "II"))
    # family I
    for P in resolvent_roots(a, b, c, h):
        if abs(1 - P * X) > 1e-14:
            pairs.append(Pair(P / (1 - P * X), 0.0, "I"))
    good = [P for P in pairs if _finish_pair(P, a, b, c, x0, v0, C0)]
    if not good:
        return None
    glob = [P for P in good if P.glob]
    if glob:
        glob.sort(key=lambda P: (P.family != "II", abs(np.log(P.S))))
        return Representation("pair", a, b, c, x0, v0, glob[0], True, len(good))
    good.sort(key=lambda P: (P.family != "II", abs(np.log(P.S))))
    return Representation("pair", a, b, c, x0, v0, good[0], False, len(good))


# --------------------------------------------------------------------------------------------- test driver
FIXED = [  # (a, b, c, x0, v0)
    (1, 0.4, 0.2, 0.8, 0.1), (1, -4, 3, 0.3, 0.05), (-1, 1, 0.1, 0.5, 0.3), (1, 1, -0.05, 0.5, 0.4),
    (2, -1, 0.3, -0.7, 0.9), (1, -4, 3, 1.1, 0.05), (-1, 0.5, 0.2, 1.0, 0.0), (1, 0, 1, 0.3, 0.2),
    (0.5, -1, 0.6, 0.0, 0.7), (1, 2, -1, 0.9, 0.2), (-2, 1, 0.5, -0.4, 0.6), (1, -2, 1.2, 1.05, -0.1),
]
SPECIAL = [
    (1, 0.5, 0.3, 0.0, 0.8), (-1, 1, 0.2, 0.0, 0.5), (1, 0.5, 0.3, 0.7, 0.0), (-1, 2, 0.4, -0.9, 0.0),
    (0, 1, 0.5, 0.6, 0.3), (0, -1, 0.8, 0.9, 0.1), (1, 0, 0.5, 0.6, 0.4), (-1, 0, 1, 0.7, 0.2),
    (1, 0.5, 0, 0.8, 0.3), (-1, 1, 0, 0.3, 0.2),
    (-1, 1, 0.1, 0.0, np.sqrt(2 * 0.0)), (-1, 0, 1, 1e-3, 0.0),       # equilibria / near-equilibrium
    (-1, 1, 0.0, 0.0, 1e-9), (1, -1, 0.0, 1.0, 1 / np.sqrt(2) * 0.999999),  # near separatrices
]


def run_case(a, b, c, x0, v0, t_max=40.0, x_esc=10.0):
    ts, xr, t_esc = dop853_reference(a, b, c, x0, v0, t_max=t_max, x_esc=x_esc, rtol=1e-12, atol=1e-14)
    rep = solve(a, b, c, x0, v0)
    unb = np.isfinite(t_esc)
    if rep is None:
        return dict(rep=False, unbounded=unb, err=np.nan, kind="none")
    x = rep(ts)
    if not np.all(np.isfinite(x)):
        return dict(rep=False, unbounded=unb, err=np.inf, kind=rep.kind)
    _, eh = metrics(x, xr)
    ref = "dop853"
    if eh > 1e-9:                      # re-check against a 30-digit Taylor-series reference
        xt = taylor_reference(a, b, c, x0, v0, ts)
        if np.all(np.isfinite(xt)):
            _, eh = metrics(x, xt)
            ref = "taylor30"
        else:
            ref = "dop853_taylor_nonfinite"
    fam = rep.pair.family if rep.pair else rep.kind
    return dict(rep=True, unbounded=unb, err=eh, kind=fam, ref=ref)


def main(n_random=200, seed=20261109):
    rng = np.random.default_rng(seed)
    probs = [("fixed", *p) for p in FIXED]
    for _ in range(n_random):
        a, b = rng.uniform(-2, 2, 2)
        c = rng.uniform(-1.5, 1.5)
        x0, v0 = rng.uniform(-1.2, 1.2, 2)
        probs.append(("random", a, b, c, x0, v0))
    probs += [("special", *p) for p in SPECIAL]
    rows = []
    for grp, a, b, c, x0, v0 in probs:
        r = run_case(a, b, c, x0, v0)
        r.update(group=grp, a=a, b=b, c=c, x0=x0, v0=v0)
        rows.append(r)
    for grp in ("fixed+random", "special"):
        sel = [r for r in rows if (r["group"] != "special") == (grp == "fixed+random")]
        bnd = [r for r in sel if not r["unbounded"]]
        unb = [r for r in sel if r["unbounded"]]
        errs = np.array([r["err"] for r in sel if r["rep"]])
        print(f"{grp}: {len(sel)} problems ({len(bnd)} bounded, {len(unb)} unbounded); "
              f"represented {sum(r['rep'] for r in sel)} "
              f"(bounded {sum(r['rep'] for r in bnd)}, unbounded {sum(r['rep'] for r in unb)})")
        if len(errs):
            print(f"   hat e_x: median {np.median(errs):.2e}  p95 {np.percentile(errs, 95):.2e}  max {errs.max():.2e}")
        for r in sel:
            if not r["rep"] or r["err"] > 1e-6:
                print(f"   !! {r['group']} a={r['a']:.4g} b={r['b']:.4g} c={r['c']:.4g} x0={r['x0']:.4g} "
                      f"v0={r['v0']:.4g} unbounded={r['unbounded']} rep={r['rep']} kind={r['kind']} err={r['err']:.2e} ref={r.get('ref')}")
    return rows


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    s = int(sys.argv[2]) if len(sys.argv) > 2 else 20261109
    rows = main(n, s)
    if any(not row["rep"] or row.get("ref") == "dop853_taylor_nonfinite" for row in rows):
        raise SystemExit("Some cases lacked a representation or a finite Taylor reference.")
