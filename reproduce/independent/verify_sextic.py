"""Checks of the sextic (eq:sextic) of Supplementary Appendix B, with the polynomials parsed from the
printed LaTeX (appendix_formulas.POLY):

  1  (float) the six family-II pairs given by radicals,
         lambda = s sigma/2,  mu = s^2 (sigma^2 - 4 rho)/16,
         s = 4/(2 - sigma x0^2 +- 2 sqrt(1 - sigma x0^2 + rho x0^4)),
     for the three roots P_j of (eq:resolvent) and both signs, satisfy Phi = Psi = 0 and are the roots of
     the printed sextic; every family-I value lambda = P_j/(1 - P_j x0^2) satisfies Psi(lambda, 0) = 0;
  2  (exact, Fraction) the printed sextic is, up to a constant factor (powers of C0, x0 and a number), the
     resultant Res_mu(Phi, Psi) -- i.e. it is the elimination polynomial it is said to be;
  3  (exact, Fraction) the discriminant of the printed sextic equals the printed factorization (eq:disc)
         delta = 3442605166011971360784384 C0^6 x0^24 (3 alpha beta + 4 C0 gamma)^2
                 (12 C0 - 6 x0^2 alpha - 3 x0^4 beta - 2 x0^6 gamma)^2 G^2 W^2.

Usage
  python -I verify_sextic.py [n_float] [n_exact] [--tex DIR]
--tex DIR additionally checks that the embedded LaTeX (_appendix_tex.py) equals DIR/appB_*.tex.
Pure Python + NumPy (exact parts use fractions.Fraction; no SymPy needed).
"""
from __future__ import annotations

import itertools
import os
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from appendix_formulas import POLY, KDELTA, ev, ev_abs, ev_exact   # noqa: E402


# ------------------------------------------------------------------------------------------ helpers
def resolvent_roots(a, b, c, h):
    return np.roots([h, -a, -b / 2, -c / 3])


def family_II(a, b, c, x0, C0):
    X = x0 * x0
    P = resolvent_roots(a, b, c, 2 * C0)
    out = []
    for j in range(3):
        kl = np.delete(P, j)
        sig, rho = kl.sum(), kl.prod()
        rt = np.sqrt(complex(1 - sig * X + rho * X * X))
        for sg in (1, -1):
            s = 4 / (2 - sig * X + 2 * sg * rt)
            out.append((s * sig / 2, s * s * (sig * sig - 4 * rho) / 16, sg))
    return out, P


def poly_mul(p, q):
    r = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        if x:
            for j, y in enumerate(q):
                r[i + j] += x * y
    return r


def poly_add(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def det_fraction(M):
    """Exact determinant (Gaussian elimination over Q)."""
    M = [row[:] for row in M]
    n = len(M)
    d = Fraction(1)
    for i in range(n):
        piv = next((r for r in range(i, n) if M[r][i] != 0), None)
        if piv is None:
            return Fraction(0)
        if piv != i:
            M[i], M[piv] = M[piv], M[i]
            d = -d
        d *= M[i][i]
        for r in range(i + 1, n):
            if M[r][i]:
                f = M[r][i] / M[i][i]
                for k in range(i, n):
                    M[r][k] -= f * M[i][k]
    return d


def det_poly(M):
    """Determinant of a small matrix with polynomial (list) entries, by permutation expansion."""
    n = len(M)
    tot = [Fraction(0)]
    for perm in itertools.permutations(range(n)):
        sgn = 1
        for i in range(n):
            for j in range(i + 1, n):
                if perm[i] > perm[j]:
                    sgn = -sgn
        term = [Fraction(sgn)]
        for i in range(n):
            e = M[i][perm[i]]
            if not any(e):
                term = None
                break
            term = poly_mul(term, e)
        if term is not None:
            tot = poly_add(tot, term)
    return trim(tot)


def sylvester(p, q):
    """Sylvester matrix of p, q (coefficient lists high -> low); entries may be numbers or polynomials."""
    m, n = len(p) - 1, len(q) - 1
    zero = 0 if not isinstance(p[0], list) else [Fraction(0)]
    rows = []
    for i in range(n):
        rows.append([zero] * i + list(p) + [zero] * (n - 1 - i))
    for i in range(m):
        rows.append([zero] * i + list(q) + [zero] * (m - 1 - i))
    return rows


def coeff_lambda_mu(terms, env):
    """Exact coefficients: dict mu-power -> list (in lambda, low -> high) of Fractions."""
    out = {}
    for c, d in terms:
        t = c
        for v, p in d.items():
            if v not in "lm":
                t *= env[v] ** p
        mp, lp = d.get("m", 0), d.get("l", 0)
        row = out.setdefault(mp, [])
        while len(row) <= lp:
            row.append(Fraction(0))
        row[lp] += t
    return out


# ------------------------------------------------------------------------------------------ checks
def check_float(n=300, seed=7):
    rng = np.random.default_rng(seed)
    w = {1: [0.0, 0.0], -1: [0.0, 0.0]}
    w_psiI = 0.0
    sext = []
    for _ in range(n):
        a, b, c = rng.uniform(-2, 2, 3)
        x0, v0 = rng.uniform(-1.2, 1.2, 2)
        C0 = v0 * v0 / 2 + a * x0 ** 2 / 2 + b * x0 ** 4 / 4 + c * x0 ** 6 / 6
        env = dict(a=a, b=b, c=c, x=x0, C=C0)
        e = np.array([ev(POLY[f"e{k}"], env) for k in range(7)], dtype=complex)
        pairs, P = family_II(a, b, c, x0, C0)
        for lam, mu, sg in pairs:
            env2 = dict(env, l=lam, m=mu)
            sc_phi = ev_abs(POLY["Phi"], {k: abs(v) for k, v in env2.items()})
            sc_psi = ev_abs(POLY["Psi"], {k: abs(v) for k, v in env2.items()})
            w[sg][0] = max(w[sg][0], abs(ev(POLY["Phi"], env2)) / sc_phi)
            w[sg][1] = max(w[sg][1], abs(ev(POLY["Psi"], env2)) / sc_psi)
            val = sum(e[k] * lam ** k for k in range(7))
            sc = sum(abs(e[k]) * abs(lam) ** k for k in range(7))
            sext.append(abs(val) / sc)
        for Pj in P:
            lam = Pj / (1 - Pj * x0 * x0)
            env2 = dict(env, l=lam, m=0.0)
            w_psiI = max(w_psiI, abs(ev(POLY["Psi"], env2)) / ev_abs(POLY["Psi"], {k: abs(v) for k, v in env2.items()}))
    sext = np.array(sext)
    print(f"1  float, {n} problems (a,b,c in [-2,2], x0,v0 in [-1.2,1.2], seed {seed}):")
    print(f"   family II, '+' branch (3 per problem): max rel. residual Phi {w[1][0]:.1e}, Psi {w[1][1]:.1e}")
    print(f"   family II, '-' branch (3 per problem): max rel. residual Phi {w[-1][0]:.1e}, Psi {w[-1][1]:.1e}"
          "  (s ~ x0^-4: ill-conditioned in double precision, cf. Section B.1)")
    print(f"   family II values in the printed sextic: median rel. residual {np.median(sext):.1e}, "
          f"95th pct {np.percentile(sext, 95):.1e}, max {sext.max():.1e}")
    print(f"   family I: max rel. residual of Psi(lambda, 0) {w_psiI:.1e}")
    return max(w[1] + [w_psiI]) < 1e-12 and np.median(sext) < 1e-13


def rand_rat(rng, lo=-3, hi=3, den=7):
    return Fraction(int(rng.integers(lo * den, hi * den + 1)), den) or Fraction(1, den)


def check_exact(n=4, seed=11):
    rng = np.random.default_rng(seed)
    ok_res = ok_disc = True
    print(f"2-3  exact rational arithmetic, {n} random rational problems (seed {seed}):")
    for it in range(n):
        a, b, c, x0, C = (rand_rat(rng) for _ in range(5))
        env = dict(a=a, b=b, c=c, x=x0, C=C)
        E = [ev_exact(POLY[f"e{k}"], env) for k in range(7)]           # low -> high
        # --- 2: resultant in mu of Phi and Psi, as a polynomial in lambda
        cphi = coeff_lambda_mu(POLY["Phi"], env)
        cpsi = coeff_lambda_mu(POLY["Psi"], env)
        p_hi = [trim(cphi.get(k, [Fraction(0)])) for k in range(max(cphi), -1, -1)]
        q_hi = [trim(cpsi.get(k, [Fraction(0)])) for k in range(max(cpsi), -1, -1)]
        R = det_poly(sylvester(p_hi, q_hi))
        # R must be a constant multiple of E (as polynomials in lambda)
        k = next(i for i in range(7) if E[i] != 0)
        ratio = R[k] / E[k] if k < len(R) else None
        prop = ratio is not None and len(R) == len(trim(E)) and all(
            R[i] == ratio * E[i] for i in range(len(R)))
        ok_res &= prop
        # --- 3: discriminant of the sextic from the resultant of E and E'
        Ehi = E[::-1]
        dE = [E[i] * i for i in range(1, 7)][::-1]
        res = det_fraction(sylvester(Ehi, dE))
        disc = (-1) ** (6 * 5 // 2) * res / Ehi[0]
        G = ev_exact(POLY["G"], env)
        W = ev_exact(POLY["W"], env)
        printed = (KDELTA * C ** 6 * x0 ** 24 * (3 * a * b + 4 * C * c) ** 2
                   * (12 * C - 6 * x0 ** 2 * a - 3 * x0 ** 4 * b - 2 * x0 ** 6 * c) ** 2 * G ** 2 * W ** 2)
        eq = disc == printed
        ok_disc &= eq
        print(f"   problem {it + 1}: (a,b,c,x0,C0) = ({a}, {b}, {c}, {x0}, {C})")
        print(f"      Res_mu(Phi, Psi) = const * sextic: {prop}   (const = {ratio})")
        print(f"      discriminant == printed factorization: {eq}"
              + ("" if eq else f"   ratio = {float(disc / printed) if printed else 'printed=0'}"))
    return ok_res, ok_disc


def check_tex(d):
    import _appendix_tex as T
    ok = True
    for name in ("appB_PhiPsi", "appB_sextic", "appB_G", "appB_W"):
        same = open(os.path.join(d, name + ".tex")).read() == getattr(T, name.upper())
        ok &= same
        print(f"   {name}.tex identical to embedded copy: {same}")
    return ok


def main():
    args = sys.argv[1:]
    tex = None
    if "--tex" in args:
        i = args.index("--tex")
        tex = args[i + 1]
        del args[i:i + 2]
    nf = int(args[0]) if len(args) > 0 else 300
    ne = int(args[1]) if len(args) > 1 else 4
    ok = True
    if tex:
        print("0  embedded LaTeX vs sources:")
        ok &= check_tex(tex)
    ok &= check_float(nf)
    r, d = check_exact(ne)
    ok &= r and d
    print("PASS" if ok else "CHECK: see above")
    return ok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
