"""Shared, package-independent helpers for the Appendix-B benchmark scripts.

Nothing here imports nonlinear_oscillator_lab. The polynomials Phi, Psi, e_0..e_6, G and W are parsed
directly from the printed LaTeX of Supplementary Appendix B (verbatim copies in _appendix_tex.py), so the
scripts test the formulas as printed, not a retyped version.

Contents
  parse_latex_poly / POLY       printed polynomials as exact (Fraction) term lists
  ev, ev_exact                  evaluation in float / Fraction arithmetic
  wp_y                          Duffing solution y(t) from the universal Weierstrass formula (eq:wp),
                                returned as numerator/denominator so that poles of y and of wp are regular
  x_from_y                      cubic-quintic solution x(t) from (eq:alg1sol), regular through poles of y
  dop853_reference              SciPy DOP853 reference with escape detection
  taylor_reference              30-digit Taylor-series reference (Python decimal; no mpmath needed)
  metrics                       e_x and hat e_x of Eq. (metrics) of the paper
"""
from __future__ import annotations

import re
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np
from scipy.integrate import solve_ivp

from _appendix_tex import APPB_PHIPSI, APPB_SEXTIC, APPB_G, APPB_W

# ----------------------------------------------------------------------------------------------------------
# 1. Printed polynomials
# ----------------------------------------------------------------------------------------------------------
_SYM = [(r"C_0", "C"), (r"x_0", "x"), (r"\alpha", "a"), (r"\beta", "b"), (r"\gamma", "c"),
        (r"\lambda", "l"), (r"\mu", "m"), ("a", "a"), ("b", "b"), ("c", "c")]


def _clean(s: str) -> str:
    s = re.sub(r"\\label\{[^}]*\}", "", s)
    s = re.sub(r"\\(begin|end)\{align\}", "", s)
    for t in (r"\nonumber", r"\\", "&", r"\quad", "{}", r"\,"):
        s = s.replace(t, " " if t != r"\," else "")
    return s


def parse_latex_poly(body: str):
    """Parse a printed polynomial such as '- 1\\,296C_0x_0^{2}\\alpha^{5} \\beta^{4} + ...'.
    Returns a list of (Fraction coefficient, {var: power})."""
    s = re.sub(r"\s+", "", _clean(body))
    terms = []
    for m in re.finditer(r"([+-]?)([^+-]+)", s):
        sign, mono = m.group(1), m.group(2)
        num = re.match(r"\d+", mono)
        coef = Fraction(int(num.group(0))) if num else Fraction(1)
        rest = mono[num.end():] if num else mono
        powers = {}
        while rest:
            for tex, v in _SYM:
                if rest.startswith(tex):
                    rest = rest[len(tex):]
                    p = 1
                    pm = re.match(r"\^\{(\d+)\}|\^(\d)", rest)
                    if pm:
                        p = int(pm.group(1) or pm.group(2))
                        rest = rest[pm.end():]
                    powers[v] = powers.get(v, 0) + p
                    break
            else:
                raise ValueError(f"cannot parse monomial {mono!r} at {rest!r}")
        terms.append((-coef if sign == "-" else coef, powers))
    return terms


def _split_defs(tex: str):
    """Split an align block 'NAME={}& ... NAME2={}& ...' into {NAME: body}."""
    s = re.sub(r"\\label\{[^}]*\}", "", tex)
    s = re.sub(r"\\(begin|end)\{align\}", "", s)
    parts = re.split(r"(\\Phi\(\\lambda,\\mu\)|\\Psi\(\\lambda,\\mu\)|e_\d|W)=\{\}&", s)
    out = {}
    for i in range(1, len(parts), 2):
        key = parts[i].replace("\\", "").replace("(lambda,mu)", "")
        out[key] = parts[i + 1]
    return out


_pp = _split_defs(APPB_PHIPSI)
_ss = _split_defs(APPB_SEXTIC)
_ww = _split_defs(APPB_W)
POLY = {
    "Phi": parse_latex_poly(_pp["Phi"]),
    "Psi": parse_latex_poly(_pp["Psi"]),
    "G": parse_latex_poly(APPB_G),
    "W": parse_latex_poly(_ww["W"]),
}
for k in range(7):
    POLY[f"e{k}"] = parse_latex_poly(_ss[f"e_{k}"])
KDELTA = 3442605166011971360784384          # leading constant of (eq:disc) as printed


def ev(terms, env):
    """Float evaluation; env maps a,b,c,C,x,l,m to numbers."""
    tot = 0.0
    for c, d in terms:
        t = float(c)
        for v, p in d.items():
            t *= env[v] ** p
        tot += t
    return tot


def ev_abs(terms, env):
    """Sum of absolute values of the terms (scale for relative residuals)."""
    tot = 0.0
    for c, d in terms:
        t = abs(float(c))
        for v, p in d.items():
            t *= abs(env[v]) ** p
        tot += t
    return tot


def ev_exact(terms, env):
    tot = Fraction(0)
    for c, d in terms:
        t = c
        for v, p in d.items():
            t *= env[v] ** p
        tot += t
    return tot


def coeffs_in(terms, var, env):
    """Coefficients (low -> high) of a polynomial in `var`, other variables evaluated (float)."""
    deg = max(d.get(var, 0) for _, d in terms)
    co = [0.0] * (deg + 1)
    for c, d in terms:
        t = float(c)
        for v, p in d.items():
            if v != var:
                t *= env[v] ** p
        co[d.get(var, 0)] += t
    return co


def C0_of(a, b, c, x0, v0):
    """Energy C0 with exactly rounded summation (math.fsum): near a minimum the terms cancel, and naive
    summation loses v0^2/2 entirely (e.g. a=1, b=-4, c=3, x0=1, v0=1e-8 gives C0 = 0 instead of 5e-17)."""
    import math
    return math.fsum([v0 * v0 / 2, a * x0 ** 2 / 2, b * x0 ** 4 / 4, c * x0 ** 6 / 6])


# ----------------------------------------------------------------------------------------------------------
# 2. Weierstrass wp and the universal formula (eq:wp)
# ----------------------------------------------------------------------------------------------------------
def wp_roots(p, q, hy):
    """Roots of 4z^3 - g2 z - g3 for the invariants (eq:g2g3) of the even quartic h_y - p y^2 - q y^4/2.
    The resolvent factors exactly:  4z^3 - g2 z - g3 = (z - p/6)(4z^2 + (2p/3) z + p^2/36 + q h_y/2),
    so e = p/6 and e = -p/12 +- sqrt(-q h_y/2)/2. Using these closed forms instead of numerical roots
    of the cubic keeps full precision when h_y is tiny (orbits near a separatrix of y)."""
    rad = -0.5 * q * hy
    return p / 6.0, -p / 12.0, 0.5 * np.sqrt(abs(rad)), rad >= 0


def jacobi_sncndn(u, m, m1):
    """sn, cn, dn (u | m) by the arithmetic-geometric mean (A&S 16.4), with the complementary parameter
    m1 = 1 - m supplied separately so that orbits within 1e-18 of a separatrix (m1 tiny) keep full
    precision. Vectorised over u; 0 <= m <= 1."""
    u = np.asarray(u, dtype=float)
    if m1 <= 0.0:                                               # m = 1
        th = np.tanh(u)
        se = 1.0 / np.cosh(u)
        return th, se, se
    if m <= 0.0:
        return np.sin(u), np.cos(u), np.ones_like(u)
    a, b, c = [1.0], [np.sqrt(m1)], [np.sqrt(m)]
    while abs(c[-1]) > 1e-17 * a[-1] and len(a) < 60:
        a.append(0.5 * (a[-1] + b[-1]))
        b.append(np.sqrt(a[-2] * b[-1]))
        c.append(0.5 * (a[-2] - b[-2]))
    N = len(a) - 1
    phi = (2.0 ** N) * a[N] * u
    phis = [phi]
    for n in range(N, 0, -1):
        phi = 0.5 * (phi + np.arcsin(np.clip(c[n] / a[n] * np.sin(phi), -1.0, 1.0)))
        phis.append(phi)
    sn, cn = np.sin(phi), np.cos(phi)
    dn = np.sqrt(np.maximum(m1 + m * cn * cn, 0.0))
    return sn, cn, dn


def _wp_parts(t, p, q, hy):
    """Return (u, w) with u = 1/wp(t), w = wp'(t)/wp(t)^2, both regular at the poles of wp.
    Classical Jacobi representations of wp (DLMF 23.6.16-23.6.17; A&S 18.9.11); the parameter m and its
    complement 1 - m are both formed from differences of the exact roots, without cancellation."""
    t = np.asarray(t, dtype=float)
    eA, c0, r, real = wp_roots(p, q, hy)
    if real:                                                    # three real e_i
        e1, e2, e3 = sorted([eA, c0 + r, c0 - r], reverse=True)
        k2 = e1 - e3
        if k2 <= 0:                                             # triple root (p = 0, q h_y = 0): wp = 1/t^2
            return t ** 2, -2 * t
        m, m1 = (e2 - e3) / k2, (e1 - e2) / k2
        sn, cn, dn = jacobi_sncndn(np.sqrt(k2) * t, m, m1)
        den = e3 * sn ** 2 + k2                                 # wp = den / sn^2
        return sn ** 2 / den, -2 * k2 ** 1.5 * cn * dn * sn / den ** 2
    e2 = eA                                                     # one real root, complex pair c0 +- i r
    d = e2 - c0                                                 # = p/4
    H2 = np.sqrt(d * d + r * r)                                 # |e2 - e1| = sqrt((e2-e1)(e2-e3))
    # m = 1/2 - 3 e2/(4 H2) = (H2 - d)/(2 H2),  1 - m = (H2 + d)/(2 H2); the smaller one via r^2
    if d >= 0:
        m, m1 = r * r / ((H2 + d) * 2 * H2), (H2 + d) / (2 * H2)
    else:
        m, m1 = (H2 - d) / (2 * H2), r * r / ((H2 - d) * 2 * H2)
    sn, cn, dn = jacobi_sncndn(2 * np.sqrt(H2) * t, m, m1)
    den = e2 * (1 - cn) + H2 * (1 + cn)                         # wp = den / (1 - cn)
    return (1 - cn) / den, -4 * H2 ** 1.5 * sn * dn / den ** 2


def wp_invariants(p, q, y0, y1):
    hy = y1 * y1 + p * y0 * y0 + 0.5 * q * y0 ** 4
    g2 = p * p / 12 - q * hy / 2
    g3 = p ** 3 / 216 + p * q * hy / 12
    return hy, g2, g3


def wp_y(t, p, q, y0, y1, hy=None):
    """Solution of y'' + p y + q y^3 = 0, y(0)=y0, y'(0)=y1 by (eq:wp), as y = Ny/Dy.
    Numerator and denominator of (eq:wp) are divided by wp^2, which makes them regular where wp has a pole;
    Dy = 0 marks a pole of y."""
    if hy is None:
        hy, _, _ = wp_invariants(p, q, y0, y1)
    u, w = _wp_parts(t, p, q, hy)
    k = (p + 3 * q * y0 * y0) / 12
    num = y1 * w + y0 * (p + q * y0 * y0) * (u + k * u * u) + 0.5 * q * y0 * y1 * y1 * u * u
    den = 2 * (1 + k * u) ** 2 + 0.25 * q * y1 * y1 * u * u
    return y0 * den - num, den


def wp_y_value(t, p, q, y0, y1):
    Ny, Dy = wp_y(t, p, q, y0, y1)
    with np.errstate(divide="ignore", invalid="ignore"):
        return Ny / Dy


def x_from_y(Ny, Dy, S, lam, mu, signed=True, sgn_x0=1.0):
    """x = sqrt(S) y / sqrt(1 + lam y^2 + mu y^4), y = Ny/Dy.
    Where y is moderate (|Dy| >= 1e-3 |Ny|) the formula is used as printed; near a pole of y it is written as
    sqrt(S) Ny Dy / sqrt(Dy^4 + lam Ny^2 Dy^2 + mu Ny^4), which is regular there when mu > 0 (Proposition 1a).
    The split also avoids the removable 0/0 points of (eq:wp) (Ny = Dy = 0 with y finite).
    With signed=False, y is replaced by sgn(x0)|y| (family I, Algorithm 1 step 7).
    Returns (x, Dfac) with Dfac > 0 inside the existence interval (sign of 1 + lam y^2 + mu y^4)."""
    Ny = np.asarray(Ny, dtype=float)
    Dy = np.asarray(Dy, dtype=float)
    near_pole = np.abs(Dy) < 1e-3 * np.abs(Ny)
    with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
        y = np.where(near_pole, 0.0, Ny / np.where(near_pole, 1.0, Dy))
        ys = y if signed else sgn_x0 * np.abs(y)
        D1 = 1 + lam * y * y + mu * y ** 4
        x1 = np.sqrt(S) * ys / np.sqrt(D1)
        Q = Dy ** 4 + lam * Ny ** 2 * Dy ** 2 + mu * Ny ** 4
        num = Ny * Dy if signed else sgn_x0 * np.abs(Ny * Dy)
        x2 = np.sqrt(S) * num / np.sqrt(Q)
    return np.where(near_pole, x2, x1), np.where(near_pole, Q, D1)


# ----------------------------------------------------------------------------------------------------------
# 3. References and metrics
# ----------------------------------------------------------------------------------------------------------
def rhs(a, b, c):
    return lambda t, z: [z[1], -(a * z[0] + b * z[0] ** 3 + c * z[0] ** 5)]


def dop853_reference(a, b, c, x0, v0, t_max=20.0, x_esc=1e3, rtol=2.5e-14, atol=1e-16, npts=401):
    """Grid t_j = j t_w/400, t_w = min(t_max, 0.9 t_esc); t_esc: first time |x| reaches x_esc."""
    f = rhs(a, b, c)

    def esc(t, z):
        return abs(z[0]) - x_esc
    esc.terminal = True
    pre = solve_ivp(f, (0, t_max), [x0, v0], method="DOP853", rtol=1e-12, atol=1e-14, events=esc)
    t_esc = float(pre.t_events[0][0]) if len(pre.t_events[0]) else np.inf
    tw = min(t_max, 0.9 * t_esc)
    ts = np.linspace(0, tw, npts)
    ref = solve_ivp(f, (0, tw), [x0, v0], method="DOP853", t_eval=ts, rtol=rtol, atol=atol)
    return ts, ref.y[0], t_esc


def taylor_reference(a, b, c, x0, v0, ts, digits=30, order=30, h=None, safety=0.25, t_stop_abs=1e6):
    """Taylor-series solution of x'' = -(a x + b x^3 + c x^5) in `digits`-digit decimal arithmetic, with
    the step chosen from the radius of convergence estimated from the last two series coefficients
    (step = safety * radius, truncation ~ safety^order). Used where DOP853 is not converged: orbits at a
    separatrix level and stiff large-amplitude orbits. Values after |x| exceeds t_stop_abs are nan."""
    getcontext().prec = digits + 8
    # Exact binary values of the double inputs: decimal strings repr(float(.)) define a
    # slightly different problem, which near a separatrix moves the solution by up to 1e-6
    A, B, C = Decimal(float(a)), Decimal(float(b)), Decimal(float(c))
    x, v = Decimal(float(x0)), Decimal(float(v0))

    def series(x, v):
        X = [x, v] + [Decimal(0)] * (order - 1)
        X2 = [Decimal(0)] * (order + 1)
        X3 = [Decimal(0)] * (order + 1)
        X5 = [Decimal(0)] * (order + 1)
        for n in range(order - 1):
            X2[n] = sum(X[i] * X[n - i] for i in range(n + 1))
            X3[n] = sum(X2[i] * X[n - i] for i in range(n + 1))
            X5[n] = sum(X3[i] * X2[n - i] for i in range(n + 1))
            X[n + 2] = -(A * X[n] + B * X3[n] + C * X5[n]) / ((n + 1) * (n + 2))
        return X

    def horner(X, dt):
        acc = Decimal(0)
        for coef in reversed(X):
            acc = acc * dt + coef
        return acc

    out = np.full(len(ts), np.nan)
    t = Decimal(0)
    ts_d = [Decimal(float(s)) for s in ts]
    j = 0
    while j < len(ts_d):
        X = series(x, v)
        if h is None:
            rad = []
            for n in (order - 1, order):
                if X[n] != 0:
                    rad.append(float(abs(X[n])) ** (-1.0 / n))
            H = Decimal(repr(safety * min(rad))) if rad else Decimal(1)
            H = min(H, Decimal(1))
        else:
            H = Decimal(repr(h))
        while j < len(ts_d) and ts_d[j] <= t + H:
            out[j] = float(horner(X, ts_d[j] - t))
            j += 1
        xn = Decimal(0)
        vn = Decimal(0)
        for n in range(len(X) - 1, -1, -1):
            xn = xn * H + X[n]
            if n >= 1:
                vn = vn * H + n * X[n]
        x, v = xn, vn
        t += H
        if abs(x) > t_stop_abs:
            break
    return out


def metrics(x, xref):
    e = float(np.max(np.abs(x - xref)))
    return e, e / max(float(np.max(np.abs(xref))), 1e-300)
