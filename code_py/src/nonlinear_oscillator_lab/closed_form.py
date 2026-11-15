"""Closed-form solutions (A. H. Salas) -- Python port of the Wolfram Language solver.

This module is a line-by-line port of ``duffingSolution``, ``cqCandidates``,
``cqBuild`` and ``cqSolution`` in ``mathematica/NonlinearOscillatorLab.wl``.

* :func:`duffing_solution` solves ``y'' + p y + q y**3 = 0`` with ``y(0)=y0``,
  ``y'(0)=y1`` in Jacobi-elliptic (or hyperbolic) form.
* :func:`exact_solution` solves ``x'' + a x + b x**3 + c x**5 = 0`` through the
  representation ``x = A y / sqrt(1 + lam*y**2 + mu*y**4)``, where ``y`` solves a
  cubic Duffing equation.

Jacobi functions are evaluated with mpmath for arbitrary real parameter ``m``
(including ``m < 0`` and ``m > 1``, which occur in this construction), so the
port reproduces the branch structure of the Wolfram implementation.  Every
candidate representation is accepted only after the same consistency checks as
in the Wolfram code: positivity of the denominator, ``x(0) = x0``, agreement with
the local Taylor expansion, and the first integral ``x'**2 = 2 (E - V(x))``.
If no candidate passes, :class:`ClosedFormNotAvailable` is raised.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np
import mpmath as mp

__all__ = ["DuffingSolution", "ClosedFormSolution", "ClosedFormNotAvailable",
           "duffing_solution", "exact_solution"]

_DPS = 40                       # digits of the construction stage (roots, parameters)
_TOL = 1e-9                     # tolerance of cqCandidates (as in the WL code)


class ClosedFormNotAvailable(ValueError):
    """No real closed-form representation passed the consistency checks."""


def _ctx():
    ctx = mp.mp.clone()
    ctx.dps = _DPS
    return ctx


def _real(z, rel=1e-8):
    """float(z) if z is (numerically) real, else nan."""
    z = complex(z)
    if abs(z.imag) <= rel * (1.0 + abs(z.real)):
        return z.real
    return float("nan")


def _clip(v, lo, hi):
    return min(max(v, lo), hi)


# --- inverse Jacobi functions, principal branches as in the Wolfram Language ---
def _inv_cn(ctx, v, m):
    return ctx.ellipf(ctx.acos(v), m)


def _inv_sn(ctx, v, m):
    return ctx.ellipf(ctx.asin(v), m)


def _inv_dn(ctx, v, m):
    return ctx.ellipf(ctx.asin(ctx.sqrt((1 - v * v) / m)), m)


def _ellipj_real(u, m):
    """sn, cn, dn for real u (array) and any real parameter m (A&S 16.10-16.11)."""
    from scipy.special import ellipj
    u = np.asarray(u, dtype=float)
    if 0.0 <= m <= 1.0:
        return ellipj(u, m)[:3]
    if m < 0.0:                                   # A&S 16.10, negative parameter
        mu = -m / (1.0 - m); r = np.sqrt(1.0 - m)
        s, c, d, _ = ellipj(u * r, mu)
        return s / (d * r), c / d, 1.0 / d
    k = np.sqrt(m)                                # A&S 16.11, parameter > 1
    s, c, d, _ = ellipj(u * k, 1.0 / m)
    return s / k, d, c


class _Jac:
    """Evaluate sn, cn, dn at u = c1 + sw*t, vectorised over t.

    Uses SciPy (with the real parameter transformations) for real arguments.  Falls
    back to mpmath, with the unrounded phase c1 and frequency sw, when c1 or sw is
    complex (principal-branch inverse functions outside their real range) or when the
    effective parameter passed to SciPy lies within 1e-8 of 1, where ``ellipj`` loses
    accuracy (this occurs for orbits extremely close to a separatrix)."""
    def __init__(self, ctx, c1, sw, m):
        self.ctx, self.m = ctx, float(m)
        self.c1_mp, self.sw_mp = c1, sw
        c1, sw = complex(c1), complex(sw)
        real = abs(c1.imag) < 1e-14 * (1 + abs(c1.real)) and abs(sw.imag) < 1e-14 * (1 + abs(sw.real))
        mf = self.m
        m_eff = mf if 0.0 <= mf <= 1.0 else (-mf / (1.0 - mf) if mf < 0.0 else 1.0 / mf)
        self.real = real and not (0.0 < 1.0 - m_eff < 1e-8)
        self.c1, self.sw = (c1.real, sw.real) if self.real else (c1, sw)
        self.m_mp = m

    def __call__(self, t):
        t = np.asarray(t, dtype=float)
        if self.real:
            return _ellipj_real(self.c1 + self.sw * t, self.m)
        ctx = self.ctx
        with ctx.workdps(25):
            out = [[complex(ctx.ellipfun(k, self.c1_mp + self.sw_mp * ctx.mpf(float(s)), m=self.m_mp))
                    for s in t.ravel()] for k in ("sn", "cn", "dn")]
        return tuple(np.array(o).reshape(t.shape) for o in out)


def _realv(z, rel=1e-8):
    """Real part of z where numerically real, nan elsewhere (vectorised)."""
    z = np.asarray(z)
    if not np.iscomplexobj(z):
        return z.astype(float)
    ok = np.abs(z.imag) <= rel * (1.0 + np.abs(z.real))
    return np.where(ok, z.real, np.nan)


@dataclass
class DuffingSolution:
    """Closed-form solution of y'' + p y + q y^3 = 0."""
    case: str
    y: Callable[[float], float]
    m: Optional[float] = None
    omega: Optional[float] = None
    c1: Optional[complex] = None
    period: Optional[float] = None
    blow_up_time: Optional[float] = None
    amplitude: Optional[float] = None

    @property
    def bounded(self) -> bool:
        return self.blow_up_time is None


def _sgn(v):
    return 1 if v > 0 else (-1 if v < 0 else 0)


def duffing_solution(p, q, y0, y1) -> DuffingSolution:
    """Closed-form solution of ``y'' + p y + q y**3 = 0``, ``y(0)=y0``, ``y'(0)=y1``.

    Port of ``duffingSolution`` (Wolfram Language).  The parameters (roots, modulus,
    frequency, phase) are computed with ``_DPS``-digit mpmath arithmetic, because near a
    separatrix they follow from nearly coalescing roots; the solution itself is evaluated
    in double precision.  Raises :class:`ClosedFormNotAvailable` in the degenerate case
    not covered by the Wolfram solver.
    """
    ctx = _ctx()
    p, q, y0, y1 = (ctx.mpf(v) for v in (p, q, y0, y1))
    fy0, fy1 = float(y0), float(y1)
    h = y1 * y1 + p * y0 * y0 + q * y0 ** 4 / 2
    X = y0 * y0
    if q == 0:
        if p == 0:
            return DuffingSolution("Linear", lambda t: fy0 + fy1 * np.asarray(t, float))
        sp = np.sqrt(complex(float(p)))
        return DuffingSolution(
            "Linear", lambda t: np.real(fy0 * np.cos(sp * np.asarray(t, float))
                                        + fy1 * np.sin(sp * np.asarray(t, float)) / sp))
    dd = p * p + 2 * q * h
    roots = [] if dd < 0 else sorted([(-p + ctx.sqrt(dd)) / q, (-p - ctx.sqrt(dd)) / q])
    if y1 == 0 and roots:
        k = min(range(len(roots)), key=lambda i: abs(roots[i] - X))
        roots[k] = X
    pos = [r for r in roots if r > 0]

    # orbit around 0: A cn
    if h > 0 and pos and X <= min(pos):
        A2 = min(pos); A = ctx.sqrt(A2); w = p + q * A2; m = q * A2 / (2 * w)
        c1 = (-1 if y1 >= 0 else 1) * _inv_cn(ctx, _clip(y0 / A, -1, 1), m)
        sw = ctx.sqrt(w); J = _Jac(ctx, c1, sw, m); fA = float(A)
        return DuffingSolution(
            "cn: orbit around 0", lambda t: _realv(fA * J(t)[1]),
            m=float(m), omega=float(w), c1=complex(c1), period=_real(4 * ctx.ellipk(m) / sw), amplitude=fA)
    if q > 0:
        if h == 0:                                          # separatrix
            if y0 == 0:
                return DuffingSolution("Equilibrium", lambda t: np.zeros_like(np.asarray(t, float)))
            a = ctx.sqrt(-2 * p / q)
            c1 = (-1 if y0 * y1 > 0 else 1) * ctx.acosh(a / abs(y0))
            sq = float(ctx.sqrt(-p)); fa = _sgn(y0) * float(a); fc1 = float(c1)
            return DuffingSolution(
                "Separatrix (sech)",
                lambda t: fa / np.cosh(sq * np.asarray(t, float) + fc1), c1=complex(c1))
        # one well: B dn
        B2 = max(roots); b2 = min(roots); B = ctx.sqrt(B2)
        m = 1 - b2 / B2; w = q * B2 / 2
        c1 = (-1 if y0 * y1 > 0 else 1) * _inv_dn(ctx, _clip(abs(y0) / B, ctx.sqrt(max(1 - m, 0)), 1), m)
        sw = ctx.sqrt(w); fsB = _sgn(y0) * float(B); J = _Jac(ctx, c1, sw, m)
        return DuffingSolution(
            "dn: oscillation inside one well",
            lambda t: _realv(fsB * J(t)[2]),
            m=float(m), omega=float(w), c1=complex(c1), period=_real(2 * ctx.ellipk(m) / sw), amplitude=float(B))
    # q < 0
    k = -q
    if h > 0 and not pos:
        R = ctx.sqrt(2 * h / k); m = ctx.mpf(1) / 2 + p / (2 * k * R); w = 2 * k * R
        sg = 1 if y1 >= 0 else -1
        c1 = _sgn(y0) * _inv_cn(ctx, (R - X) / (R + X), m)
        sw = ctx.sqrt(w); sR = float(ctx.sqrt(R)); J = _Jac(ctx, c1, sg * sw, m)

        def y(t):
            sn, cn, _ = J(t)
            with np.errstate(all="ignore"):
                return _realv(sR * sn / (1 + cn))
        return DuffingSolution("Crosses 0 and blows up", y, m=float(m), omega=float(w), c1=complex(c1),
                               blow_up_time=_real((2 * ctx.ellipk(m) - sg * c1) / sw))
    sel = [r for r in pos if r <= X]
    if not sel:
        raise ClosedFormNotAvailable("Degenerate case")
    ua = max(sel)
    below = sorted(r for r in roots + [ctx.zero] if r < ua)
    if len(below) < 2:
        raise ClosedFormNotAvailable("Degenerate case")
    r, s = below[-1], below[-2]
    m = (r - s) / (ua - s); w = k / 2 * (ua - s)
    sn0 = ctx.sqrt(_clip((X - ua) / (X - r), 0, 1))
    c1 = (1 if y0 * y1 >= 0 else -1) * _inv_sn(ctx, sn0, m)
    sw = ctx.sqrt(w); sg = _sgn(y0); J = _Jac(ctx, c1, sw, m); fua, fr = float(ua), float(r)

    def y(t):
        sn, cn, _ = J(t)
        with np.errstate(all="ignore"):
            return _realv(sg * np.sqrt((fua - fr * sn ** 2).astype(complex)) / cn)
    return DuffingSolution("Outside barrier: blows up", y, m=float(m), omega=float(w), c1=complex(c1),
                           blow_up_time=_real((ctx.ellipk(m) - c1) / sw))


# ------------------------------------------------------------------ cubic-quintic

def _wl_order(ctx, roots):
    """Order polynomial roots as the Wolfram implementation enumerates them:
    real roots in decreasing order, then complex roots.  With this order the
    port selects the same representation as the Wolfram solver in all 200
    cases of the published randomised campaign."""
    roots = [ctx.mpc(r) for r in roots]
    real = sorted((r for r in roots if abs(r.imag) < _TOL), key=lambda z: -z.real)
    cplx = sorted((r for r in roots if abs(r.imag) >= _TOL), key=lambda z: (z.real, z.imag))
    return real + cplx


def _polyroots(ctx, coeffs):
    """All complex roots of a polynomial (highest degree first) at ctx precision."""
    coeffs = list(coeffs)
    while coeffs and coeffs[0] == 0:
        coeffs.pop(0)
    if len(coeffs) < 2:
        return []
    return list(ctx.polyroots(coeffs, maxsteps=200, extraprec=4 * ctx.prec))


def _candidates(al, be, ga, x0, v0):
    """Candidate representations (port of cqCandidates), computed at _DPS digits."""
    ctx = _ctx()
    al, be, ga, x0, v0 = (ctx.mpf(v) for v in (al, be, ga, x0, v0))
    h = v0 * v0 + al * x0 ** 2 + be * x0 ** 4 / 2 + ga * x0 ** 6 / 3
    X = x0 * x0
    roots = _wl_order(ctx, _polyroots(ctx, [h, -al, -be / 2, -ga / 3]))
    if ga == 0:
        roots = _wl_order(ctx, _polyroots(ctx, [h, -al, -be / 2])) + [ctx.mpc(0)]
    out = []
    for j, rj in enumerate(roots):
        if abs(rj.imag) >= _TOL:
            continue
        pj = rj.real
        if abs(1 - pj * X) > _TOL:                                    # family I (mu = 0)
            s = 1 / (1 - pj * X)
            out.append(dict(family="I", s=s, lam=s * pj, mu=ctx.zero, lm_disc=(s * pj) ** 2,
                            p=al - 3 * h * pj, q=s * (be + 4 * al * pj - 6 * h * pj * pj), y1=s * v0))
        kl = [r for i, r in enumerate(roots) if i != j]               # family II (mu != 0)
        sig = ctx.fsum(kl) if kl else ctx.mpc(0)
        pp = ctx.fprod(kl) if kl else ctx.mpc(1)
        if abs(sig.imag) < _TOL and abs(pp.imag) < _TOL and 1 - sig.real * X + pp.real * X * X >= 0:
            # sigma^2 - 4 pi = (k1 - k2)^2, evaluated without cancellation
            disc = (kl[0] - kl[1]) ** 2 if len(kl) == 2 else sig * sig - 4 * pp
            sig, pp, disc = sig.real, pp.real, ctx.re(disc)
            s = 4 / (2 - sig * X + 2 * ctx.sqrt(1 - sig * X + pp * X * X))
            lam = s * sig / 2; mu = s * s * disc / 16
            if abs(1 - mu * X * X) > _TOL:
                # lam^2 - 4 mu = s^2 * pi exactly (Vieta), used to factor the denominator
                out.append(dict(family="II", s=s, lam=lam, mu=mu, lm_disc=s * s * pp, p=h * (pj - sig / 2),
                                q=-h * s * disc / 8, y1=s * v0 / (1 - mu * X * X)))
    return out


@dataclass
class ClosedFormSolution:
    """Closed-form solution x(t) = A y(t)/sqrt(1 + lam y^2 + mu y^4)."""
    representation: str
    A: float
    lam: float
    mu: float
    p: float
    q: float
    y0: float
    y1: float
    duffing: DuffingSolution = field(repr=False)
    _x: Callable[[float], float] = field(repr=False)
    time_shift: float = 0.0

    @property
    def bounded_y(self) -> bool:
        return self.duffing.bounded

    @property
    def family(self) -> str:
        return self.representation.split(" ")[0]

    def x(self, t):
        """Evaluate x(t); accepts a scalar or an array."""
        out = self._x(np.asarray(t, dtype=float) + self.time_shift)
        return float(out) if np.ndim(t) == 0 else out

    def __call__(self, t):
        return self.x(t)


def _build(al, be, ga, x0, v0, c):
    try:
        ys = duffing_solution(c["p"], c["q"], x0, c["y1"])
    except (ClosedFormNotAvailable, ZeroDivisionError, ValueError, FloatingPointError):
        return None
    yf = ys.y
    if not c["s"] > 0:
        return None
    lam, mu, sqs = float(c["lam"]), float(c["mu"]), float(c["s"] ** 0.5)
    signed = ys.case.startswith("cn") or ys.case.startswith("Linear") or (lam == 0 and mu == 0)
    sx = 1.0 if signed else np.sign(x0)

    # D(y) = 1 + lam y^2 + mu y^4 is evaluated in factored form (1 + al1 y^2)(1 + al2 y^2),
    # al1 + al2 = lam, al1 al2 = mu (factors from 40-digit arithmetic).  Near separatrices
    # |lam|, |mu| become very large and D nearly a perfect square, so the expanded form
    # would lose up to all significant digits by cancellation.
    lam_mp, disc = c["lam"], c["lm_disc"]
    if disc >= 0:
        r_ = disc ** 0.5
        big = (lam_mp + r_) / 2 if lam_mp >= 0 else (lam_mp - r_) / 2
        small = c["mu"] / big if big != 0 else lam_mp - big
        fa1, fa2 = float(big), float(small)

        def den(y):
            y2 = y * y
            return (1 + fa1 * y2) * (1 + fa2 * y2)
    else:
        fre, fim = float(lam_mp / 2), float((-disc) ** 0.5 / 2)

        def den(y):
            y2 = y * y
            return (1 + fre * y2) ** 2 + (fim * y2) ** 2

    def xf(t):
        with np.errstate(all="ignore"):
            y = yf(t)
            d = den(y)
            num = y if signed else np.abs(y)
            return np.where(d > 0, sx * sqs * num / np.sqrt(np.where(d > 0, d, 1.0)), np.nan)

    tw = 6.0
    if lam == 0 and mu == 0 and ys.blow_up_time is not None and np.isfinite(ys.blow_up_time):
        tw = min(tw, 0.9 * ys.blow_up_time)
    h = v0 * v0 + al * x0 ** 2 + be * x0 ** 4 / 2 + ga * x0 ** 6 / 3

    def Pf(u):
        return h - al * u - be * u * u / 2 - ga * u ** 3 / 3
    ts = np.arange(1, 401) * tw / 400
    with np.errstate(all="ignore"):
        yv = np.asarray(yf(ts), dtype=float)
        good = np.isfinite(yv) & (np.abs(yv) < 1e6)
        if not good.any() or np.min(den(yv[good])) <= 0:
            return None
        x00 = float(xf(0.0))
        if not np.isfinite(x00) or abs(x00 - x0) > 1e-8:
            return None
        dl = 1e-3
        acc = -al * x0 - be * x0 ** 3 - ga * x0 ** 5
        tay = x0 + v0 * dl + acc * dl * dl / 2
        if abs(tay - x0) > 0 and not (abs(float(xf(dl)) - tay) <= 0.05 * abs(tay - x0)):
            return None
        dt = 1e-6
        xv = xf(ts); xp = (xf(ts + dt) - xf(ts - dt)) / (2 * dt)
        ok = np.isfinite(xp) & np.isfinite(xv) & good
        if np.any(good & ~np.isfinite(xv)):
            return None
        if ok.any():
            P = Pf(xv[ok] ** 2)
            if np.max(np.abs(xp[ok] ** 2 - P) / (1 + np.abs(P))) > 1e-3:
                return None
    return ClosedFormSolution(f"{c['family']} / y: {ys.case}", sqs, lam, mu, float(c["p"]), float(c["q"]),
                              x0, float(c["y1"]), ys, xf)


def exact_solution(model_or_a, *args) -> ClosedFormSolution:
    """Closed-form solution of x'' + a x + b x^3 + c x^5 = 0 (A. H. Salas).

    Call as ``exact_solution(model, x0, v0)`` or ``exact_solution(a, b, c, x0, v0)``.
    Returns a :class:`ClosedFormSolution` whose ``x(t)`` evaluates the representation,
    or raises :class:`ClosedFormNotAvailable`.
    """
    if hasattr(model_or_a, "a"):
        al, be, ga = model_or_a.a, model_or_a.b, model_or_a.c
        x0, v0 = args
    else:
        al = model_or_a
        be, ga, x0, v0 = args
    al, be, ga, x0, v0 = (float(v) for v in (al, be, ga, x0, v0))
    if v0 == 0 and al * x0 + be * x0 ** 3 + ga * x0 ** 5 == 0:          # exact equilibrium
        eq = DuffingSolution("Equilibrium", lambda t: np.full(np.shape(t), x0, dtype=float))
        return ClosedFormSolution("Equilibrium", abs(x0), 0.0, 0.0, 0.0, 0.0, x0, 0.0, eq,
                                  lambda t: np.full(np.shape(t), x0, dtype=float))
    try:
        cands = _candidates(al, be, ga, x0, v0)
    except (ZeroDivisionError, ValueError) as exc:
        raise ClosedFormNotAvailable("Degenerate parameters: %s" % type(exc).__name__)
    good = [s for s in (_build(al, be, ga, x0, v0, c) for c in cands) if s]
    if good:
        good.sort(key=lambda s: 0 if s.bounded_y else 1)            # stable, as SortBy
        return good[0]
    # start at rest at a turning point that is a pole of every y: start from the
    # other end of the motion and shift time
    if v0 == 0:
        h = al * x0 ** 2 + be * x0 ** 4 / 2 + ga * x0 ** 6 / 3

        def Pf(w):
            return h - al * w - be * w * w / 2 - ga * w ** 3 / 3
        coeffs = [ga / 3, be / 2, al, -h]
        while len(coeffs) > 1 and coeffs[0] == 0:
            coeffs.pop(0)
        roots = [float(np.real(r)) for r in np.roots(coeffs) if abs(np.imag(r)) < 1e-10]
        others = sorted([u for u in roots if u > 0 and abs(u - x0 * x0) > 1e-9 and Pf((u + x0 * x0) / 2) > 0],
                        key=lambda u: abs(u - x0 * x0))
        for uo in others:
            xo = np.sign(x0) * np.sqrt(uo)
            so = [s for s in (_build(al, be, ga, xo, 0.0, c) for c in _candidates(al, be, ga, xo, 0.0))
                  if s and not s.bounded_y]
            if so:
                sol = so[0]
                sol.representation += " (time-shifted)"
                sol.time_shift = float(sol.duffing.blow_up_time)
                return sol
    raise ClosedFormNotAvailable("No candidate representation passed the acceptance checks")
