"""Energy methods for x'' + a*x + b*x**3 + c*x**5 = 0.

No assertion of analytic closed-form trajectories is made. Period is computed
by numerical quadrature of the factored, singularity-free integrand
obtained from E - V(x) = (x - x_-)(x_+ - x) Q(x) and x = mid + half*sin(theta).
"""
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from math import isfinite, sqrt
import numpy as np
from scipy.integrate import solve_ivp, quad

@dataclass(frozen=True)
class Model:
    a: float = 1.0
    b: float = 0.0
    c: float = 0.0

    def potential(self, x):
        x = np.asarray(x)
        # Fraction coordinates arrive in object arrays.
        if x.dtype.kind in "biuO":
            x = x.astype(float)
        c = float(self.c) if isinstance(self.c, Fraction) else self.c
        with np.errstate(over="ignore", invalid="ignore"):
            value = 0.5*self.a*x*x + 0.25*self.b*x**4 + c*x**6/6
        if value.ndim == 0 and np.isfinite(value):
            return value
        retry = ~np.isfinite(value) & np.isfinite(x)
        if not np.any(retry) or not all(isfinite(k) for k in (self.a, self.b, self.c)):
            return value

        # Recover finite potentials when an unscaled power overflows before
        # multiplication by its coefficient. Round the binary-input polynomial once.
        a, b, c = map(lambda k: Fraction(float(k)), (self.a, self.b, self.c))
        value = np.array(value, dtype=np.result_type(value.dtype, np.float64), copy=True)
        for index in np.flatnonzero(retry):
            y = Fraction(float(x.flat[index]))**2
            exact = y*(a/2 + y*(b/4 + y*c/6))
            try:
                value.flat[index] = float(exact)
            except OverflowError:
                value.flat[index] = np.inf if exact > 0 else -np.inf
        return value[()] if value.ndim == 0 else value

    def force(self, x):
        return -(self.a*x+self.b*x**3+self.c*x**5)

    def energy(self, x, v):
        return 0.5*np.asarray(v)**2+self.potential(x)


def _stationary_sq_roots(a, b, c, number, square_root):
    """Pairs (y, d(a + by + cy²)/dy) for the real roots in y."""
    if c == 0:
        return [(number(-a / b), number(b))] if b else []
    discriminant = b*b - 4*a*c
    if discriminant < 0:
        return []
    if discriminant == 0:
        return [(number(-b / (2*c)), number(0))]

    # Take the sign with exact coefficients; scale before converting, then
    # recover the second root by Vieta to avoid subtracting close numbers.
    scale = max(abs(a), abs(b), abs(c))
    d = square_root(number(discriminant / scale**2))
    sign = 1 if b >= 0 else -1
    q = -(number(b / scale) + sign*d)/2
    return [(q / number(c / scale), -sign*d*number(scale)),
            (number(a / scale) / q, sign*d*number(scale))]


def _stationary_geometry(model, tol=1e-10):
    a, b, c = (Fraction(float(x)) for x in (model.a, model.b, model.c))
    try:
        roots = _stationary_sq_roots(a, b, c, float, sqrt)
        if not all(isfinite(y) and isfinite(slope) for y, slope in roots):
            raise OverflowError("Stationary geometry outside binary64 scale")
        positive = [(sqrt(y), 2*y*slope) for y, slope in roots if y > tol]
    except (OverflowError, ZeroDivisionError):
        # An extreme coefficient ratio may underflow in the scaled quadratic
        # even when a physical stationary position is finite.
        with localcontext() as ctx:
            ctx.prec = 40

            def decimal(value):
                return Decimal(value.numerator) / Decimal(value.denominator)

            roots = _stationary_sq_roots(a, b, c, decimal, ctx.sqrt)
            positive = [(float(y.sqrt()), float(2*y*slope)) for y, slope in roots if y > tol]

    points = [(0.0, float(a))]
    for x, curvature in positive:
        # At a stationary point with x != 0, V''(x) = 2 y (b + 2 c y).
        points.extend([(-x, curvature), (x, curvature)])
    return sorted(points)


def stationary_points(model, tol=1e-10):
    """Return sorted stationary positions, including zero and each repeated root once.

    Nonzero points are retained when their squared position exceeds tol.
    """
    return np.array([x for x, _ in _stationary_geometry(model, tol)])


def classify_energy(model, x0, v0, tol=1e-9):
    """Label energy surface relative to equilibrium barrier energies.

    Labels are descriptive; `bounded_component` is not inferred merely from E.
    """
    E=float(model.energy(x0,v0))
    if not isfinite(E):
        raise ValueError("Initial energy outside the finite range of double precision")
    eq = _stationary_geometry(model)
    barriers, minima, degenerate = [], [], []
    for x, curvature in eq:
        target = barriers if curvature < -tol else minima if curvature > tol else degenerate
        eb = float(model.potential(x))
        if not isfinite(eb):
            raise ValueError("Stationary energy outside the finite range of double precision")
        target.append((x, eb))
    at_barrier = any(abs(E-eb) < tol*(1+abs(eb)) for _, eb in barriers + degenerate)
    return {'energy': E, 'equilibria': [x for x, _ in eq], 'maxima': barriers,
            'minima': minima, 'degenerate': degenerate, 'on_separatrix_energy': at_barrier}


def turning_points(model, E, tol=1e-9):
    """All real turning points V(x)=E, without assuming confinement."""
    E=float(E)
    coef=[model.c/6,model.b/4,model.a/2,-E]
    while len(coef)>1 and abs(coef[0])<1e-14: coef.pop(0)
    yy=np.roots(coef) if len(coef)>1 else np.array([])
    xs=[]
    for y in yy:
        if abs(y.imag) < tol*(1+abs(y.real)) and y.real >= -tol:
            z=np.sqrt(max(0.,float(y.real)))
            xs.extend([-z,z] if z>tol else [0.])
    # Newton polishing of each root of V(x) - E (skipped at near-double roots)
    polished=[]
    for z in xs:
        for _ in range(3):
            dv=-float(model.force(z))
            if abs(dv) < 1e-8*(1+abs(z)): break
            z=z-(float(model.potential(z))-E)/dv
        polished.append(float(z))
    out=[]
    for z in sorted(polished):
        if not out or abs(z-out[-1]) > 1e-10*(1+abs(z)): out.append(z)
    return np.array(out)


_EPS = np.finfo(float).eps


def _gap_positive(model, E, x):
    """True if E - V(x) > 0 beyond the rounding error of its evaluation."""
    x = float(x)
    scale = abs(E) + abs(model.a)*x*x/2 + abs(model.b)*x**4/4 + abs(model.c)*x**6/6
    return E - float(model.potential(x)) > 64*_EPS*scale


def _admissible_interval(model, x0, E, tol):
    """Consecutive turning points [x-, x+] enclosing x0 with E > V inside, or None."""
    roots = turning_points(model, E)
    pairs = [(float(l), float(r)) for l, r in zip(roots[:-1], roots[1:])
             if _gap_positive(model, E, (l + r)/2) and l - tol <= x0 <= r + tol]
    return pairs[0] if len(pairs) == 1 else None


def _near_minimum(model, x0, E, tol):
    """Non-degenerate minimum x_m whose energy lies within tol of E, with x0 inside the
    corresponding tiny well (oscillation amplitude below about sqrt(tol)); else None."""
    for xm, k in _stationary_geometry(model):
        if k > tol:
            Vm = float(model.potential(xm))
            if -tol*(1+abs(Vm)) <= E - Vm <= tol*(1+abs(Vm)) and \
                    abs(x0 - xm) <= 2*np.sqrt(2*tol*(1+abs(Vm))/k) + tol:
                return float(xm), float(k)
    return None


def _deflate(model, p, left, right, E):
    """Quotient Q of E - V(x) (coefficients p) by (x - left)(right - x).

    Forward deflation (np.polydiv) is used unless it fails to reproduce E - V on the
    interval to 1e-12 relative to the size of its terms.  This happens for wide wells
    (|x| >> 1, e.g. a=1, b=-1, c=1e-8), where forward deflation loses the low-order
    coefficients; backward deflation (from the constant term) is then used if it is
    more accurate (v0.2.1).
    """
    D = np.array([-1., left + right, -left*right])
    qf = np.polydiv(p, D)[0]
    if p[-1] == 0 or left*right == 0:
        return qf
    xs = (left + right)/2 + (right - left)/2*np.sin(np.linspace(-np.pi/2, np.pi/2, 33)[1:-1])
    scale = abs(E) + abs(model.a)*xs**2/2 + abs(model.b)*xs**4/4 + abs(model.c)*xs**6/6

    def resid(q):
        d = np.polyval(q, xs)*(xs - left)*(right - xs) - (E - model.potential(xs))
        return float(np.max(np.abs(d)/scale))
    rf = resid(qf)
    if rf <= 1e-12:
        return qf
    qb = np.polydiv(p[::-1], D[::-1])[0][::-1]
    return qb if resid(qb) < rf else qf


def _resolved_period_quotient(q, left, right):
    """Screen the endpoint and interior minima of the computed quotient."""
    nodes = [left, right]
    try:
        roots = np.roots(np.polyder(q))
    except (ValueError, np.linalg.LinAlgError):
        return False
    for z in roots:
        if left < z.real < right and abs(z.imag) <= sqrt(_EPS) * (1 + abs(z.real)):
            nodes.append(float(z.real))
    q_at_nodes = np.polyval(q, nodes)
    scale = np.polyval(np.abs(q), np.abs(nodes))
    return np.all(np.isfinite(scale) & (q_at_nodes > 64*_EPS*scale))


def period(model, x0, v0, *, tol=1e-9):
    """Period in the connected bounded admissible component containing x0.

    Returns np.inf for equilibrium and raises ValueError for non-closed
    or separatrix trajectories. An unresolved or unconverged double-precision
    quadrature raises ValueError with a referral to period_mp.
    """
    if abs(v0)<tol and abs(float(model.force(x0)))<tol:
        return float('inf')
    geometry = classify_energy(model, x0, v0, tol)
    if geometry['on_separatrix_energy']:
        raise ValueError('Separatrix energy: infinite period')
    E = geometry['energy']
    comp = _admissible_interval(model, x0, E, tol)
    if comp is None:
        nm = _near_minimum(model, x0, E, tol)
        if nm is not None:
            # amplitude below ~sqrt(tol): harmonic limit, relative error O(tol)
            return float(2*np.pi/np.sqrt(nm[1]))
        raise ValueError('Initial state is not in a verified bounded component')
    left, right = comp
    center, half = (right+left)/2, (right-left)/2
    stationary = geometry['equilibria']
    for xt in (left, right):
        if any(abs(xt-xe) <= tol*(1+abs(xe)) for xe in stationary):
            raise ValueError('Degenerate (multiple) turning point: period not finite')

    # E - V = (x - left)(right - x) Q(x), with dx/sqrt(2(E-V)) = dtheta/sqrt(2 Q).
    p = np.array([-model.c/6, 0., -model.b/4, 0., -model.a/2, 0., E])
    while len(p) > 1 and p[0] == 0:
        p = p[1:]
    q = _deflate(model, p, left, right, E)
    resolution_error = 'Period unresolved in double precision; use period_mp with the same inputs'
    if not np.all(np.isfinite(q)) or not _resolved_period_quotient(q, left, right):
        raise ValueError(resolution_error)

    def integrand(theta):
        qx = float(np.polyval(q, center + half*np.sin(theta)))
        if qx <= 0:
            raise ValueError(resolution_error)
        return 1.0/np.sqrt(2*qx)

    # Split at potential extrema near separatrices. QUADPACK adds a message
    # to full_output when a call does not converge to the requested tolerance.
    cuts = sorted(float(np.arcsin((xe-center)/half)) for xe in stationary if left < xe < right)
    nodes = [-np.pi/2] + cuts + [np.pi/2]
    total = 0.0
    for lo, hi in zip(nodes[:-1], nodes[1:]):
        if hi > lo:
            result = quad(integrand, lo, hi, epsabs=1e-14, epsrel=1e-13, limit=400, full_output=1)
            if len(result) != 3:
                raise ValueError(resolution_error)
            total += result[0]
    if not isfinite(total) or total <= 0:
        raise ValueError(resolution_error)
    return float(2*total)


def simulate(model, x0, v0, t_end=40., n=2001, rtol=1e-10, atol=1e-12):
    """Integrate IVP using adaptive DOP853; returns scipy OdeResult."""
    if t_end<=0 or n<2: raise ValueError('t_end must be positive and n>=2')
    def f(t,z): return [z[1],float(model.force(z[0]))]
    t=np.linspace(0,t_end,n)
    result=solve_ivp(f,(0.,t_end),[x0,v0],method='DOP853',t_eval=t,rtol=rtol,atol=atol)
    if not result.success: raise RuntimeError(result.message)
    return result


def energy_drift(model, solution):
    energies=model.energy(solution.y[0],solution.y[1])
    E0=float(energies[0]); return float(np.max(np.abs(energies-E0))/(1+abs(E0)))


def motion_type(model, x0, v0, tol=1e-9):
    """Classify the motion through (x0, v0) from the energy geometry.

    Returns one of ``'equilibrium'``, ``'separatrix'``, ``'bounded'`` or
    ``'unbounded'``.  ``'bounded'`` means that x0 lies in an admissible interval
    [x-, x+] bounded by two consecutive simple turning points (closed orbit), or in
    an oscillation around a non-degenerate minimum too small to resolve the turning
    points in double precision (amplitude below about sqrt(tol)).
    """
    if abs(v0) < tol and abs(float(model.force(x0))) < tol:
        return 'equilibrium'
    if classify_energy(model, x0, v0, tol)['on_separatrix_energy']:
        return 'separatrix'
    E = float(model.energy(x0, v0))
    if _admissible_interval(model, x0, E, tol) is not None or _near_minimum(model, x0, E, tol) is not None:
        return 'bounded'
    return 'unbounded'


def period_mp(model, x0, v0, *, dps=50, tol=1e-9):
    """Arbitrary-precision period (mpmath), the Python analogue of
    ``NOLPeriod[..., WorkingPrecision -> dps]``.

    Coefficients and initial data may be given as ``str`` (e.g. ``'0.4'``, read as an
    exact decimal), ``fractions.Fraction``, ``int`` or ``float`` (taken as the exact
    binary value).  Uses the same factored integrand as :func:`period`, with 20 guard
    digits, and tanh-sinh quadrature split at the stationary points.  Returns an
    ``mpmath.mpf`` correct to about ``dps`` digits.
    """
    import mpmath as mp
    ctx = mp.mp.clone()
    ctx.dps = int(dps) + 20

    def num(v):
        if isinstance(v, Fraction):
            return ctx.mpf(v.numerator) / v.denominator
        return ctx.mpf(v) if isinstance(v, (str, int)) else ctx.mpf(float(v))
    coefficients = tuple(Fraction(v if isinstance(v, (str, int, Fraction)) else float(v))
                         for v in (model.a, model.b, model.c))
    a, b, c = map(num, coefficients)
    x0, v0 = num(x0), num(v0)

    def V(x):
        return a * x ** 2 / 2 + b * x ** 4 / 4 + c * x ** 6 / 6
    E = v0 ** 2 / 2 + V(x0)
    if abs(v0) < tol and abs(a * x0 + b * x0 ** 3 + c * x0 ** 5) < tol:
        return ctx.inf

    def real_sq_roots(coeffs):                     # coefficients in y = x^2, highest first
        while len(coeffs) > 1 and coeffs[0] == 0:
            coeffs = coeffs[1:]
        if len(coeffs) < 2:
            return []
        ys = ctx.polyroots(coeffs, maxsteps=200, extraprec=2 * ctx.prec)
        xs = []
        for y in ys:
            y = ctx.mpc(y)
            if abs(y.imag) <= ctx.mpf(10) ** (-dps) * (1 + abs(y.real)) and y.real >= 0:
                z = ctx.sqrt(y.real)
                xs.extend([-z, z] if z > 0 else [ctx.zero])
        return sorted(xs)
    equilibria = [(ctx.zero, a)]
    for y, slope in _stationary_sq_roots(*coefficients, num, ctx.sqrt):
        if y > 0:
            x = ctx.sqrt(y)
            equilibria.extend([(-x, 2*y*slope), (x, 2*y*slope)])
    eq = sorted(x for x, _ in equilibria)
    for xe, curvature in equilibria:
        if curvature <= tol and abs(E - V(xe)) <= tol * (1 + abs(V(xe))):
            raise ValueError('Separatrix energy: infinite period')
    roots = real_sq_roots([c / 6, b / 4, a / 2, -E])
    gap_tol = ctx.mpf(10) ** (-int(dps) - 10)
    pairs = [(l, r) for l, r in zip(roots[:-1], roots[1:])
             if E - V((l + r) / 2) > gap_tol * (1 + abs(E)) and l - tol <= x0 <= r + tol]
    if len(pairs) != 1:
        raise ValueError('Initial state is not in a verified bounded component')
    left, right = pairs[0]
    ctr, half = (left + right) / 2, (right - left) / 2
    # E - V(x) = (x-left)(right-x) Q(x): synthetic division of -V + E by -(x-left)(x-right)
    P = [-c / 6, 0, -b / 4, 0, -a / 2, 0, E]
    while P[0] == 0:
        P = P[1:]
    D = [-1, left + right, -left * right]
    Q = list(P)
    quot = []
    for i in range(len(P) - 2):
        coef = Q[i] / D[0]
        quot.append(coef)
        for j in range(3):
            Q[i + j] -= coef * D[j]

    def integrand(th):
        x = ctr + half * ctx.sin(th)
        return 1 / ctx.sqrt(2 * ctx.polyval(quot, x))
    nodes = [-ctx.pi / 2] + sorted(ctx.asin((xe - ctr) / half) for xe in eq if left < xe < right) + [ctx.pi / 2]
    T = 2 * ctx.quad(integrand, nodes)
    out = mp.mp.clone(); out.dps = int(dps)
    return out.mpf(T)
