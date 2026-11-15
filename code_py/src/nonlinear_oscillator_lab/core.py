"""Energy methods for x'' + a*x + b*x**3 + c*x**5 = 0.

No assertion of analytic closed-form trajectories is made. Period is computed
by numerical quadrature of the factored, singularity-free integrand
obtained from E - V(x) = (x - x_-)(x_+ - x) Q(x) and x = mid + half*sin(theta).
"""
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp, quad

@dataclass(frozen=True)
class Model:
    a: float = 1.0
    b: float = 0.0
    c: float = 0.0

    def potential(self, x):
        x = np.asarray(x)
        return 0.5*self.a*x*x + 0.25*self.b*x**4 + self.c*x**6/6

    def force(self, x):
        return -(self.a*x+self.b*x**3+self.c*x**5)

    def energy(self, x, v):
        return 0.5*np.asarray(v)**2+self.potential(x)


def stationary_points(model, tol=1e-10):
    """Return sorted real roots of dV/dx, including zero."""
    roots=[0.0]
    # c y^2 + b y + a = 0, with y=x^2>0
    coefficients=[model.c,model.b,model.a]
    while len(coefficients)>1 and abs(coefficients[0])<tol:
        coefficients.pop(0)
    if len(coefficients)>1:
        for y in np.roots(coefficients):
            if abs(y.imag)<tol and y.real>tol:
                z=float(np.sqrt(y.real));roots.extend([-z,z])
    return np.array(sorted(set(round(x,12) for x in roots)))


def classify_energy(model, x0, v0, tol=1e-9):
    """Label energy surface relative to equilibrium barrier energies.

    Labels are descriptive; `bounded_component` is not inferred merely from E.
    """
    E=float(model.energy(x0,v0))
    eq=stationary_points(model)
    curv=lambda x: model.a+3*model.b*x*x+5*model.c*x**4
    barriers=[(float(x),float(model.potential(x))) for x in eq if curv(x) < -tol]
    minima=[(float(x),float(model.potential(x))) for x in eq if curv(x) > tol]
    # degenerate equilibria (V''=0, e.g. inflection points when b^2 = 4ac) are
    # approached asymptotically like saddles: their energy level is also critical
    degenerate=[(float(x),float(model.potential(x))) for x in eq if abs(curv(x)) <= tol]
    at_barrier=any(abs(E-eb)<tol*(1+abs(eb)) for _,eb in barriers+degenerate)
    return {'energy': E,'equilibria':eq.tolist(),'maxima':barriers,'minima':minima,
            'degenerate':degenerate,'on_separatrix_energy': bool(at_barrier)}


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
    for xm in stationary_points(model):
        k = model.a + 3*model.b*xm*xm + 5*model.c*xm**4
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


def period(model, x0, v0, *, tol=1e-9):
    """Period in the connected bounded admissible component containing x0.

    Returns np.inf for equilibrium and raises ValueError for non-closed
    or separatrix trajectories. Period=2*integral dx/sqrt(2*(E-V(x))).
    """
    E=float(model.energy(x0,v0))
    if abs(v0)<tol and abs(float(model.force(x0)))<tol:
        return float('inf')
    if classify_energy(model,x0,v0,tol)['on_separatrix_energy']:
        raise ValueError('Separatrix energy: infinite period')
    comp=_admissible_interval(model,x0,E,tol)
    if comp is None:
        nm=_near_minimum(model,x0,E,tol)
        if nm is not None:
            # amplitude below ~sqrt(tol): harmonic limit, relative error O(tol)
            return float(2*np.pi/np.sqrt(nm[1]))
        raise ValueError('Initial state is not in a verified bounded component')
    pairs=[comp]
    left,right=pairs[0]; center=(right+left)/2; half=(right-left)/2
    # simple turning points are required: a turning point that coincides with an
    # equilibrium is a multiple root of E - V (infinite period)
    for xt in (left,right):
        if any(abs(xt-xe) <= tol*(1+abs(xe)) for xe in stationary_points(model)):
            raise ValueError('Degenerate (multiple) turning point: period not finite')
    # Factor E - V(x) = (x - left)(right - x) Q(x).  With x = center + half*sin(theta),
    # dx/sqrt(2(E-V)) = dtheta/sqrt(2 Q(x)): the integrand is smooth and nonsingular.
    p=np.array([-model.c/6,0.,-model.b/4,0.,-model.a/2,0.,E])
    while len(p)>1 and p[0]==0: p=p[1:]
    q=_deflate(model,p,left,right,E)
    def integrand(theta):
        qx=float(np.polyval(q,center+half*np.sin(theta)))
        if qx<=0:
            raise ValueError('Degenerate turning point: period quadrature not regular')
        return 1.0/np.sqrt(2*qx)
    # Near a separatrix the integrand peaks where the orbit crosses a barrier
    # (a local maximum of V); split the quadrature there.
    cuts=sorted(float(np.arcsin((xe-center)/half)) for xe in stationary_points(model)
                if left<xe<right)
    nodes=[-np.pi/2]+cuts+[np.pi/2]
    val=sum(quad(integrand,lo,hi,epsabs=1e-14,epsrel=1e-13,limit=400)[0]
            for lo,hi in zip(nodes[:-1],nodes[1:]) if hi>lo)
    if not np.isfinite(val) or val<=0:
        raise ValueError('Period quadrature failed')
    return float(2*val)


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
    from fractions import Fraction
    ctx = mp.mp.clone()
    ctx.dps = int(dps) + 20

    def num(v):
        if isinstance(v, Fraction):
            return ctx.mpf(v.numerator) / v.denominator
        return ctx.mpf(v) if isinstance(v, (str, int)) else ctx.mpf(float(v))
    a, b, c = (num(getattr(model, k)) for k in ('a', 'b', 'c'))
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
    eq = sorted(set(real_sq_roots([c, b, a]) + [ctx.zero]))
    for xe in eq:
        if a + 3 * b * xe ** 2 + 5 * c * xe ** 4 <= tol and abs(E - V(xe)) <= tol * (1 + abs(V(xe))):
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
