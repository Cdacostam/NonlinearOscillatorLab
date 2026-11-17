"""Directed period and energy checks for the nearly degenerate family.

The benchmark reference uses monotone bisection in squared position and
analytically factors the resulting quadratic integrand, without package calls.
Run from the bundle root. Output: period_boundary_check.json in NOL_OUTPUTS
or reproduce/outputs/.
"""
import json
import math
import os
import platform
import sys
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy

from nonlinear_oscillator_lab import Model, classify_energy, motion_type, period, period_mp


def reference(c, v, ctx):
    c, v = ctx.mpf(c), ctx.mpf(v)
    energy = v*v/2

    def potential_sq(u):
        return u/2 - u*u/2 + c*u**3/6

    # dV(sqrt(u))/du = (1 - 2u + cu^2)/2 > 0 for c > 1.
    left, right = ctx.zero, ctx.mpf(2)
    assert potential_sq(left) < energy < potential_sq(right)
    for _ in range(4*ctx.dps):
        mid = (left + right)/2
        if potential_sq(mid) < energy:
            left = mid
        else:
            right = mid
    upper = (left + right)/2
    u_star = (3/c - upper)/2
    minimum = ctx.mpf(1)/2 - (u_star+upper)/2 + c*(u_star**2+u_star*upper+upper**2)/6
    width = ctx.sqrt(6*minimum/c)

    def integrand(theta):
        u = upper*ctx.sin(theta)**2
        return 4/ctx.sqrt(2*(minimum + c*(u-u_star)**2/6))

    # Split the independent quadrature around the narrow peak.
    cuts = [u_star + k*width for k in (-512,-64,-8,-1,0,1,8,64,512)]
    angles = [ctx.zero] + [ctx.asin(ctx.sqrt(u/upper)) for u in cuts if 0 < u < upper]
    return ctx.quad(integrand, angles + [ctx.pi/2])


def main():
    ctx = mp.mp.clone()
    ctx.dps = 85
    v = math.sqrt(1/3)
    rows = []
    for c in (math.nextafter(1., math.inf), 1 + 2**-49, 1 + 2**-40, 1.001):
        model = Model(1, -2, c)
        ref = reference(c, v, ctx)
        mp50 = period_mp(model, 0, v, dps=50)
        mp80 = period_mp(model, 0, v, dps=80)
        diff = max(abs(ctx.mpf(t) - ref)/ref for t in (mp50, mp80))
        assert diff < ctx.mpf("1e-12")
        assert motion_type(model, 0, v) == "bounded"
        try:
            double = period(model, 0, v)
        except ValueError as exc:
            assert "use period_mp" in str(exc) and c != 1.001
            double = str(exc)
        else:
            assert abs(ctx.mpf(double)-ref)/ref < ctx.mpf("1e-6")
        rows.append(dict(c_binary64=c.hex(), v0_binary64=v.hex(), reference=ctx.nstr(ref, 50),
                         period_mp50=ctx.nstr(mp50, 46), period_mp80=ctx.nstr(mp80, 77),
                         max_relative_error=float(diff), period_double=double))

    wide = classify_energy(Model(1, 0, -1e-220), 0, 1)["maxima"]
    expected = ctx.mpf(1)/(3*ctx.sqrt(-ctx.mpf(float(-1e-220))))
    assert len(wide) == 2 and all(abs(ctx.mpf(e)-expected)/expected < ctx.mpf("1e-14") for _, e in wide)
    result = dict(python=sys.version.split()[0], platform=platform.platform(), mpmath=mp.__version__,
                  numpy=np.__version__, scipy=scipy.__version__, periods=rows, wide_barriers=wide)
    folder = Path(os.environ.get("NOL_OUTPUTS", Path(__file__).resolve().parent / "outputs"))
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "period_boundary_check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
