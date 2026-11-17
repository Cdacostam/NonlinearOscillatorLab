"""Coverage and accuracy campaigns for the closed-form solver (Python port).

Usage: python reproduce/closed_form_campaign.py [--workers N]
The default run reads the fixed round-trippable binary64 inputs in
data/inputs/closed_form_campaign.csv. Use --generate to exercise the sampling design.

Writes reproduce/outputs/closed_form_campaign_py.csv and a JSON summary.

Campaigns (all seeds fixed, see SEED):
  R200      the 200 cases of the published Wolfram campaign (tables/closed_form_campaign.csv)
  RWIDE     N random problems, a,b ~ U[-12,12], c ~ U[-1,1], x0,v0 ~ U[-2,2]
            (includes c < 0, i.e. unbounded and escaping motion)
  SEP       energies E = E_s + s*|E_s|*10^-k (k = 1..12, s = -1,+1) near the barriers of
            three multi-well potentials, plus E = E_s + s*10^-k when E_s = 0
  REP       degenerate parameters: b^2 = 4ac (coalescing equilibria, inflection points),
            states at minima (equilibria) and near-repeated turning points
  UNB       escape boundaries of unbounded potentials (c < 0, or c = 0 and b < 0):
            E = E_barrier*(1 +- 10^-k) and far above the barrier

Error metric (identical to the manuscript):
  e_x = max_{t in grid} |x_cf(t) - x_ref(t)| on the uniform grid of 401 points in
  [0, t_w], t_w = min(20, 0.9 t_esc), where t_esc is the time at which the reference
  leaves |x| <= 1e3 (t_esc = inf for bounded motion).
Reference:
  x_ref from SciPy DOP853 with rtol = 2.5e-14, atol = 1e-16.  Every case with
  e_x > 1e-10 is recomputed against a 25-digit mpmath Taylor-series reference started from the
  exact binary inputs
  (mpmath.odefun); the reported e_x is then the one against that reference and
  the column ref = "mp25".
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import multiprocessing as mpc
import platform
import sys
import time
from pathlib import Path

import mpmath
import numpy as np
import scipy
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parents[1]                       # package root
OUTPUTS = Path(__import__("os").environ.get("NOL_OUTPUTS", ROOT / "reproduce" / "outputs"))
sys.path.insert(0, str(ROOT / "code_py" / "src"))
from nonlinear_oscillator_lab import Model, classify_energy, stationary_points  # noqa: E402
from nonlinear_oscillator_lab.closed_form import exact_solution, ClosedFormNotAvailable  # noqa: E402
from nonlinear_oscillator_lab.core import motion_type  # noqa: E402

SEED = 20261103
OUT = OUTPUTS
CANONICAL = ROOT / "data" / "inputs" / "closed_form_campaign.csv"


def mp_reference(a, b, c, x0, v0, ts, dps=25):
    import mpmath as mp
    ctx = mp.mp.clone(); ctx.dps = dps
    # v0.2.1: exact binary values of the inputs and grid times. v0.2.0 used mpf(repr(float(v))), the
    # shortest decimal string, which near a barrier energy moved the reference by up to 8.1e-7.
    a, b, c, x0, v0 = (ctx.mpf(float(v)) for v in (a, b, c, x0, v0))
    F = ctx.odefun(lambda t, z: [z[1], -(a * z[0] + b * z[0] ** 3 + c * z[0] ** 5)], 0, [x0, v0])
    return np.array([float(F(ctx.mpf(float(t)))[0]) for t in ts])


def reference(a, b, c, x0, v0):
    def f(t, z):
        return [z[1], -(a * z[0] + b * z[0] ** 3 + c * z[0] ** 5)]

    def esc(t, z):
        return abs(z[0]) - 1e3
    esc.terminal = True
    pre = solve_ivp(f, (0, 20), [x0, v0], method="DOP853", rtol=1e-12, atol=1e-14, events=esc)
    t_esc = float(pre.t_events[0][0]) if len(pre.t_events[0]) else np.inf
    tw = min(20.0, 0.9 * t_esc)
    ts = np.linspace(0, tw, 401)
    ref = solve_ivp(f, (0, tw), [x0, v0], method="DOP853", t_eval=ts, rtol=2.5e-14, atol=1e-16)
    if not ref.success or len(ref.t) != len(ts):
        raise RuntimeError(f"Incomplete DOP853 reference: {ref.message}")
    return ts, ref.y[0], t_esc


def run_case(args):
    group, a, b, c, x0, v0, note = args
    m = Model(a, b, c)
    E = float(m.energy(x0, v0))
    mt = motion_type(m, x0, v0)
    on_sep = classify_energy(m, x0, v0)["on_separatrix_energy"]
    t0 = time.perf_counter()
    try:
        s = exact_solution(a, b, c, x0, v0)
        ok, rep, reason = True, s.representation, ""
    except ClosedFormNotAvailable as exc:
        s, ok, rep, reason = None, False, "", str(exc)
    except Exception as exc:                                   # record, never hide
        s, ok, rep, reason = None, False, "", "exception: %s" % type(exc).__name__
    build = time.perf_counter() - t0
    err, refkind, t_esc = np.nan, "", np.nan
    if ok:
        ts, xr, t_esc = reference(a, b, c, x0, v0)
        with np.errstate(all="ignore"):
            xc = s.x(ts)
        err = float(np.max(np.abs(xc - xr))) if np.all(np.isfinite(xc)) else np.inf
        refkind = "dop853"
        if not err <= 1e-10:
            try:
                xr2 = mp_reference(a, b, c, x0, v0, ts)
                err = float(np.max(np.abs(xc - xr2))) if np.all(np.isfinite(xc)) else np.inf
                refkind = "mp25"
            except Exception:
                refkind = "dop853 (mp failed)"
        if not np.isfinite(err) and ok:
            reason = "non-finite values on the comparison window"
    return dict(group=group, note=note, a=a, b=b, c=c, x0=x0, v0=v0, E=E, motion=mt,
                on_separatrix=bool(on_sep), ok=ok, representation=rep, reason=reason,
                err=err, ref=refkind, t_esc=t_esc, build_s=build)


# --------------------------------------------------------------------- case sets

def cases_r200():
    rows = list(csv.DictReader(open(ROOT / "data" / "wolfram" / "closed_form_campaign.csv")))
    return [("R200", float(r["a"]), float(r["b"]), float(r["c"]), float(r["x0"]), float(r["v0"]),
             r["rep"]) for r in rows]


def cases_rwide(n):
    rng = np.random.default_rng(SEED)
    out = []
    for _ in range(n):
        a, b = rng.uniform(-12, 12, 2); c = rng.uniform(-1, 1)
        x0, v0 = rng.uniform(-2, 2, 2)
        out.append(("RWIDE", float(a), float(b), float(c), float(x0), float(v0), ""))
    return out


def _v0_for(m, x0, E):
    d = 2 * (E - float(m.potential(x0)))
    return float(np.sqrt(d)) if d >= 0 else None


def cases_sep(kmax=12):
    out = []
    pots = {  # name: (a, b, c, barrier x, wells' sample x0)
        "double well (c=0)": (-1.0, 1.0, 0.0, 0.0, [1.0, 0.0]),
        "triple well": (1.0, -4.0, 3.0, 3 ** -0.5, [0.0, 1.0]),
        "cq double well": (-1.0, -0.5, 0.5, 0.0, [1.4, 0.0]),
    }
    for name, (a, b, c, xb, x0s) in pots.items():
        m = Model(a, b, c); Es = float(m.potential(xb))
        for k in range(1, kmax + 1):
            for sgn in (-1, 1):
                E = Es + sgn * (abs(Es) if Es else 1.0) * 10.0 ** -k
                for x0 in x0s:
                    v0 = _v0_for(m, x0, E)
                    if v0 is not None:
                        out.append(("SEP", a, b, c, x0, v0, "%s, E-Es=%+.0e" % (name, E - Es)))
    return out


def cases_rep():
    out = []
    # b^2 = 4ac: V'(x) = x (a + b x^2 + c x^4) has double roots (inflection equilibria)
    for a, b, c in [(1.0, -2.0, 1.0), (2.0, -4.0, 2.0), (0.5, -1.0, 0.5)]:
        for x0, v0 in [(0.3, 0.0), (0.5, 0.2), (0.9, 0.0), (1.2, 0.0), (0.0, 0.8), (0.0, 1.2)]:
            out.append(("REP", a, b, c, x0, v0, "b^2=4ac"))
    # equilibria (E at a minimum) and almost-equilibria
    for a, b, c in [(1.0, 0.4, 0.2), (-1.0, 1.0, 0.0), (1.0, -4.0, 3.0)]:
        m = Model(a, b, c)
        for xe in stationary_points(m):
            if a + 3 * b * xe ** 2 + 5 * c * xe ** 4 > 0:
                for dv in (0.0, 1e-8, 1e-4):
                    out.append(("REP", a, b, c, float(xe), dv, "near minimum, v0=%g" % dv))
    # c = 0 and b = 0 limits (harmonic and pure quintic / pure cubic)
    for a, b, c in [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 1.0)]:
        out.append(("REP", a, b, c, 0.7, 0.3, "degenerate coefficients"))
    return out


def cases_unb(kmax=10):
    out = []
    for a, b, c in [(1.0, 1.0, -0.5), (1.0, -1.0, 0.0), (-1.0, 2.0, -0.6), (1.0, 0.5, -1.0)]:
        m = Model(a, b, c)
        maxima = [x for x in stationary_points(m) if a + 3 * b * x * x + 5 * c * x ** 4 < 0 and x > 0]
        if not maxima:
            continue
        xb = max(maxima); Eb = float(m.potential(xb))
        for k in range(1, kmax + 1, 3):
            for sgn in (-1, 1):
                E = Eb * (1 + sgn * 10.0 ** -k)
                for x0 in (0.0, 0.5 * xb):
                    v0 = _v0_for(m, x0, E)
                    if v0 is not None:
                        out.append(("UNB", a, b, c, x0, v0, "E/Eb-1=%+.0e" % (E / Eb - 1)))
        for x0, v0 in [(1.5 * xb, 0.0), (1.5 * xb, -0.5), (0.0, 3.0 * np.sqrt(abs(Eb)) + 1)]:
            out.append(("UNB", a, b, c, x0, v0, "beyond barrier / far above"))
    return out


def _run_indexed(item):
    i, case = item
    r = run_case(case)
    r["case_id"] = i + 1
    return r


def read_cases(path=CANONICAL):
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if [int(row["case_id"]) for row in rows] != list(range(1, len(rows) + 1)):
        raise ValueError("Campaign case IDs must be consecutive and start at 1")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--generate", action="store_true", help="generate a new sample from the design")
    ap.add_argument("--quick", action="store_true", help="generate a small RWIDE sample")
    ap.add_argument("--n", type=int, default=1000, help="RWIDE size with --generate")
    ap.add_argument("--workers", type=int, default=min(4, mpc.cpu_count()))
    args = ap.parse_args()
    if args.workers < 1 or args.n < 1:
        ap.error("--workers and --n must be positive")
    generated = args.generate or args.quick
    if generated:
        n = 100 if args.quick else args.n
        cases = cases_r200() + cases_rwide(n) + cases_sep() + cases_rep() + cases_unb()
    else:
        cases = [(r["group"], *(float(r[k]) for k in ("a", "b", "c", "x0", "v0")), r["note"])
                 for r in read_cases()]
        n = sum(c[0] == "RWIDE" for c in cases)

    OUT.mkdir(parents=True, exist_ok=True)
    # Resume only if inputs, solver code and numerical environment match.
    part = OUT / "closed_form_campaign_py.partial.jsonl"
    signature = OUT / "closed_form_campaign_py.partial.sha256"
    sources = sorted((ROOT / "code_py" / "src" / "nonlinear_oscillator_lab").glob("*.py"))
    canonical_hash = "" if generated else hashlib.sha256(CANONICAL.read_bytes()).hexdigest()
    environment = (platform.platform(), sys.version, np.__version__, scipy.__version__, mpmath.__version__)
    digest = hashlib.sha256(repr((cases, environment, canonical_hash)).encode())
    for path in [Path(__file__), *sources]:
        digest.update(path.read_bytes())
    stamp = digest.hexdigest()
    if part.exists():
        if not signature.exists() or signature.read_text().strip() != stamp:
            raise RuntimeError(f"Checkpoint differs from this run; use a fresh NOL_OUTPUTS path: {part}")
        with part.open() as f:
            done = {r["case_id"]: r for r in map(json.loads, f)}
    else:
        signature.write_text(stamp + "\n")
        done = {}
    todo = [(i, case) for i, case in enumerate(cases) if i + 1 not in done]
    print(f"cases {len(cases)}, already done {len(done)}", flush=True)
    t0 = time.time()
    with mpc.Pool(args.workers) as pool, part.open("a") as fh:
        for k, r in enumerate(pool.imap_unordered(_run_indexed, todo, chunksize=1)):
            fh.write(json.dumps(r, default=float) + "\n")
            fh.flush()
            done[r["case_id"]] = r
            if k % 50 == 0:
                print(f"{len(done)}/{len(cases)}  {time.time() - t0:.0f}s", flush=True)
    rows = [done[i] for i in range(1, len(cases) + 1)]
    with (OUT / "closed_form_campaign_py.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, ["case_id", *(k for k in rows[0] if k != "case_id")])
        w.writeheader()
        w.writerows(rows)
    summary = {}
    for g in ["R200", "RWIDE", "SEP", "REP", "UNB"]:
        R = [r for r in rows if r["group"] == g]
        okR = [r for r in R if r["ok"]]
        e = np.array([r["err"] for r in okR if np.isfinite(r["err"])])
        reasons = {}
        for r in R:
            if not r["ok"]:
                key = "%s | motion=%s" % (r["reason"], r["motion"])
                reasons[key] = reasons.get(key, 0) + 1
        by_motion = {}
        for r in R:
            d = by_motion.setdefault(r["motion"], [0, 0]); d[0] += 1; d[1] += int(r["ok"])
        summary[g] = dict(
            cases=len(R), represented=len(okR),
            failure_reasons=reasons, by_motion={k: {"cases": v[0], "represented": v[1]} for k, v in by_motion.items()},
            err_median=float(np.median(e)) if len(e) else None,
            err_p95=float(np.percentile(e, 95)) if len(e) else None,
            err_max=float(e.max()) if len(e) else None,
            frac_below_1e8=float(np.mean(e < 1e-8)) if len(e) else None,
            nonfinite=int(sum(1 for r in okR if not np.isfinite(r["err"]))),
            mp25_rechecks=int(sum(1 for r in okR if r["ref"] == "mp25")),
            build_median_ms=float(1e3 * np.median([r["build_s"] for r in R])))
        if g == "R200":
            summary[g]["labels_agree_with_wolfram"] = int(sum(r["representation"] == r["note"] for r in R))
    summary["_meta"] = dict(seed=SEED, rwide_n=n, wall_s_last_run=time.time() - t0, python=sys.version.split()[0],
                            numpy=np.__version__, scipy=scipy.__version__, mpmath=mpmath.__version__,
                            platform=platform.platform(), workers=args.workers, canonical_sha256=canonical_hash)
    with (OUT / "closed_form_campaign_py_summary.json").open("w") as f:
        json.dump(summary, f, indent=2)
        f.write("\n")
    print(json.dumps(summary, indent=1, default=str))


if __name__ == "__main__":
    main()
