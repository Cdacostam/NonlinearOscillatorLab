"""Re-reference the 25-digit mpmath rows of the archived closed-form campaign (v7, finding H1).

The script closed_form_campaign.py of release v0.2.0 built its 25-digit
reference from ctx.mpf(repr(float(v))), i.e. from the shortest decimal string of each double,
not from the binary value that exact_solution receives. Near a separatrix a one-ulp change of
the data moves the solution by up to 1e-6, so that reference solves a different problem.

This script keeps everything else identical (cases, closed-form solver, grid, DOP853 window,
the rule "recheck with mp25 when e_x > 1e-10") and recomputes, for every row with ref = mp25:
  err_dec  : reference from decimal strings (must reproduce the archived err)
  err_bin  : reference from the exact binary inputs and grid times, ctx.mpf(float(v))
Usage: python reproduce/rereference_campaign.py [--workers N]
"""
from __future__ import annotations
import argparse, csv, json, multiprocessing as mpc, sys, time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "code_py" / "src")); sys.path.insert(0, str(HERE))
import closed_form_campaign as cfc  # noqa: E402
from nonlinear_oscillator_lab.closed_form import exact_solution  # noqa: E402

IN = HERE.parent / "data" / "results" / "closed_form_campaign_py.csv"
OUT = Path(__import__("os").environ.get("NOL_OUTPUTS", HERE / "outputs")) / "closed_form_campaign_py_exactref.csv"


def mp_reference(a, b, c, x0, v0, ts, dps=25, exact=True):
    import mpmath as mp
    ctx = mp.mp.clone(); ctx.dps = dps
    conv = (lambda v: ctx.mpf(float(v))) if exact else (lambda v: ctx.mpf(repr(float(v))))
    a, b, c, x0, v0 = (conv(v) for v in (a, b, c, x0, v0))
    F = ctx.odefun(lambda t, z: [z[1], -(a * z[0] + b * z[0] ** 3 + c * z[0] ** 5)], 0, [x0, v0])
    return np.array([float(F(conv(t))[0]) for t in ts])


def work(row):
    a, b, c, x0, v0 = (float(row[k]) for k in ("a", "b", "c", "x0", "v0"))
    s = exact_solution(a, b, c, x0, v0)
    ts, _, _ = cfc.reference(a, b, c, x0, v0)
    with np.errstate(all="ignore"):
        xc = s.x(ts)
    e_dec = float(np.max(np.abs(xc - mp_reference(a, b, c, x0, v0, ts, exact=False))))
    e_bin = float(np.max(np.abs(xc - mp_reference(a, b, c, x0, v0, ts, exact=True))))
    return dict(row, err_archived=row["err"], err_dec=e_dec, err_bin=e_bin)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=mpc.cpu_count())
    args = ap.parse_args()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(open(IN)))
    todo = [r for r in rows if r["ref"] == "mp25"]
    t0 = time.time()
    with mpc.Pool(args.workers) as pool:
        res = pool.map(work, todo, chunksize=1)
    by = {(r["group"], r["a"], r["b"], r["c"], r["x0"], r["v0"]): r for r in res}
    keys = list(rows[0].keys()) + ["err_archived", "err_dec", "err_bin"]
    out = []
    for r in rows:
        k = (r["group"], r["a"], r["b"], r["c"], r["x0"], r["v0"])
        if k in by:
            q = by[k]; q = dict(q); q["err"] = repr(q["err_bin"]); q["ref"] = "mp25bin"
        else:
            q = dict(r, err_archived=r["err"], err_dec="", err_bin="")
        out.append(q)
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, keys); w.writeheader(); w.writerows(out)
    rel = [abs(r["err_dec"] - float(r["err_archived"])) / max(float(r["err_archived"]), 1e-300) for r in res]
    print("rechecked %d mp25 rows in %.0f s; max rel. diff err_dec vs archived %.2e" % (len(res), time.time() - t0, max(rel)))


if __name__ == "__main__":
    main()
