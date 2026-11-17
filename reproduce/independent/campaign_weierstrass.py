"""Algorithm 1 on the 1436 initial states of the closed-form campaign.

The cases (a, b, c, x0, v0, motion class) are read from the published campaign file
    data/results/v0.2.2_campaign/closed_form_campaign_py.csv (R200, RWIDE, SEP, REP, UNB)
so that Algorithm 1 is run on exactly the states on which the packages' Algorithm 2 was evaluated.

Reference and metric (the initial states are shared with Algorithm 2):
  DOP853, rtol 2.5e-14, atol 1e-16, on t_j = j t_w/400, t_w = min(20, 0.9 t_esc), t_esc the first time
  |x| reaches 1e3;  hat e_x = max_j |x - x_ref| / max_j |x_ref|.
  Cases at a separatrix level whose error against DOP853 exceeds 1e-8 are compared instead with a 30-digit
  Taylor-series reference on t <= 3 (DOP853 leaves the separatrix; Supplement, Algorithm-1 paragraph).
  Other cases with error > 1e-8 use the 30-digit adaptive Taylor reference on the full window when finite.
  The Taylor result takes precedence regardless of the direction of the discrepancy. Both errors are written.

Usage
  python -I campaign_weierstrass.py [path/to/closed_form_campaign_py.csv] [--workers N] [--out FILE]
Writes reproduce/outputs/campaign_weierstrass.csv and prints the summary by motion class.
NumPy and SciPy only; does not import nonlinear_oscillator_lab.
"""
from __future__ import annotations

import argparse
import csv
import multiprocessing as mpc
import os
import sys
from pathlib import Path

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from algorithm_weierstrass import solve                                    # noqa: E402
from appendix_formulas import dop853_reference, taylor_reference, metrics  # noqa: E402


def run(row):
    a, b, c, x0, v0 = (float(row[k]) for k in ("a", "b", "c", "x0", "v0"))
    on_sep = row["on_separatrix"] == "True"
    motion = "separatrix" if on_sep else row["motion"]
    out = dict(group=row["group"], a=a, b=b, c=c, x0=x0, v0=v0, motion=motion, on_separatrix=on_sep,
               represented=False, family="", err=np.nan, ref="dop853", err_dop853=np.nan,
               err_taylor30=np.nan, alg2_represented=row["ok"])
    try:
        ts, xr, _ = dop853_reference(a, b, c, x0, v0)
        rep = solve(a, b, c, x0, v0)
        if rep is None:
            return out
        x = rep(ts)
        if not np.all(np.isfinite(x)):
            out["family"] = "nonfinite"
            return out
        _, eh = metrics(x, xr)
        out["err_dop853"] = eh
        if eh > 1e-8 and on_sep:
            # DOP853 leaves the separatrix: 30-digit Taylor reference on t <= 3 (as in the Supplement)
            ts3 = np.linspace(0, min(3.0, ts[-1]), 401)
            xt = taylor_reference(a, b, c, x0, v0, ts3)
            if np.all(np.isfinite(xt)):
                _, eh = metrics(rep(ts3), xt)
                out["err_taylor30"] = eh
                out["ref"] = "taylor30_t3"
            else:
                out["ref"] = "dop853_taylor_nonfinite"
        elif eh > 1e-8:
            # Prefer the higher-precision reference on the same grid.
            xt = taylor_reference(a, b, c, x0, v0, ts)
            if np.all(np.isfinite(xt)):
                _, eh = metrics(x, xt)
                out["err_taylor30"] = eh
                out["ref"] = "taylor30"
            else:
                out["ref"] = "dop853_taylor_nonfinite"
        out.update(represented=True, family=rep.pair.family if rep.pair else rep.kind, err=eh)
    except Exception as exc:                      # report, do not hide
        out["family"] = f"error: {exc}"
    return out


def main():
    results = Path(HERE).parents[1] / "data" / "results"
    outputs = Path(HERE).parent / "outputs"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path, default=results / "v0.2.2_campaign" / "closed_form_campaign_py.csv")
    parser.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1))
    parser.add_argument("--out", type=Path, default=outputs / "campaign_weierstrass.csv")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be at least 1")
    with args.input.open(newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        parser.error("the input contains no campaign cases")
    with mpc.Pool(args.workers) as pool:
        res = pool.map(run, rows, chunksize=8)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(res[0].keys()))
        w.writeheader()
        w.writerows(res)

    print(f"{len(res)} cases from {args.input.name}")
    for cls in ("bounded", "unbounded", "equilibrium", "separatrix"):
        sel = [r for r in res if r["motion"] == cls]
        n_rep = sum(r["represented"] for r in sel)
        n_alg2 = sum(r["alg2_represented"] == "True" for r in sel)
        print(f"  {cls:11s} {len(sel):5d} cases   Algorithm 1 represents {n_rep:5d}   (packages/Algorithm 2: {n_alg2})")
    errs = np.array([r["err"] for r in res if r["represented"]])
    print(f"  hat e_x over represented cases: median {np.median(errs):.1e}  p95 {np.percentile(errs, 95):.1e}"
          f"  max {errs.max():.1e}")
    print(f"  reported against the 30-digit Taylor reference: separatrix t<=3: "
          f"{sum(r['ref'] == 'taylor30_t3' for r in res)}, full window: {sum(r['ref'] == 'taylor30' for r in res)}")
    worst = sorted([r for r in res if r["represented"]], key=lambda r: -r["err"])[:5]
    print("  largest deviations:")
    for r in worst:
        print(f"    {r['group']:5s} {r['motion']:11s} a={r['a']:.6g} b={r['b']:.6g} c={r['c']:.6g} "
              f"x0={r['x0']:.6g} v0={r['v0']:.6g}  fam={r['family']}  hat e_x={r['err']:.2e} ({r['ref']})")
    miss = [r for r in res if not r["represented"]]
    if miss:
        print("  not represented:")
        for r in miss:
            print(f"    {r['group']:5s} {r['motion']:11s} a={r['a']:.6g} b={r['b']:.6g} c={r['c']:.6g} "
                  f"x0={r['x0']:.6g} v0={r['v0']:.6g}  {r['family']}")
    print(f"written: {args.out.name}")
    if any(r["family"].startswith("error:") or r["ref"] == "dop853_taylor_nonfinite" for r in res):
        raise SystemExit("Some cases failed during evaluation or lacked a finite Taylor reference.")


if __name__ == "__main__":
    main()
