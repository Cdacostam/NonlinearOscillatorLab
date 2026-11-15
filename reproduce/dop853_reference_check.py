"""How much of e_x <= 1e-10 against DOP853 is reference error? (v7, observation O1)

Draws 30 represented campaign rows with ref = dop853 (random.Random(20261114), fixed before any run) and
recomputes e_x against a 32-digit mpmath reference started from the exact binary inputs.
Usage: python reproduce/dop853_reference_check.py
"""
import csv, json, random, sys, multiprocessing as mpc
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rereference_campaign as rr  # noqa: E402
from nonlinear_oscillator_lab.closed_form import exact_solution  # noqa: E402


def job(r):
    a, b, c, x0, v0 = (float(r[k]) for k in ("a", "b", "c", "x0", "v0"))
    ts, _, _ = rr.cfc.reference(a, b, c, x0, v0)
    xc = exact_solution(a, b, c, x0, v0).x(ts)
    e = float(np.max(np.abs(xc - rr.mp_reference(a, b, c, x0, v0, ts, dps=32, exact=True))))
    return dict(group=r["group"], motion=r["motion"], err_dop853=float(r["err"]), err_mp32bin=e)


if __name__ == "__main__":
    rows = [r for r in csv.DictReader(open(HERE.parent / "data" / "results" / "closed_form_campaign_py_exactref.csv")) if r["ref"] == "dop853"]
    sample = random.Random(20261114).sample(rows, 30)
    with mpc.Pool(2) as p:
        res = p.map(job, sample, chunksize=1)
    ratio = np.array([q["err_dop853"] / q["err_mp32bin"] for q in res if q["err_mp32bin"] > 0])
    out = dict(n=len(res), ratio_median=float(np.median(ratio)), ratio_min=float(ratio.min()), ratio_max=float(ratio.max()),
               err_dop853_median=float(np.median([q["err_dop853"] for q in res])),
               err_mp32bin_median=float(np.median([q["err_mp32bin"] for q in res])), rows=res)
    o = Path(__import__("os").environ.get("NOL_OUTPUTS", HERE / "outputs")); o.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(o / "dop853_reference_check.json", "w"), indent=1)
    print({k: v for k, v in out.items() if k != "rows"})
