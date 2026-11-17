"""Compare DOP853 with 32-digit Taylor on 30 fixed campaign cases.

The IDs were selected with Random(20261114) among rows that used DOP853 in the
published campaign. Freezing the IDs prevents platform-sensitive threshold
crossings from changing the comparison sample.
"""
import json
import multiprocessing as mpc
import os
import platform
import sys
from pathlib import Path

import mpmath
import numpy as np
import scipy

from closed_form_campaign import exact_solution, mp_reference, read_cases, reference


def job(row):
    a, b, c, x0, v0 = (float(row[k]) for k in ("a", "b", "c", "x0", "v0"))
    ts, dop853, _ = reference(a, b, c, x0, v0)
    x = exact_solution(a, b, c, x0, v0).x(ts)
    taylor = mp_reference(a, b, c, x0, v0, ts, dps=32)
    return dict(case_id=int(row["case_id"]), group=row["group"],
                err_dop853=float(np.max(np.abs(x - dop853))),
                err_mp32bin=float(np.max(np.abs(x - taylor))),
                ref_disagreement=float(np.max(np.abs(dop853 - taylor))))


def main():
    sample = [r for r in read_cases() if r["reference_check"] == "1"]
    if len(sample) != 30:
        raise ValueError("Expected 30 fixed reference-check cases")
    with mpc.Pool(2) as pool:
        rows = pool.map(job, sample, chunksize=1)
    ratios = [r["err_dop853"] / r["err_mp32bin"] for r in rows if r["err_mp32bin"] > 0]
    out = dict(n=len(rows), ratio_n=len(ratios), ratio_median=float(np.median(ratios)),
               ratio_min=min(ratios), ratio_max=max(ratios),
               err_dop853_median=float(np.median([r["err_dop853"] for r in rows])),
               err_mp32bin_median=float(np.median([r["err_mp32bin"] for r in rows])),
               max_ref_disagreement=max(r["ref_disagreement"] for r in rows), rows=rows,
               _meta=dict(seed=20261114, python=sys.version.split()[0], platform=platform.platform(),
                          scipy=scipy.__version__, numpy=np.__version__, mpmath=mpmath.__version__))
    folder = Path(os.environ.get("NOL_OUTPUTS", Path(__file__).parent / "outputs"))
    folder.mkdir(parents=True, exist_ok=True)
    with (folder / "dop853_reference_check.json").open("w") as f:
        json.dump(out, f, indent=2)
        f.write("\n")
    print({k: v for k, v in out.items() if k != "rows"})


if __name__ == "__main__":
    main()
