"""Generate the LaTeX tables and number macros of the v0.2.0 manuscript from the result files.

Usage:  python reproduce/make_revision_tables.py
Reads   data/results/{closed_form_campaign_py_summary.json, baseline_summary.json}
        (or the folder in $NOL_RESULTS, e.g. reproduce/outputs after a new run)
        and runs the test suite once (pytest + coverage) to count tests by type.
Writes  reproduce/outputs/tables/numbers_v020.tex, reproduce/outputs/tables/campaign_py.tex, reproduce/outputs/tables/baseline_periods.tex,
        reproduce/outputs/tables/baseline_traj.tex, reproduce/outputs/tables/tests_v020.tex, reproduce/outputs/test_summary.json
"""
import json
import math
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RES = Path(__import__("os").environ.get("NOL_RESULTS", ROOT / "data" / "results"))
TAB = Path(__import__("os").environ.get("NOL_OUTPUTS", ROOT / "reproduce" / "outputs")) / "tables"


def sci(x, digits=1):
    """LaTeX scientific notation, e.g. 3.9\\times10^{-14}."""
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "---"
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    m = x / 10 ** e
    if round(m, digits) >= 10:
        m /= 10; e += 1
    return r"%.*f\times10^{%d}" % (digits, m, e)


def S(x, digits=1):
    return r"\ensuremath{%s}" % sci(x, digits)


def pct(x):
    return r"%.1f\%%" % (100 * x) if x is not None else "---"


def run_tests():
    out = {}
    for mk in ["unit", "integration", "regression", "failure", "slow"]:
        r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", "--co", "-m", mk,
                            str(ROOT / "code_py" / "tests")], capture_output=True, text=True, cwd=ROOT / "code_py")
        m = re.search(r"(\d+)/(\d+) tests collected", r.stdout) or re.search(r"(\d+) tests? collected", r.stdout)
        out[mk] = int(m.group(1))
    cov_json = RES / "coverage.json"
    r = subprocess.run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider",
                        "--cov=nonlinear_oscillator_lab", "--cov-branch", f"--cov-report=json:{cov_json}",
                        str(ROOT / "code_py" / "tests")], capture_output=True, text=True, cwd=ROOT / "code_py")
    m = re.search(r"(\d+) passed", r.stdout); f = re.search(r"(\d+) failed", r.stdout)
    t = re.search(r"in ([\d.]+)s", r.stdout)
    cov = json.load(open(cov_json))
    out.update(passed=int(m.group(1)) if m else 0, failed=int(f.group(1)) if f else 0,
               wall_s=float(t.group(1)) if t else None,
               coverage_total=cov["totals"]["percent_covered"],
               coverage_by_module={Path(k).name: v["summary"]["percent_covered"] for k, v in cov["files"].items()})
    out["total"] = out["passed"] + out["failed"]
    TAB.mkdir(parents=True, exist_ok=True); json.dump(out, open(TAB / "test_summary.json", "w"), indent=1)
    return out


def main():
    camp = json.load(open(RES / "closed_form_campaign_py_summary.json"))
    base = json.load(open(RES / "baseline_summary.json"))
    tests = run_tests()
    N = {}

    # ---------------- campaign table
    desc = {
        "R200": r"published Wolfram campaign (Fig.~\ref{fig:camp}): $a,b\sim U[-12,12]$, $c\sim U[0,1]$, $x_0\sim U[-2,0]$, $v_0\sim U[-2,2]$",
        "RWIDE": r"$a,b\sim U[-12,12]$, $c\sim U[-1,1]$, $x_0,v_0\sim U[-2,2]$ (seed %d)" % camp["_meta"]["seed"],
        "SEP": r"$E=E_s\pm|E_s|10^{-k}$, $k=1,\dots,12$; three multi-well potentials",
        "REP": r"$b^2=4ac$, states at/near minima, $a$, $b$ or $c=0$",
        "UNB": r"$c<0$ or $b<0=c$: $E=E_b(1\pm10^{-k})$, beyond/above barrier",
    }
    rows = []
    for g in ["R200", "RWIDE", "SEP", "REP", "UNB"]:
        s = camp[g]
        bm = s["by_motion"]

        def frac(kind):
            d = bm.get(kind)
            return "%d/%d" % (d["represented"], d["cases"]) if d else "---"
        rows.append(" & ".join([g, desc[g], str(s["cases"]), frac("bounded"), frac("unbounded"),
                                frac("separatrix") if "separatrix" in bm else frac("equilibrium"),
                                "$%s$" % sci(s["err_median"]), "$%s$" % sci(s["err_p95"]),
                                "$%s$" % sci(s["err_max"])]) + r"\\")
        N[g] = s
    (TAB / "campaign_py.tex").write_text("\n".join(rows) + "\n")

    # ---------------- baseline tables
    rows = []
    for meth, s in base["periods"].items():
        mem = "%.0f" % s["peak_mem_kib"] if s["peak_mem_kib"] is not None else "---"
        name = meth.replace("NOL ", r"\texttt{").replace("period_mp(30)", r"period\_mp}, 30 digits") \
                   .replace("period", "period}") if meth.startswith("NOL") else meth
        if meth == "NOL period":
            name = r"\texttt{period} (this work)"
        elif meth.startswith("NOL period_mp"):
            name = r"\texttt{period\_mp}, 30 digits (this work)"
        rows.append(" & ".join([name, "%d/%d" % (s["n"] - s["failed"], s["n"]), "$%s$" % sci(s["err_median"]),
                                "$%s$" % sci(s["err_max"]), pct(s["frac_below_1e10"]),
                                "%.2f" % s["time_median_ms"], mem]) + r"\\")
    (TAB / "baseline_periods.tex").write_text("\n".join(rows) + "\n")
    rows = []
    for meth, s in base["trajectories"].items():
        name = r"\texttt{exact\_solution} (this work)" if meth.startswith("NOL") else "SciPy " + meth
        rows.append(" & ".join([name, "$%s$" % sci(s["err_median"]), "$%s$" % sci(s["err_max"]),
                                "%.1f" % s["time_median_ms"], "%.0f" % s["peak_mem_kib"]]) + r"\\")
    (TAB / "baseline_traj.tex").write_text("\n".join(rows) + "\n")

    # ---------------- tests table
    rows = [
        r"Unit & %d & single functions against analytic, elliptic or arbitrary-precision values\\" % tests["unit"],
        r"Regression & %d & published 40-digit periods; Wolfram representation labels\\" % tests["regression"],
        r"Failure mode & %d & separatrix, degenerate, unbounded and invalid input raise explicit errors\\" % tests["failure"],
        r"Integration & %d & complete workflows, example and benchmark scripts\\" % tests["integration"],
    ]
    (TAB / "tests_v020.tex").write_text("\n".join(rows) + "\n")

    # ---------------- macros
    P = base["periods"]; T = base["trajectories"]; R = camp["R200"]; W = camp["RWIDE"]
    bounded_all = sum(camp[g]["by_motion"].get("bounded", {}).get("cases", 0) for g in camp if not g.startswith("_"))
    bounded_rep = sum(camp[g]["by_motion"].get("bounded", {}).get("represented", 0) for g in camp if not g.startswith("_"))
    total_cases = sum(camp[g]["cases"] for g in camp if not g.startswith("_"))
    macros = {
        "pyTestsTotal": tests["total"], "pyTestsPassed": tests["passed"],
        "pyTestsUnit": tests["unit"], "pyTestsReg": tests["regression"], "pyTestsFail": tests["failure"],
        "pyTestsInt": tests["integration"], "pyCoverage": "%.0f\\%%" % tests["coverage_total"],
        "pyTestTime": "%.0f~s" % tests["wall_s"],
        "campTotal": total_cases, "campBoundedAll": bounded_all, "campBoundedRep": bounded_rep,
        "campLabelsAgree": R["labels_agree_with_wolfram"], "campRN": R["cases"],
        "campRMed": S(R["err_median"]), "campRMax": S(R["err_max"]),
        "campWN": W["cases"], "campWRep": W["represented"], "campWMed": S(W["err_median"]),
        "campWPninetyfive": S(W["err_p95"]), "campWMax": S(W["err_max"]),
        "campBuildMs": "%.0f~ms" % R["build_median_ms"],
        "campMpRechecks": sum(camp[g]["mp25_rechecks"] for g in camp if not g.startswith("_")),
        "benchNOLMed": S(P["NOL period"]["err_median"]), "benchNOLMax": S(P["NOL period"]["err_max"]),
        "benchNOLms": "%.2f~ms" % P["NOL period"]["time_median_ms"],
        "benchMpMax": S(P["NOL period_mp(30)"]["err_max"]),
        "benchMpms": "%.0f~ms" % P["NOL period_mp(30)"]["time_median_ms"],
        "benchQuadMed": S(P["SciPy quad"]["err_median"]), "benchQuadMax": S(P["SciPy quad"]["err_max"]),
        "benchQuadms": "%.2f~ms" % P["SciPy quad"]["time_median_ms"],
        "benchQAWSMed": S(P["SciPy quad, QAWS weight"]["err_median"]),
        "benchQAWSMax": S(P["SciPy quad, QAWS weight"]["err_max"]),
        "benchQAWSms": "%.1f~ms" % P["SciPy quad, QAWS weight"]["time_median_ms"],
        "benchEvMed": S(P["SciPy events"]["err_median"]), "benchEvMax": S(P["SciPy events"]["err_max"]),
        "benchEvms": "%.1f~ms" % P["SciPy events"]["time_median_ms"],
        "benchEvDefMax": S(P["SciPy events (default tol)"]["err_max"]),
        "benchEvDefFailed": P["SciPy events (default tol)"]["failed"],
        "campSepMax": S(camp["SEP"]["err_max"]), "campUnbMax": S(camp["UNB"]["err_max"]),
        "campSepRep": camp["SEP"]["represented"], "campSepN": camp["SEP"]["cases"],
        "campUnbRep": camp["UNB"]["represented"], "campUnbN": camp["UNB"]["cases"],
        "campUnboundedAll": sum(camp[g]["by_motion"].get("unbounded", {}).get("cases", 0) for g in camp if not g.startswith("_")),
        "campUnboundedRep": sum(camp[g]["by_motion"].get("unbounded", {}).get("represented", 0) for g in camp if not g.startswith("_")),
        "campSepLevelAll": sum(camp[g]["by_motion"].get("separatrix", {}).get("cases", 0) for g in camp if not g.startswith("_")),
        "campSepLevelRep": sum(camp[g]["by_motion"].get("separatrix", {}).get("represented", 0) for g in camp if not g.startswith("_")),
        "trajCFMax": S(T["NOL closed form"]["err_max"]), "trajCFms": "%.1f~ms" % T["NOL closed form"]["time_median_ms"],
        "trajDopTenMax": S(T["DOP853 rtol=1e-10"]["err_max"]),
        "trajDopThirteenMax": S(T["DOP853 rtol=1e-13"]["err_max"]),
        "trajDopThirteenms": "%.1f~ms" % T["DOP853 rtol=1e-13"]["time_median_ms"],
        "benchPlatform": base["_meta"]["platform"].replace("_", r"\_"),
        "benchPython": base["_meta"]["python"], "benchNumpy": base["_meta"]["numpy"],
        "benchScipy": base["_meta"]["scipy"], "benchMpmath": base["_meta"]["mpmath"],
    }
    import csv as _csv
    _rows = list(_csv.DictReader(open(RES / "closed_form_campaign_py.csv")))
    _bmax = max(float(r["err"]) for r in _rows if r["ok"] == "True" and r["motion"] == "bounded")
    _allmax = max(float(r["err"]) for r in _rows if r["ok"] == "True")
    macros["campBoundedMax"] = S(_bmax); macros["campAllMax"] = S(_allmax)
    macros["campEqAll"] = sum(camp[g]["by_motion"].get("equilibrium", {}).get("cases", 0) for g in camp if not g.startswith("_"))
    macros["campEqRep"] = sum(camp[g]["by_motion"].get("equilibrium", {}).get("represented", 0) for g in camp if not g.startswith("_"))
    lines = ["%% generated by benchmarks/make_revision_tables.py -- do not edit"]
    for k, v in macros.items():
        lines.append(r"\newcommand{\%s}{%s}" % (k, v))
    (TAB / "numbers_v020.tex").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
