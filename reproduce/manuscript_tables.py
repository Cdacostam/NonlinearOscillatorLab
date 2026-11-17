"""Generate the revised manuscript's shared counters and table rows from saved runs.

From the bundle root: python reproduce/manuscript_tables.py ../source
"""
import argparse
import csv
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from fractions import Fraction

import numpy as np


DATA = Path(__file__).resolve().parents[1] / "data" / "results"
PERIOD_NAMES = {
    "NOL period": r"\texttt{period}",
    "NOL period_mp(30)": r"\texttt{period\_mp}, 30 digits",
    "SciPy quad": r"SciPy \texttt{quad}",
    "SciPy QAWS, direct": "SciPy QAWS, direct",
    "SciPy QAWS, divided differences": "SciPy QAWS DD",
    "SciPy events": "SciPy events",
    "SciPy events (default tol)": "SciPy events, default",
}
DESIGNS = {
    "R200": r"$a,b\sim U[-12,12]$, $c\sim U[0,1]$, $x_0\sim U[-2,0]$, $v_0\sim U[-2,2]$",
    "RWIDE": r"$a,b\sim U[-12,12]$, $c\sim U[-1,1]$, $x_0,v_0\sim U[-2,2]$",
    "SEP": r"Near separatrices of three potentials",
    "REP": r"$b^2=4ac$, near minima, vanishing coefficients",
    "UNB": r"Escape thresholds and exterior states",
}


def read_json(name):
    return json.loads((DATA / name).read_text())


def sci(v):
    if v == 0:
        return "$0$"
    number, exponent = f"{v:.1e}".split("e")
    return rf"\ensuremath{{{number}\times10^{{{int(exponent)}}}}}"


def write_rows(out, filename, rows):
    (out / filename).write_text("".join(" & ".join(map(str, row)) + r"\\" + "\n" for row in rows))


def main(out):
    period = read_json("v0.2.2_period/baseline_summary.json")
    periods = period["periods"]
    write_rows(out, "baseline_rows.tex", [
        (label, f"{d['n']-d['failed']}/{d['n']}", sci(d["err_median"]), sci(d["err_max"]),
         rf"{100*d['frac_below_1e10']:.1f}\%", f"{d['time_median_ms']:.2f}",
         f"{d['peak_mem_kib']:.0f}" if d["peak_mem_kib"] is not None else "---")
        for key, label in PERIOD_NAMES.items() for d in (periods[key],)
    ])
    campaign = read_json("v0.2.2_campaign/closed_form_campaign_py_summary.json")
    campaign_rows = []
    for g, design in DESIGNS.items():
        d = campaign[g]
        motion = d["by_motion"]
        counts = [
            f"{motion[k]['represented']}/{motion[k]['cases']}" if k in motion else "---"
            for k in ("bounded", "unbounded", "equilibrium" if g == "REP" else "separatrix")
        ]
        campaign_rows.append([g, design, d["cases"], *counts,
                              *(sci(d[m]) for m in ("err_median", "err_p95", "err_max"))])
    write_rows(out, "campaign_rows.tex", campaign_rows)

    period_folder = DATA / "v0.2.2_period"
    with (period_folder / "python_periods.csv").open(newline="") as f:
        examples = {r["id"]: r for r in csv.DictReader(f)}
    with (DATA.parent / "wolfram" / "benchmark_periods.csv").open(newline="") as f:
        wolfram_examples = list(csv.DictReader(f))
    assert len(examples) == len(wolfram_examples) == 9

    def tex_input(text):
        value = Fraction(text)
        if value.denominator == 1:
            return "$" + str(value.numerator) + "$"
        sign = "-" if value < 0 else ""
        return "$" + rf"{sign}\frac{{{abs(value.numerator)}}}{{{value.denominator}}}" + "$"

    write_rows(out, "period_examples_rows.tex", [
        (r["case"].replace("_", r"\_"), *(tex_input(r[k]) for k in ("a", "b", "c", "x0", "v0")),
         r["T_40digit"], sci(float(examples[r["case"]]["relative_difference_binary"])),
         sci(float(r["relerr_wolfram"])), sci(float(r["closed_form_maxerr_0_20"])))
        for r in wolfram_examples
    ])
    examples_max = max(float(r["relative_difference_binary"]) for r in examples.values())
    with (period_folder / "period_accuracy_panel_a.csv").open(newline="") as f:
        panel_a = list(csv.DictReader(f))
    py_errors = [float(r["error_python"]) for r in panel_a]
    assert len(py_errors) == 320

    tests = DATA / "v0.2.2_tests"
    suite = ET.parse(tests / "pytest.xml").getroot().find("testsuite").attrib
    assert all(int(suite[x]) == 0 for x in ("failures", "errors", "skipped"))
    classes = dict(re.findall(r"(\w+)\n(\d+)/\d+ tests collected", (tests / "categories.log").read_text()))
    assert sum(map(int, classes.values())) == int(suite["tests"])
    cov = read_json("v0.2.2_tests/coverage.json")["totals"]
    folder = DATA / "v0.2.2_campaign"
    with (folder / "closed_form_campaign_py.csv").open(newline="") as f:
        rows = list(csv.DictReader(f))
    ref = read_json("v0.2.2_campaign/dop853_reference_check.json")
    accepted = [r for r in rows if r["ok"] == "True"]
    rejected = [r for r in rows if r["ok"] == "False"]
    assert len(rows) == 1436 and len(accepted) + len(rejected) == len(rows)
    assert all(r["ref"] in {"dop853", "mp25"} and np.isfinite(float(r["err"])) for r in accepted)
    assert all(r["motion"] == "unbounded" and r["reason"] == "No candidate representation passed the acceptance checks"
               for r in rejected)
    classes_summary = {}
    for name in ("bounded", "separatrix", "equilibrium", "unbounded"):
        group = [r for r in rows if r["motion"] == name]
        errors = [float(r["err"]) for r in group if r["ok"] == "True"]
        classes_summary[name] = dict(total=len(group), represented=len(errors),
                                     rejected=len(group)-len(errors), error_median=float(np.median(errors)),
                                     error_p95=float(np.percentile(errors, 95)), error_max=max(errors))
    (out / "campaign_by_motion.json").write_text(json.dumps(classes_summary, indent=2) + "\n")

    cf = period["trajectories"]["NOL closed form"]
    ode = period["trajectories"]["DOP853 rtol=1e-13"]
    stats = {
        "PythonTests": suite["tests"], "UnitTests": classes["unit"],
        "FailureTests": classes["failure"], "RegressionTests": classes["regression"],
        "IntegrationTests": classes["integration"],
        "ExamplesPeriodMax": sci(examples_max),
        "FigurePeriodMedian": sci(float(np.median(py_errors))),
        "FigurePeriodMax": sci(max(py_errors)),
        "CombinedCoverage": rf"{cov['percent_covered']:.2f}\%",
        "BranchCoverage": rf"{cov['percent_branches_covered']:.2f}\%",
        "BranchCount": f"{cov['covered_branches']}/{cov['num_branches']}",
        "PeriodMax": sci(periods["NOL period"]["err_max"]),
        "QawsMax": sci(periods["SciPy QAWS, divided differences"]["err_max"]),
        "PeriodTime": f"{periods['NOL period']['time_median_ms']:.2f}",
        "QawsTime": f"{periods['SciPy QAWS, divided differences']['time_median_ms']:.2f}",
        "QawsRawFailed": periods["SciPy QAWS, direct"]["failed"],
        "DefaultFailed": periods["SciPy events (default tol)"]["failed"],
        "TrajectoryCFMax": sci(cf["err_max"]), "TrajectoryCFTime": f"{cf['time_median_ms']:.2f}",
        "TrajectoryODEMax": sci(ode["err_max"]), "TrajectoryODETime": f"{ode['time_median_ms']:.2f}",
        "BoundedTotal": classes_summary["bounded"]["total"], "CriticalTotal": classes_summary["separatrix"]["total"],
        "EquilibriumTotal": classes_summary["equilibrium"]["total"], "UnboundedTotal": classes_summary["unbounded"]["total"],
        "UnboundedRepresented": classes_summary["unbounded"]["represented"],
        "RejectedTotal": len(rejected), "RejectedRWIDE": sum(r["group"] == "RWIDE" for r in rejected),
        "RejectedUNB": sum(r["group"] == "UNB" for r in rejected),
        "ExteriorID": next(r["case_id"] for r in rejected if r["a"] == "9.269312955746553"),
        "LabelsAgree": campaign["R200"]["labels_agree_with_wolfram"],
        "CampaignMax": sci(max(float(r["err"]) for r in accepted)),
        "TaylorCount": sum(r["ref"] == "mp25" for r in accepted),
        "RefDOPMedian": sci(ref["err_dop853_median"]), "RefMPMedian": sci(ref["err_mp32bin_median"]),
        "RefRatio": f"{ref['ratio_median']:.1f}",
    }
    (out / "revision_numbers.tex").write_text(
        "% Generated from the bundled runs by reproduce/manuscript_tables.py.\n"
        + "".join(rf"\newcommand{{\{k}}}{{{v}}}" + "\n" for k, v in stats.items())
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tex_dir", type=Path)
    main(parser.parse_args().tex_dir)
