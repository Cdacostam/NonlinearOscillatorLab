"""Figure of the period accuracy (Figure 3 of the paper, Supplementary Figure S5).

Panel (a): relative error of the machine-precision period for 8 families x 40 energies,
Python `period` and Wolfram `NOLPeriodFromEnergy`, both against 40-digit references.
The Python series is recomputed on the energies recorded in
data/wolfram/period_accuracy.csv, which supplies the references and the archived
Wolfram series. Errors below 1e-17 are shown at 1e-17.

Panel (b): separatrix asymptotics of the double well x'' - x + x^3 = 0 (E_s = 0).
For dE = 10^-k, k = 1, 1.5, ..., 16, the script computes with `period_mp` at 40 digits
    inner orbit (x0 = 1, E = -dE):  T - ln(1/dE)    -> ln 16
    outer orbit (x0 = 0, E = +dE):  T - 2 ln(1/dE)  -> 2 ln 16
The initial velocity is passed as a 60-digit decimal string, so that the energy is
-dE or +dE to about 60 digits.

Run from the archive root after `python -m pip install ./code_py matplotlib`:
    python reproduce/period_accuracy_figure.py            # computes panel (b), a few seconds
    python reproduce/period_accuracy_figure.py --from-csv data/results/period_accuracy_panel_b.csv
Outputs (reproduce/outputs/ or $NOL_OUTPUTS): period_accuracy_panel_a.csv, period_accuracy_panel_b.csv,
fig6_period_accuracy.pdf and fig6_period_accuracy.png. Archived copies:
data/results/period_accuracy_panel_b.csv and data/figures/fig6_period_accuracy.{pdf,png}.
"""
import argparse
import csv
import os
from fractions import Fraction
from pathlib import Path

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
FAMILIES = ["hardening", "softening", "dw_inner", "dw_outer",
            "cq_demo", "tw_center", "tw_side", "tw_outer"]
NAMES = {"hardening": "hardening", "softening": "softening",
         "dw_inner": "double well, inner", "dw_outer": "double well, outer",
         "cq_demo": "a=1, b=2/5, c=1/5", "tw_center": "triple well, centre",
         "tw_side": "triple well, side", "tw_outer": "triple well, outer"}
BLUE, ORANGE, GREEN, PURPLE = "#1f5fa8", "#d1661a", "#2a8c55", "#7a4fa3"
KS = [Fraction(j, 2) for j in range(2, 33)]          # k = 1, 1.5, ..., 16
DPS = 40


def as_float(s):
    """Values of period_accuracy.csv: decimals, or rationals such as 1/100000000000000000."""
    s = s.strip().strip('"')
    return float(Fraction(s)) if "/" in s else float(s)


def panel_a_data(path):
    from baseline_comparison import FAMILIES as orbit_specs
    from nonlinear_oscillator_lab import Model, period

    py = {f: [] for f in FAMILIES}
    wl = {f: [] for f in FAMILIES}
    rows = []
    ctx = mp.mp.clone()
    ctx.dps = 50
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            family = row["family"]
            a, b, c, x0, _ = orbit_specs[family]
            model = Model(a, b, c)
            energy = as_float(row["E"])
            v0 = float(np.sqrt(2 * (energy - model.potential(x0))))
            ref = ctx.mpf(row["T_reference"].split(chr(96))[0].replace("*^", "e"))
            t = period(model, x0, v0)
            py_err = float(abs(ctx.mpf(t) - ref)/ref)
            wl_err = as_float(row["relerr_wolfram_factored"])
            py[family].append(py_err)
            wl[family].append(wl_err)
            rows.append(dict(family=family, E=energy, x0=x0, v0=v0, period_python=t,
                             period_reference=ctx.nstr(ref, 42), error_python=py_err,
                             error_wolfram_archived=wl_err))
    if any(len(py[f]) != 40 for f in FAMILIES):
        raise ValueError("Expected 40 records in each period family")
    return py, wl, rows


def panel_b_data():
    from nonlinear_oscillator_lab import Model, period_mp
    m = Model(-1.0, 1.0, 0.0)
    ctx = mp.mp.clone()
    ctx.dps = 80
    rows = []
    for k in KS:
        dE = ctx.mpf(10) ** (-ctx.mpf(k.numerator) / k.denominator)
        v_in = ctx.sqrt(2 * (ctx.mpf(1) / 4 - dE))   # V(1) = -1/4, so E = -dE
        v_out = ctx.sqrt(2 * dE)                      # V(0) = 0, so E = +dE
        t_in = period_mp(m, 1, ctx.nstr(v_in, 60), dps=DPS, tol=1e-30)
        t_out = period_mp(m, 0, ctx.nstr(v_out, 60), dps=DPS, tol=1e-30)
        ln = ctx.log(1 / dE)
        rows.append({"k": str(float(k)), "dE": ctx.nstr(dE, 20),
                     "T_inner": mp.nstr(t_in, DPS), "T_outer": mp.nstr(t_out, DPS),
                     "inner_T_minus_ln": ctx.nstr(t_in - ln, 30),
                     "outer_T_minus_2ln": ctx.nstr(t_out - 2 * ln, 30)})
    return rows


def write_csv(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def draw(py, wl, rows, stem):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"],
                         "mathtext.fontset": "stix", "font.size": 10, "axes.linewidth": 0.7,
                         "pdf.fonttype": 42, "svg.hashsalt": "nol"})
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw={"wspace": 0.28})

    rng = np.random.default_rng(1)
    for off, data, col, lab in ((-0.14, py, ORANGE, "Python"),
                                (0.14, wl, BLUE, "Wolfram Language")):
        xs, ys = [], []
        for k, f in enumerate(FAMILIES, start=1):
            v = np.maximum(np.asarray(data[f]), 1e-17)
            xs.extend(k + off + rng.uniform(-0.08, 0.08, v.size))
            ys.extend(np.log10(v))
        ax.scatter(xs, ys, s=4, color=col, linewidths=0, label=lab)
    ax.set_xlim(0.5, 8.5)
    ax.set_ylim(-17.3, -4)          # values floored at 1e-17 stay visible
    ax.set_xticks(range(1, 9))
    ax.set_xticklabels([NAMES[f] for f in FAMILIES], rotation=30, ha="right", fontsize=8.5)
    ax.set_ylabel(r"$\log_{10}$ relative error of $T$")
    ax.set_title("(a) relative period error (8 families × 40 energies)", loc="left",
                 fontsize=10, fontweight="bold")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.36), ncol=2, frameon=False,
              fontsize=9, markerscale=2.5, handletextpad=0.2)

    k = np.array([float(r["k"]) for r in rows])
    t_in = np.array([float(r["inner_T_minus_ln"]) for r in rows])
    t_out = np.array([float(r["outer_T_minus_2ln"]) for r in rows])
    bx.plot(-k, t_in, "-o", color=GREEN, ms=3.5, lw=1.0)
    bx.plot(-k, t_out, "-s", color=PURPLE, ms=3.5, lw=1.0)
    for y in (np.log(16), 2 * np.log(16)):
        bx.axhline(y, color="gray", ls="--", lw=0.8)
    bx.set_xlim(-16.5, -0.5)
    bx.set_xticks([-15, -10, -5])
    bx.set_ylim(2, 6.2)
    bx.set_xlabel(r"$\log_{10}|E-E_s|$")
    bx.set_ylabel(r"$T - k\,\ln(1/|E-E_s|)$")
    bx.text(0.55, 0.30, r"inner orbit, $k=1$ $\to$ ln 16", transform=bx.transAxes,
            color=GREEN, ha="center", fontsize=9.5)
    bx.text(0.55, 0.72, r"outer orbit, $k=2$ $\to$ 2 ln 16", transform=bx.transAxes,
            color=PURPLE, ha="center", fontsize=9.5)
    bx.set_title("(b) separatrix asymptotics, double well", loc="left",
                 fontsize=10, fontweight="bold")
    for a in (ax, bx):
        a.tick_params(direction="in", top=True, right=True, width=0.7, labelsize=9)

    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", metadata={"CreationDate": None})
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight", dpi=300, metadata={"Software": None})
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--accuracy-csv", default=ROOT / "data" / "wolfram" / "period_accuracy.csv")
    p.add_argument("--from-csv", help="read panel (b) from this CSV instead of computing it")
    args = p.parse_args()
    out = Path(os.environ.get("NOL_OUTPUTS", ROOT / "reproduce" / "outputs"))
    out.mkdir(parents=True, exist_ok=True)

    py, wl, panel_a = panel_a_data(args.accuracy_csv)
    write_csv(panel_a, out / "period_accuracy_panel_a.csv")
    if args.from_csv:
        rows = read_csv(args.from_csv)
    else:
        rows = panel_b_data()
        write_csv(rows, out / "period_accuracy_panel_b.csv")
        print("wrote", out / "period_accuracy_panel_b.csv")
    draw(py, wl, rows, out / "fig6_period_accuracy")
    allpy = np.concatenate([py[f] for f in FAMILIES])
    allwl = np.concatenate([np.maximum(wl[f], 1e-17) for f in FAMILIES])
    print(f"panel (a): Python median {np.median(allpy):.1e}, max {allpy.max():.1e}; "
          f"Wolfram median {np.median(allwl):.1e}, max {allwl.max():.1e}")
    print(f"panel (b): k=16 inner {rows[-1]['inner_T_minus_ln'][:14]} (ln 16 = {np.log(16):.12f}), "
          f"outer {rows[-1]['outer_T_minus_2ln'][:14]} (2 ln 16 = {2 * np.log(16):.12f})")
    print("wrote", out / "fig6_period_accuracy.pdf", "and .png")


if __name__ == "__main__":
    main()
