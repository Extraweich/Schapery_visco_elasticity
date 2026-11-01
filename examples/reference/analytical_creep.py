#!/usr/bin/env python3
"""Closed-form Schapery creep solution for a step load (thesis Eq. 52).

    eps(t) = [ g0(s) / E0 + g1g2(s) * sum_i (1 - exp(-t / (a_sigma(s) * tau_i))) / E_i ] * s

Used as reference for examples 01 and 03 (uniaxial creep). The material parameters are
read from the *User Material block of an Abaqus input file, so the reference always uses
exactly the constants of the simulation.

Usage:
    python analytical_creep.py --stress 30
    python analytical_creep.py --stress 30 --inp ../01_creep_3D/creep_3D.inp --csv ref_30MPa.csv
    python analytical_creep.py --stress 30 --abaqus U3_history.csv   # compare with Abaqus output

The optional Abaqus file is a two-column CSV (time, eps33), e.g. exported from the
history output U3 of the top nodes (unit cube: eps33 = U3). Only numpy is required;
matplotlib is used for a plot if it is installed.
"""
import argparse
import math
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def read_props(inp):
    """Return the PROPS of the first *User Material block of an input file."""
    values, active = [], False
    with open(inp) as f:
        for line in f:
            s = line.strip()
            if s.lower().startswith("*user material"):
                n = int(re.search(r"constants\s*=\s*(\d+)", s, re.I).group(1))
                active = True
                continue
            if active:
                if s.startswith("**") or not s:
                    continue
                if s.startswith("*"):
                    break
                values += [float(v) for v in s.split(",") if v.strip()]
                if len(values) >= n:
                    break
    return np.array(values[:n])


def nlfun(props, k, s):
    """Nonlinearizing function k (0: g0, 1: g1*g2, 2: a_sigma) at von Mises stress s."""
    ftype = int(round(props[8 + 5 * k]))
    p1, p2, p3, p4 = props[9 + 5 * k: 13 + 5 * k]
    if ftype == 1:
        return p4 + (p1 - p4) / (1.0 + (s / p3) ** p2)
    if ftype == 2:
        return p1 + p2 * math.exp(-((s - p3) ** 2) / (2.0 * p4 ** 2))
    return 1.0


def creep_strain(props, stress, t):
    E = props[[0, 2, 4, 6]]
    tau = props[[3, 5, 7]]
    g0, g12, a = (nlfun(props, k, stress) for k in range(3))
    t = np.asarray(t, dtype=float)
    transient = sum((1.0 - np.exp(-t / (a * tau[i]))) / E[i + 1] for i in range(3))
    return (g0 / E[0] + g12 * transient) * stress


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stress", type=float, required=True, help="creep stress (MPa)")
    ap.add_argument("--inp", default=os.path.join(HERE, "..", "01_creep_3D", "creep_3D.inp"))
    ap.add_argument("--tend", type=float, default=576000.0, help="end time (s), default 160 h")
    ap.add_argument("--csv", help="write the reference curve to this CSV file")
    ap.add_argument("--abaqus", help="CSV (time, eps) from Abaqus for comparison")
    args = ap.parse_args()

    props = read_props(args.inp)
    t = np.concatenate([[0.0], np.logspace(-2, math.log10(args.tend), 400)])
    eps = creep_strain(props, args.stress, t)
    print(f"PROPS read from {os.path.relpath(args.inp)}")
    for th in (1, 10, 100, 160):
        if th * 3600 <= args.tend:
            print(f"  eps({th:3d} h) = {100 * creep_strain(props, args.stress, th * 3600.0):.4f} %")
    if args.csv:
        np.savetxt(args.csv, np.column_stack([t, eps]), delimiter=",", header="time_s,eps", comments="")
        print(f"reference written to {args.csv}")

    ref = None
    if args.abaqus:
        ref = np.loadtxt(args.abaqus, delimiter=",", ndmin=2)
        load_ramp = 5.0                       # load applied within 5 s in the examples
        mask = ref[:, 0] > 10 * load_ramp
        ana = creep_strain(props, args.stress, ref[mask, 0] - load_ramp / 2)
        dev = np.max(np.abs(ref[mask, 1] - ana) / ana)
        print(f"max. relative deviation Abaqus vs. Eq. 52 for t > {10 * load_ramp:.0f} s: {dev:.2e}")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    plt.semilogx(t[1:] / 3600, 100 * eps[1:], label=f"Eq. 52, {args.stress:g} MPa")
    if ref is not None:
        plt.semilogx(ref[:, 0] / 3600, 100 * ref[:, 1], "o", ms=3, label="Abaqus")
    plt.xlabel("time (h)")
    plt.ylabel("creep strain (%)")
    plt.legend()
    plt.grid(True, which="both", alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"creep_{args.stress:g}MPa.png", dpi=150)
    print(f"plot written to creep_{args.stress:g}MPa.png")


if __name__ == "__main__":
    main()
