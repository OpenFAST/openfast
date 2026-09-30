#!/usr/bin/env python3
"""Compare rotor 1 against rotor 2 inside one multi-rotor OpenFAST run.

The glue code prefixes per-rotor channels with R1/R2. When rotor 1 is clockwise
and rotor 2 is mirrored, and the platform, substructure, mooring and sea are all
symmetric about y = 0, the whole system maps onto itself under y -> -y, so R2
must be the mirror image of R1 within the same solve.

Blade spacing follows the rotation sense, so a mirrored blade compares directly
against the clockwise blade of the same number.

Shared (unprefixed) channels belong to the single platform. Under a symmetric
solution the antisymmetric ones must vanish, which is a much sharper statement
than any pairwise comparison, so they are reported separately.
"""
import argparse
import os
import sys
from collections import Counter

import numpy as np

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
from compare_mirror import classify  # noqa: E402

sys.path.insert(0, os.path.join(REPO_ROOT, "reg_tests", "lib"))
import pass_fail  # noqa: E402

# Platform quantities that must be zero in a solution symmetric about y = 0.
ANTISYM_SHARED = ("HydroFyi", "HydroMxi", "HydroMzi",
                  "RBTDYss", "RBRDXss", "RBRDZss")
SYM_SHARED = ("HydroFxi", "HydroFzi", "HydroMyi",
              "RBTDXss", "RBTDZss", "RBRDYss")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outb")
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--floor", type=float, default=1e-8)
    ap.add_argument("--start", type=float, default=0.0)
    ap.add_argument("--show", action="store_true",
                    help="list every channel and its observed class")
    a = ap.parse_args()

    data, info, _ = pass_fail.readFASTOut(a.outb)
    names = list(info["attribute_names"])
    t = data[:, 0]
    k = t >= a.start
    col = {n: data[k, i] for i, n in enumerate(names)}

    r1 = {n[2:]: n for n in names if n.startswith("R1")}
    r2 = {n[2:]: n for n in names if n.startswith("R2")}

    counts, unresolved, rows = Counter(), [], []
    for base in sorted(r1):
        if base not in r2:
            counts["missing"] += 1
            continue
        cls, res = classify(col[r1[base]], col[r2[base]], a.tol, a.floor)
        counts[cls] += 1
        rows.append((base, cls, res))
        if cls == "?":
            unresolved.append((base, res,
                               np.abs(col[r1[base]]).max(),
                               np.abs(col[r2[base]]).max()))

    label = {"S": "same", "F": "flipped", "A": "angle-wrapped",
             "negligible": "negligible", "?": "unresolved"}
    print(f"--- rotor 1 vs rotor 2  ({a.outb})")
    for key in ("S", "F", "A", "negligible", "?", "missing"):
        if counts[key]:
            print(f"{label.get(key, key):16s} {counts[key]}")

    if a.show:
        for base, cls, res in rows:
            print(f"  {base:22s} {label.get(cls, cls):14s} {res:.3e}")

    if unresolved:
        print("\nunresolved:")
        for base, res, p1, p2 in unresolved:
            print(f"  {base:20s} res={res:.3e} "
                  f"peakR1={p1:.6g} peakR2={p2:.6g}")

    print("\n--- shared platform channels")
    print("    (a solution symmetric about y=0 has zero antisymmetric response)")
    for n in ANTISYM_SHARED:
        if n in col:
            print(f"  {n:12s} peak {np.abs(col[n]).max():.6e}   expect ~0")
    for n in SYM_SHARED:
        if n in col:
            print(f"  {n:12s} peak {np.abs(col[n]).max():.6e}   expect nonzero")

    return 1 if unresolved else 0


if __name__ == "__main__":
    sys.exit(main())
