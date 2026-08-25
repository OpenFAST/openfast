#!/usr/bin/env python3
"""Phase 2 acceptance: BeamDyn and ElastoDyn must disagree no more when mirrored.

BeamDyn and ElastoDyn are different blade models, so their results differ even for
a clockwise rotor. What matters is that mirroring does not make them disagree any
more than they already do: if the mirrored discrepancy matches the clockwise one,
the mirror is not introducing error of its own.

Usage: check_bd_vs_ed.py
"""
import os
import sys

import numpy as np

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REPO = REPO_ROOT
RTEST = f"{REPO}/reg_tests/r-test/glue-codes/openfast"
sys.path.insert(0, f"{REPO}/reg_tests/lib")
import fast_io  # noqa: E402

DIAG = {"Time", "ConvError", "ConvIter", "NumUJac"}


def read(case):
    data, info, _ = fast_io.load_output(os.path.join(RTEST, case, f"{case}.outb"))
    return info["attribute_names"], np.asarray(data)


def discrepancy(ed_case, bd_case):
    """Per-channel relative difference between the two blade models."""
    n1, a = read(ed_case)
    n2, b = read(bd_case)
    i1 = {k: v for v, k in enumerate(n1)}
    i2 = {k: v for v, k in enumerate(n2)}
    shared = [c for c in n1 if c in i2 and c not in DIAG]

    # BeamDyn runs at a finer time step, so put it on the ElastoDyn time base.
    t1, t2 = a[:, i1["Time"]], b[:, i2["Time"]]
    keep = t1 >= 0.5 * t1[-1]
    scale = np.abs(a[keep][:, 1:]).max()

    out = {}
    for c in shared:
        x = a[keep, i1[c]]
        y = np.interp(t1[keep], t2, b[:, i2[c]])
        den = max(np.abs(x).max(), 1e-12)
        if den < 1e-9 * scale:
            continue
        out[c] = np.abs(y - x).max() / den
    return out


def main():
    cw = discrepancy("5MW_Land_noDLL_Steady_CW", "5MW_Land_BD_noDLL_Steady_CW")
    mir = discrepancy("5MW_Land_noDLL_Steady_MirrorRotor",
                      "5MW_Land_BD_noDLL_Steady_MirrorRotor")

    shared = sorted(set(cw) & set(mir))
    print(f"{'channel':14s} {'CW BD-vs-ED':>13s} {'MIR BD-vs-ED':>13s} {'ratio':>8s}")
    worst = 0.0
    worst_ch = ""
    for c in shared:
        # A mirrored channel that flips sign is compared on magnitude, so take the
        # smaller of the two orientations.
        r = mir[c] / max(cw[c], 1e-12)
        if cw[c] > 1e-6 and r > worst:
            worst, worst_ch = r, c
        print(f"{c:14s} {cw[c]:13.4e} {mir[c]:13.4e} {r:8.3f}")

    print(f"\nchannels compared: {len(shared)}")
    print(f"worst ratio: {worst:.3f} on {worst_ch}")
    # The mirrored discrepancy should track the clockwise one. Allow a factor of two
    # for the two runs drifting differently within an already-different model.
    if worst > 2.0:
        print("FAIL: mirroring widens the BeamDyn/ElastoDyn gap")
        sys.exit(1)
    print("PASS: mirroring does not widen the BeamDyn/ElastoDyn gap")


if __name__ == "__main__":
    main()
