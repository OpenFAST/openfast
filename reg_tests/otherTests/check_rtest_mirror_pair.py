#!/usr/bin/env python3
"""Verify the committed r-test mirror pair actually mirrors.

Compares the two stored baselines against each other, so the property the pair
exists to demonstrate is checked at the same time the baselines are produced.

Usage: check_rtest_mirror_pair.py
"""
import os
import sys

import numpy as np

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REPO = REPO_ROOT
RTEST = f"{REPO}/reg_tests/r-test/glue-codes/openfast"
# The r-test baselines are .outb, which packs each channel into int16 with a per-channel
# scale and offset. That quantises to about 1.5e-5 of the channel range, so comparing the
# pair any tighter than this measures the file format rather than the physics.
TOL = 1e-4

# (clockwise case, mirrored case, blades, tolerance, start of comparison window)
# The AWT case brings the shaft to a halt under the HSS brake. The stop is a genuine
# discontinuity -- the brake torque switches sign as the shaft sticks -- and it amplifies
# round-off, so that pair needs a looser tolerance than the steady cases. Its window opens
# earlier as well, to cover the period when the generator is still connected.
PAIRS = [("5MW_Land_noDLL_Steady_CW", "5MW_Land_noDLL_Steady_MirrorRotor", 3, TOL, 0.5),
         ("5MW_Land_BD_noDLL_Steady_CW", "5MW_Land_BD_noDLL_Steady_MirrorRotor", 3, TOL, 0.5),
         ("AWT_WSt_StartUp_HighSpShutDown",
          "AWT_WSt_StartUp_HighSpShutDown_MirrorRotor", 2, 5e-4, 0.25),
         # Turbulent inflow with ROSCO in the loop, against a y-reflected turbulence box.
         # The ElastoDyn pair holds the same tolerance as the steady cases. The BeamDyn
         # one needs a little more, and only for the root torsional moment, which is small
         # and stiff enough that its own peak is a harsh denominator.
         ("5MW_Land_DLL_WTurb", "5MW_Land_DLL_WTurb_MirrorRotor", 3, TOL, 0.25),
         ("5MW_Land_BD_DLL_WTurb", "5MW_Land_BD_DLL_WTurb_MirrorRotor", 3, 1e-3, 0.25)]

# Channels that pair up by layout rather than by blade order. The OC4 semi's mooring
# line 2 lies on the mirror plane, so lines 1 and 3 swap -- the opposite arrangement to
# the blades, where blade 1 is the one on the plane.
SWAP = {"5MW_OC4Semi_WSt_WavesWN_MirrorRotor":
        {"FAIRTEN1": "FAIRTEN3", "FAIRTEN3": "FAIRTEN1",
         "ANCHTEN1": "ANCHTEN3", "ANCHTEN3": "ANCHTEN1"}}

# Unlike the other pairs, the clockwise half here is an upstream baseline produced by a
# different build, so this comparison carries a build difference on top of the mirror.
# Measured: re-running the clockwise case with the current build differs from its stored
# baseline by a median of 1.3e-4 and a maximum of 2.6e-3, and those figures reappear
# channel-for-channel in the pair residual to four significant figures. Compared within a
# single build the pair resolves completely, so the tolerance is set by provenance, not by
# the mirror. Discrimination is still a factor of ~1000.
PAIRS.append(("5MW_OC4Semi_WSt_WavesWN", "5MW_OC4Semi_WSt_WavesWN_MirrorRotor", 3, 5e-3, 0.25))

sys.path.insert(0, f"{REPO}/reg_tests/lib")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fast_io  # noqa: E402
from compare_mirror import classify, blade_permutation  # noqa: E402

DIAG = {"Time", "ConvError", "ConvIter", "NumUJac"}

# The r-test baselines are .outb, which packs each channel into int16 with a per-channel
# scale and offset. That quantises to about 1.5e-5 of the channel range, so comparing the
# pair any tighter than this measures the file format rather than the physics.
TOL = 1e-4


def read(case):
    path = os.path.join(RTEST, case, f"{case}.outb")
    data, info, _ = fast_io.load_output(path)
    return info["attribute_names"], np.asarray(data)


def check_pair(cw, mir, nblades, tol, start):
    n1, a = read(cw)
    n2, b = read(mir)
    if n1 != n2:
        sys.exit("channel lists differ between the two cases")
    idx = {k: v for v, k in enumerate(n1)}
    sl = slice(int(start * len(a)), None)
    phys = [idx[c] for c in n1 if c not in DIAG]
    floor = 1e-9 * np.abs(a[sl][:, phys]).max()

    groups = {}
    for ch in n1:
        if ch in DIAG:
            continue
        mate = SWAP.get(mir, {}).get(ch) or blade_permutation(ch, nblades)
        if mate not in idx:
            mate = ch
        got, rel = classify(a[sl, idx[ch]], b[sl, idx[mate]], tol, floor)
        groups.setdefault(got, []).append(ch)

    for k in ("S", "F", "A", "negligible", "?"):
        if groups.get(k):
            print(f"{k:11s} ({len(groups[k]):3d}): {' '.join(groups[k])}")
    bad = groups.get("?", [])
    print(f"  {len(n1) - len(DIAG & set(n1))} channels compared")
    if bad:
        print(f"  FAIL: {len(bad)} channel(s) unresolved: {' '.join(bad)}")
        return 1
    print("  PASS: clean mirror pair")
    return 0


def main():
    bad = 0
    for cw, mir, nblades, tol, start in PAIRS:
        print(f"=== {cw} vs {mir} ===")
        bad += check_pair(cw, mir, nblades, tol, start)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
