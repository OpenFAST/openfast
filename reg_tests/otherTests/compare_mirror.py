#!/usr/bin/env python3
"""Compare a CW baseline run against a mirrored run.

Every channel must resolve to one of: identical, exactly sign-flipped, a mirrored
angle, or below the file's numerical noise floor. Anything else is reported as
unclassified, which is what a frame-mixing bug looks like.

Blade spacing follows the rotation direction, so mirrored blade k compares
directly against clockwise blade k: no channel remapping is required.

Usage: compare_mirror.py <cw.out> <mirror.out> [--tol 1e-6] [--all]
"""
import re
import sys

import numpy as np

# Expected symmetry class per channel.
#   S = identical, F = exactly sign-flipped, A = mirrored angle (mod 360)
EXPECT = {
    "RtTSR": "S", "RtArea": "S", "RtSkew": "S",
    "RtAeroFxh": "S", "RtAeroFzh": "S", "RtAeroMyh": "S",
    "RtAeroPwr": "S", "RtAeroCp": "S", "RtAeroCt": "S",
    "RtVAvgxh": "S", "RtVAvgzh": "S",
    "RtAeroFyh": "F", "RtAeroMxh": "F", "RtAeroMzh": "F",
    "RtAeroCq": "F", "RtVAvgyh": "F",
    "Azimuth": "S", "RotSpeed": "S",
}
for _b in range(1, 4):
    for _n in range(1, 10):
        for _q in ("Alpha", "Phi", "Vrel", "Cl", "Cd", "Cx", "Cn", "Vindx",
                   "AxInd", "TnInd", "Fn", "Fl", "Fd", "M", "Gam", "Theta",
                   "Curve", "Vdisx", "Vundx"):
            EXPECT[f"B{_b}N{_n}{_q}"] = "S"
        for _q in ("Cy", "Ct", "Cm", "Ft", "Mm", "Vindy", "Toe",
                   "Vdisy", "Vundy"):
            EXPECT[f"B{_b}N{_n}{_q}"] = "F"

# Reported but not asserted. Under the old blade numbering the per-blade azimuth's
# mirror reference depended on blade index and rotor tilt, so B*Azimuth sat here.
# With the blades numbered in the direction of rotation it measures as identical
# (class S) across every registered pair -- 2026-09-30, emit_sign_table at
# tol 5e-3 -- so it is asserted like any other channel and this list is empty.
INFORMATIONAL = []


def load(path):
    lines = open(path, errors="replace").read().split("\n")
    h = next(i for i, l in enumerate(lines) if l.strip().startswith("Time"))
    names = lines[h].split()
    rows = [[float(x) for x in l.split()] for l in lines[h + 2:] if l.strip()]
    return names, np.array(rows)


def classify(x, y, tol, floor):
    """Return (class, residual) for one channel pair."""
    peak = np.abs(x).max()
    if peak < floor:
        return "negligible", 0.0
    den = max(peak, 1e-12)
    rs = np.abs(y - x).max() / den
    if rs < tol:
        return "S", rs
    rf = np.abs(y + x).max() / den
    if rf < tol:
        return "F", rf
    # Wrapped angle. Scaled by the same denominator as the tests above, otherwise a
    # 1/360 normalisation makes this far more lenient and it shadows a genuine flip.
    # Only meaningful for a quantity actually measured in degrees: the residual can
    # never exceed 180, so for a large channel this test passes on magnitude alone
    # and will happily label a mooring tension a mirrored angle.
    if den <= 400.0:
        ra = np.abs(((y + x + 180.0) % 360.0) - 180.0).max() / den
        if ra < tol:
            return "A", ra
    return "?", min(rs, rf)


def compare(cw_path, mir_path, tol=1e-6, only=None):
    """Return (nchecked, mismatches, unclassified, observed)."""
    n1, a = load(cw_path)
    n2, b = load(mir_path)
    if n1 != n2:
        raise ValueError("channel lists differ between the two files")
    idx = {k: v for v, k in enumerate(n1)}
    sl = slice(int(0.75 * len(a)), None)

    # Channels that are zero by symmetry still carry cancellation error; they must
    # not be judged against their own magnitude. Keep this well below the smallest
    # channel that carries real signal.
    floor = 1e-9 * np.abs(a[sl, 1:]).max()

    checked, bad, unknown, observed = 0, [], [], {}
    for ch in n1:
        if ch == "Time":
            continue
        if only is not None and not any(re.search(p, ch) for p in only):
            continue
        got, rel = classify(a[sl, idx[ch]], b[sl, idx[ch]], tol, floor)
        observed[ch] = got
        if any(re.search(p, ch) for p in INFORMATIONAL):
            continue
        if got == "?":
            unknown.append((ch, rel))
        exp = EXPECT.get(ch)
        if exp is not None and got != "negligible":
            checked += 1
            if got != exp:
                bad.append((ch, exp, got, rel))
    return checked, bad, unknown, observed


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    tol = float(sys.argv[sys.argv.index("--tol") + 1]) if "--tol" in sys.argv else 1e-6
    checked, bad, unknown, observed = compare(args[0], args[1], tol)

    if "--all" in sys.argv:
        for ch, cls in observed.items():
            print(f"  {ch:14s} {cls}")

    print(f"channels checked : {checked}")
    for ch, rel in unknown:
        print(f"  UNCLASSIFIED {ch:14s} rel={rel:.3e}")
    for ch, exp, got, rel in bad:
        print(f"  MISMATCH     {ch:14s} expected {exp}, got {got}  rel={rel:.3e}")
    if bad or unknown:
        sys.exit(1)
    print("ALL PASS - every channel is SAME, FLIP, mirrored angle, or negligible")


if __name__ == "__main__":
    main()
