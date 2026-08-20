#!/usr/bin/env python3
"""Measure the observed mirror sign of every output channel across all registered pairs.

The sign table in mirror_rotor.rst claims to be measured rather than asserted, so
it should be generated rather than edited by hand. This reads each clockwise /
mirrored pair, classifies every channel, and reports the aggregate. A channel seen
in more than one pair must agree across them; disagreements are reported rather
than silently resolved.
"""
import argparse
import os
import re
import sys
from collections import Counter, defaultdict

import numpy as np

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "reg_tests", "lib"))
from compare_mirror import classify  # noqa: E402
import pass_fail  # noqa: E402

R = os.path.join(REPO_ROOT, "reg_tests", "r-test")
GC = f"{R}/glue-codes/openfast"

# label, clockwise case, mirrored case, blades
FILE_PAIRS = [
    ("ElastoDyn + AeroDyn", "5MW_Land_noDLL_Steady_CW",
     "5MW_Land_noDLL_Steady_MirrorRotor", 3),
    ("BeamDyn blades", "5MW_Land_BD_noDLL_Steady_CW",
     "5MW_Land_BD_noDLL_Steady_MirrorRotor", 3),
    ("MHK buoyancy", "MHK_RM1_Floating_Steady_CW",
     "MHK_RM1_Floating_Steady_MirrorRotor", 2),
    ("AeroDisk, yawed", "5MW_Land_ADsk_SED_Yaw_CW",
     "5MW_Land_ADsk_SED_Yaw_MirrorRotor", 3),
]

# label, single output file, per-rotor channel prefixes, blades
WITHIN_PAIRS = [
    ("AeroDyn nodal outputs",
     os.path.join(REPO_ROOT, "build-docker-double", "reg_tests", "modules",
                   "aerodyn", "ad_MultipleHAWT_MirrorRotor", "ad_driver.T1.outb"),
     os.path.join(REPO_ROOT, "build-docker-double", "reg_tests", "modules",
                   "aerodyn", "ad_MultipleHAWT_MirrorRotor", "ad_driver.T2.outb"),
     None, 3),
    ("Twin-rotor semisubmersible",
     f"{GC}/5MW_MRSemi_DLL_WSt_WavesIrr_MirrorRotor/"
     "5MW_MRSemi_DLL_WSt_WavesIrr_MirrorRotor.outb",
     None, ("R1", "R2"), 3),
]


def swap_blades(name, nblades):
    """Blade 1 lies on the mirror plane; the others exchange in reverse order."""
    if nblades < 3:
        return name
    order = {"1": "1", "2": "3", "3": "2"}
    # AeroDyn module and nodal blade channels: B2N003Fn, AB2N003Fn
    m = re.match(r"^(A?)B([123])(N\d+.*)$", name)
    if m:
        return f"{m.group(1)}B{order[m.group(2)]}{m.group(3)}"
    # Any other per-blade channel named B<k><word>: B2RootFxr, B2TipTDxr,
    # B2AeroPwr.  Must follow the node pattern above, which is more specific.
    m = re.match(r"^B([123])([A-Za-z].*)$", name)
    if m:
        return f"B{order[m.group(1)]}{m.group(2)}"
    # ElastoDyn trailing-index channels: RootMxc2, OoPDefl3
    m = re.match(r"^(.*?[A-Za-z])([123])$", name)
    if m and not re.search(r"\d$", m.group(1)):
        return m.group(1) + order[m.group(2)]
    return name


# Mooring lines and connections exchange in mirror pairs as well, but the pairing
# depends on the layout. Reported separately rather than guessed at.
MOORING = re.compile(r"^(FAIRTEN|ANCHTEN|CON\d|L\d+N|M\d+N|P\d+F)")


def measure(a, na, b, nb, nblades, tol, floor):
    out = {}
    idx = {n: i for i, n in enumerate(nb)}
    for i, n in enumerate(na):
        j = idx.get(swap_blades(n, nblades), idx.get(n))
        if j is None:
            continue
        cls, res = classify(a[:, i], b[:, j], tol, floor)
        out[n] = (cls, res)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=2e-3)
    ap.add_argument("--floor", type=float, default=1e-8)
    ap.add_argument("--show-unresolved", action="store_true")
    a = ap.parse_args()

    seen = defaultdict(set)
    origin = defaultdict(set)

    for label, cw, mir, nb_ in FILE_PAIRS:
        pa, pb = f"{GC}/{cw}/{cw}.outb", f"{GC}/{mir}/{mir}.outb"
        if not (os.path.exists(pa) and os.path.exists(pb)):
            print(f"  skip {label}: missing baseline", file=sys.stderr)
            continue
        da, ia, _ = pass_fail.readFASTOut(pa)
        db, ib, _ = pass_fail.readFASTOut(pb)
        res = measure(da, list(ia["attribute_names"]), db,
                      list(ib["attribute_names"]), nb_, a.tol, a.floor)
        c = Counter(v[0] for v in res.values())
        print(f"{label:28s} same={c['S']:4d} flipped={c['F']:3d} "
              f"wrapped={c['A']:3d} negligible={c['negligible']:3d} "
              f"unresolved={c['?']:3d}")
        for n, (cl, _) in res.items():
            seen[n].add(cl)
            origin[n].add(label)

    for label, pa, pb, prefixes, nb_ in WITHIN_PAIRS:
        if not os.path.exists(pa):
            print(f"  skip {label}: missing output", file=sys.stderr)
            continue
        da, ia, _ = pass_fail.readFASTOut(pa)
        if pb:
            db, ib, _ = pass_fail.readFASTOut(pb)
            na, nbn = list(ia["attribute_names"]), list(ib["attribute_names"])
        else:
            db, ib = da, ia
            na = nbn = list(ia["attribute_names"])
        if prefixes:
            p1, p2 = prefixes
            na2 = [n for n in na if n.startswith(p1)]
            res = {}
            idx = {n: i for i, n in enumerate(nbn)}
            for n in na2:
                base = n[len(p1):]
                j = idx.get(p2 + swap_blades(base, nb_), idx.get(p2 + base))
                if j is None:
                    continue
                cls, r = classify(da[:, na.index(n)], db[:, j], a.tol, a.floor)
                res[base] = (cls, r)
        else:
            res = measure(da, na, db, nbn, nb_, a.tol, a.floor)
        c = Counter(v[0] for v in res.values())
        print(f"{label:28s} same={c['S']:4d} flipped={c['F']:3d} "
              f"wrapped={c['A']:3d} negligible={c['negligible']:3d} "
              f"unresolved={c['?']:3d}")
        for n, (cl, _) in res.items():
            seen[n].add(cl)
            origin[n].add(label)

    # A channel measured in several pairs must agree; 'negligible' carries no
    # information, so it never contradicts a real observation.
    groups = defaultdict(list)
    conflicts = []
    mooring = []
    for n, classes in sorted(seen.items()):
        if MOORING.match(n):
            mooring.append(n)
            continue
        real = classes - {"negligible", "?"}
        if len(real) > 1:
            conflicts.append((n, sorted(real), sorted(origin[n])))
        elif len(real) == 1:
            groups[real.pop()].append(n)
        elif "?" in classes:
            groups["?"].append(n)
        else:
            groups["negligible"].append(n)

    print("\n=== channel counts by observed behaviour")
    for k, lbl in (("S", "identical"), ("F", "sign-flipped"),
                   ("A", "mirrored angle"), ("negligible", "always ~zero"),
                   ("?", "unresolved")):
        print(f"  {lbl:16s} {len(groups[k])}")
    print(f"  mooring, set aside {len(mooring)}")

    if conflicts:
        print("\n=== CONFLICTS (same channel, different behaviour between pairs)")
        for n, cl, src in conflicts:
            print(f"  {n:20s} {cl}  seen in {src}")

    if a.show_unresolved and groups["?"]:
        print("\n=== unresolved")
        print("  " + ", ".join(groups["?"]))

    for k, lbl in (("S", "IDENTICAL"), ("F", "SIGN-FLIPPED"),
                   ("A", "MIRRORED ANGLE")):
        print(f"\n=== {lbl}")
        print("  " + ", ".join(groups[k]))


if __name__ == "__main__":
    main()
