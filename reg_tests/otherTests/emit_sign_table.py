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

# label, clockwise case, mirrored case
FILE_PAIRS = [
    ("ElastoDyn + AeroDyn", "5MW_Land_noDLL_Steady_CW",
     "5MW_Land_noDLL_Steady_MirrorRotor"),
    ("BeamDyn blades", "5MW_Land_BD_noDLL_Steady_CW",
     "5MW_Land_BD_noDLL_Steady_MirrorRotor"),
    ("MHK buoyancy", "MHK_RM1_Floating_Steady_CW",
     "MHK_RM1_Floating_Steady_MirrorRotor"),
    ("AeroDisk, yawed", "5MW_Land_ADsk_SED_Yaw_CW",
     "5MW_Land_ADsk_SED_Yaw_MirrorRotor"),
    ("OLAF free wake",
     f"{R}/modules/aerodyn/ad_B1n2_OLAF_CW/ad_driver.outb",
     f"{R}/modules/aerodyn/ad_B1n2_OLAF_MirrorRotor/ad_driver.outb"),
]

# label, single output file, per-rotor channel prefixes
WITHIN_PAIRS = [
    ("AeroDyn nodal outputs",
     os.path.join(REPO_ROOT, "build-docker-double", "reg_tests", "modules",
                   "aerodyn", "ad_MultipleHAWT_MirrorRotor", "ad_driver.T1.outb"),
     os.path.join(REPO_ROOT, "build-docker-double", "reg_tests", "modules",
                   "aerodyn", "ad_MultipleHAWT_MirrorRotor", "ad_driver.T2.outb"),
     None),
    ("Twin-rotor semisubmersible",
     f"{GC}/5MW_MRSemi_DLL_WSt_WavesIrr_MirrorRotor/"
     "5MW_MRSemi_DLL_WSt_WavesIrr_MirrorRotor.outb",
     None, ("R1", "R2")),
]


# Mooring lines and connections exchange in mirror pairs as well, but the pairing
# depends on the layout. Reported separately rather than guessed at.
MOORING = re.compile(r"^(FAIRTEN|ANCHTEN|CON\d|L\d+N|M\d+N|P\d+F)")


def measure(a, na, b, nb, tol, floor):
    out = {}
    idx = {n: i for i, n in enumerate(nb)}
    for i, n in enumerate(na):
        j = idx.get(n)
        if j is None:
            continue
        cls, res = classify(a[:, i], b[:, j], tol, floor)
        out[n] = (cls, res)
    return out



# The behaviour codes, in the order they are presented everywhere else.
BEHAVIOUR = (
    ("S", "identical", "v' = v"),
    ("F", "sign-flipped", "v' = -v"),
    ("A", "mirrored angle", "v' = -v, wrapped"),
    ("negligible", "below the noise floor", "indistinguishable from zero everywhere"),
    ("?", "unresolved", "measured but not classified"),
)


def emit_map(groups, mooring, origin, conflicts, tol, floor):
    """Write the per-channel sign map into docs/source.

    One artefact from one measurement: a YAML file for anything that wants to
    consume the map as data. It is generated -- never hand-edited, because the
    map is measured on every run and a hand edit would be silently overwritten
    or, worse, silently disagree with what the code actually does. The rendered
    per-channel page was retired on 2026-09-30; a grouped summary lives in
    mirror_rotor.rst and the full listing travels with the pull request text.
    """
    docs = os.path.join(REPO_ROOT, "docs", "source", "user", "glue-code")
    labels = {k: lbl for k, lbl, _ in BEHAVIOUR}

    rows = []
    for k, lbl, _ in BEHAVIOUR:
        for n in sorted(groups.get(k, [])):
            rows.append((n, k, lbl, sorted(origin[n])))
    for n in sorted(mooring):
        rows.append((n, "mooring", "set aside, layout-specific",
                     sorted(origin[n])))
    rows.sort(key=lambda r: r[0].upper())

    yml = os.path.join(docs, "mirror_rotor_sign_map.yaml")
    with open(yml, "w") as fh:
        fh.write("# Measured mirror sign map -- GENERATED, do not edit by hand.\n")
        fh.write("# Regenerate with:\n")
        fh.write("#   python3 reg_tests/otherTests/emit_sign_table.py "
                 f"--tol {tol} --emit-map\n")
        fh.write("#\n")
        fh.write("# Every entry is observed by running a registered clockwise/mirrored\n")
        fh.write("# pair and comparing the channel against its counterpart.  Nothing here\n")
        fh.write("# is inferred from reading the source.\n")
        fh.write(f"tolerance: {tol}\n")
        fh.write(f"noise_floor: {floor}\n")
        fh.write("behaviours:\n")
        for k, lbl, meaning in BEHAVIOUR:
            # "?" is a YAML indicator character; quote the key so every parser
            # reads it as the scalar it is meant to be.
            key = f'"{k}"' if k == "?" else k
            fh.write(f'  {key}: {{name: "{lbl}", meaning: "{meaning}"}}\n')
        fh.write('  mooring: {name: "set aside", '
                 'meaning: "pairing is layout-specific"}\n')
        fh.write("channels:\n")
        for n, k, lbl, src in rows:
            fh.write(f'  {n}: {{behaviour: {k}, measured_in: [{", ".join(src)}]}}\n')

    print(f"\nwrote {yml}")
    print(f"  {len(rows)} channels")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tol", type=float, default=2e-3)
    ap.add_argument("--floor", type=float, default=1e-8)
    ap.add_argument("--show-unresolved", action="store_true")
    ap.add_argument("--emit-map", action="store_true",
                    help="write the per-channel sign map into docs/source as YAML "
                         "and as a reStructuredText table")
    a = ap.parse_args()

    seen = defaultdict(set)
    origin = defaultdict(set)

    for label, cw, mir in FILE_PAIRS:
        if cw.startswith("/") or "/" in cw:
            pa, pb = cw, mir          # explicit paths, for the module-level driver cases
        else:
            pa, pb = f"{GC}/{cw}/{cw}.outb", f"{GC}/{mir}/{mir}.outb"
        if not (os.path.exists(pa) and os.path.exists(pb)):
            print(f"  skip {label}: missing baseline", file=sys.stderr)
            continue
        da, ia, _ = pass_fail.readFASTOut(pa)
        db, ib, _ = pass_fail.readFASTOut(pb)
        res = measure(da, list(ia["attribute_names"]), db,
                      list(ib["attribute_names"]), a.tol, a.floor)
        c = Counter(v[0] for v in res.values())
        print(f"{label:28s} same={c['S']:4d} flipped={c['F']:3d} "
              f"wrapped={c['A']:3d} negligible={c['negligible']:3d} "
              f"unresolved={c['?']:3d}")
        for n, (cl, _) in res.items():
            seen[n].add(cl)
            origin[n].add(label)

    for label, pa, pb, prefixes in WITHIN_PAIRS:
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
                j = idx.get(p2 + base)
                if j is None:
                    continue
                cls, r = classify(da[:, na.index(n)], db[:, j], a.tol, a.floor)
                res[base] = (cls, r)
        else:
            res = measure(da, na, db, nbn, a.tol, a.floor)
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

    if a.emit_map:
        emit_map(groups, mooring, origin, conflicts, a.tol, a.floor)


if __name__ == "__main__":
    main()
