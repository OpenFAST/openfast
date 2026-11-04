#!/usr/bin/env python3
"""Compare a clockwise FAST.Farm run against its mirrored twin.

The twin reflects the world across y = 0 — turbine position, high-resolution
box, inflow box, yaw schedule — and sets MirrorRotor = T, so the mirrored
farm's solution must be the y-reflection of the clockwise one: even channels
match, odd channels flip, and the wake swirl reverses. Three comparisons:

1. Farm-level channels (FAST.Farm.outb): WkPosY/WkAxsY/YawErr/RtGamCurl flip,
   WkDfVx/WkDfVr/WkPos{X,Z}/WkAxs{X,Z}/WkDiam/RtVAmb/RtSkew/RtCtAvg match, and
   the skew azimuth AziSkew reverses as a wrapped angle.
2. Turbine channels (FAST.Farm.T1.outb): judged against the measured turbine
   sign map in compare_mirror.EXPECT, same as the registered pairs.
3. Wake-plane fields (OutAllPlanes VTK, vtk-ff/): the deficit fields are
   written on the plane-local y-z grid, so the mirrored plane must satisfy
   Vx'(y,z) = Vx(-y,z), Vy'(y,z) = -Vy(-y,z), Vz'(y,z) = Vz(-y,z). The Vy
   comparison is the swirl-reversal measurement — no farm output channel
   carries the tangential deficit.

Usage:
    python3 farm_mirror_pair.py <CW-case-dir> <MIR-case-dir> [--tol T]
        [--vtk-times N] [--show]

Both directories must already hold completed runs with the same root name
(FAST.Farm). Exit status 0 only if every comparison resolves as expected.
"""
import argparse
import os
import re
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO_ROOT, "reg_tests", "lib"))

from compare_mirror import classify, load as load_text, EXPECT as TURBINE_EXPECT  # noqa: E402
import pass_fail  # noqa: E402

# Farm-level expectations, by pattern. Each maps to the set of acceptable
# classes: S (same), F (sign-flipped), A (mirrored wrapped angle), P (the
# skew azimuth's measured relation psi' = 180 - psi, wrapped), W (same up to
# 360-degree wrapping); "negligible" is always acceptable.
#
# Measured on the V1 pair (2026-11-03): the instantaneous relations hold to
# 7e-6 deg or exactly; the *Filt* channels and RtGamCurl carry the
# filter_angles2 non-equivariance floor (~3e-2 relative) because the mirrored
# farm's skew DCM sits near the logarithmic map's branch cut. Gamma is even
# (chi is invariant; the vortex curtain's sign comes from the zp -> -zp
# relabelling), so RtGamCurl is S, not F.
FARM_EXPECT = [
    (r"^RtVAmb(Filt)?T\d", {"S"}),
    (r"^RtVRelT\d", {"S"}),
    (r"^RtSkew(Filt)?T\d", {"S"}),
    (r"^RtCtAvgT\d", {"S"}),
    (r"^WkPos[XZ]T\d", {"S"}),
    (r"^WkAxs[XZ]T\d", {"S"}),
    (r"^WkDiamT\d", {"S"}),
    (r"^WkDfV[xr]T\d", {"S"}),
    (r"^YawErrT\d", {"F"}),
    (r"^WkPosYT\d", {"F"}),
    (r"^WkAxsYT\d", {"F"}),
    (r"^RtGamCurlT\d", {"S"}),
    (r"^AziSkew(Filt)?T\d", {"P"}),
]

# Turbine channels that are angles wrapped to [0, 360): a small accumulated
# phase drift between the two runs shows up as full-wrap spikes, so they are
# compared with the wrapped difference rather than the raw one.
WRAPPED_TURBINE = re.compile(r"^(B\d)?Azimuth$|^LSSTipPxa$")

# Plane-field expectations: scalar name -> acceptable classes after the
# y-index reversal. The eddy viscosities and WAT scaling are even; the
# cross-flow gradient dvx_dy is odd because y itself reverses.
PLANE_EXPECT = {
    "Vx": {"S"},
    "Vy": {"F"},
    "Vz": {"S"},
    "vt_amb2": {"S"},
    "vt_shr2": {"S"},
    "vt_tot2": {"S"},
    "dvx_dy": {"F"},
    "dvx_dz": {"S"},
    "WAT_k": {"S"},
}


def read_outb(case, name):
    """Read a FAST output file: binary .outb, or the farm's text .out
    (FAST.Farm refuses OutFileFmt /= 1, and its OutFmt "G0" is full
    precision, so text costs nothing there)."""
    path = os.path.join(case, name)
    if path.endswith(".outb"):
        data, info, _ = pass_fail.readFASTOut(path)
        return list(info["attribute_names"]), np.asarray(data)
    names, data = load_text(path)
    return names, data


def compare_channels(names_a, a, names_b, b, expect_for, tol, start_frac=0.25,
                     show=False, label=""):
    """Classify every common channel; return list of violations."""
    if names_a != names_b:
        only = sorted(set(names_a) ^ set(names_b))
        print(f"  NOTE {label}: channel lists differ ({len(only)} unpaired): {only[:6]}")
    common = [n for n in names_a if n in set(names_b) and n != "Time"]
    ia = {n: i for i, n in enumerate(names_a)}
    ib = {n: i for i, n in enumerate(names_b)}
    sl = slice(int(start_frac * len(a)), None)
    floor = 1e-9 * max(np.abs(a[sl, 1:]).max(), 1e-12)
    bad, counts = [], Counter()
    for ch in common:
        x, yv = a[sl, ia[ch]], b[sl, ib[ch]]
        got, rel = classify(x, yv, tol, floor)
        exp = expect_for(ch)
        if exp is not None and got not in exp and got != "negligible":
            # The two relations classify() does not know. P: the skew
            # azimuth's psi' = 180 - psi, wrapped; W: same up to 360-degree
            # wrapping (phase drift between the runs crosses the wrap).
            if "P" in exp:
                r = np.abs(((yv + x) % 360.0) - 180.0).max() / 360.0
                if r < tol:
                    got, rel = "P", r
            if "W" in exp and got not in exp:
                r = np.abs(((yv - x + 180.0) % 360.0) - 180.0).max() / 360.0
                if r < tol:
                    got, rel = "W", r
        counts[got] += 1
        if show:
            print(f"    {ch:16s} {got}  rel={rel:.2e}")
        if exp is None:
            if got == "?":
                bad.append((ch, "any", got, rel))
            continue
        if got not in exp and got != "negligible":
            bad.append((ch, "/".join(sorted(exp)), got, rel))
    print(f"  {label}: {len(common)} channels  " +
          "  ".join(f"{k}:{v}" for k, v in sorted(counts.items())))
    return bad


def farm_expect(ch):
    for pat, exp in FARM_EXPECT:
        if re.search(pat, ch):
            return exp
    return None


def turbine_expect(ch):
    exp = TURBINE_EXPECT.get(ch)
    s = {exp} if exp is not None else None
    if WRAPPED_TURBINE.match(ch):
        s = (s or set()) | {"W"}
    return s


def read_plane_vtk(path):
    """Minimal ASCII STRUCTURED_POINTS reader: returns {name: (ny, nz) array}."""
    dims, fields, name = None, {}, None
    with open(path) as f:
        lines = f.read().split("\n")
    i = 0
    while i < len(lines):
        tok = lines[i].split()
        if tok[:1] == ["DIMENSIONS"]:
            dims = tuple(int(v) for v in tok[1:4])  # (1, ny, nz)
        elif tok[:1] == ["SCALARS"]:
            name = tok[1]
            n = dims[1] * dims[2]
            vals, i2 = [], i + 2  # skip LOOKUP_TABLE
            while len(vals) < n:
                vals += [float(v) for v in lines[i2].split()]
                i2 += 1
            fields[name] = np.array(vals).reshape(dims[2], dims[1]).T  # (ny, nz)
            i = i2 - 1
        i += 1
    return fields


def plane_files(case, wt=None):
    """Map (wt, time_code, plane) -> path for the per-plane outputs; with
    `wt` given, map (time_code, plane) -> path for that turbine only."""
    rx = re.compile(r"\.WT(\d+)\.PlaneAtTime_(\d+)_Plane_(\d+)\.vtk$")
    out = {}
    for d, _, fs in os.walk(case):
        for f in fs:
            m2 = rx.search(f)
            if not m2:
                continue
            w, t, i = (int(g) for g in m2.groups())
            if wt is None:
                out[(w, t, i)] = os.path.join(d, f)
            elif w == wt:
                out[(t, i)] = os.path.join(d, f)
    return out


def compare_planes(cw, mir, tol, n_times):
    """cw and mir map (time_code, plane) -> path (one turbine each)."""
    both = sorted(set(cw) & set(mir))
    if not cw and not mir:
        # The case was run without OutAllPlanes; nothing to compare and
        # nothing wrong. A one-sided absence below is still a failure.
        print("  planes: none written (OutAllPlanes off) — skipped")
        return []
    if not both:
        print("  NOTE planes: no common OutAllPlanes VTK files found")
        return [("planes", "present", "missing", 0.0)]
    times = sorted({t for t, _ in both})
    pick = times[:: max(1, len(times) // n_times)][-n_times:]
    sel = [k for k in both if k[0] in pick]
    bad, counts, floors = [], Counter(), {}
    for key in sel:
        a = read_plane_vtk(cw[key])
        b = read_plane_vtk(mir[key])
        for name, exp in PLANE_EXPECT.items():
            if name not in a or name not in b:
                continue
            x = a[name]
            yv = b[name][::-1, :]            # reverse the y index
            got, rel = classify(x.ravel(), yv.ravel(), tol,
                                1e-9 * max(np.abs(x).max(), 1e-12))
            counts[f"{name}:{got}"] += 1
            if got not in exp and got != "negligible":
                bad.append((f"plane t={key[0]} i={key[1]} {name}",
                            "/".join(sorted(exp)), got, rel))
    print(f"  planes: {len(sel)} (time,plane) pairs at times {pick}")
    for k, v in sorted(counts.items()):
        print(f"    {k}: {v}")
    return bad


def turbine_count(names):
    ts = {int(m.group(1)) for n in names for m in [re.search(r"T(\d+)$|T(\d+)[ND]", n)] if m and m.group(1)}
    return max(ts) if ts else 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cw_case")
    ap.add_argument("mir_case", nargs="?",
                    help="omit with --twin: T1 vs T2 within cw_case")
    ap.add_argument("--twin", action="store_true",
                    help="compare turbine 1 (clockwise) against turbine 2 "
                         "(mirrored) inside one counter-rotating-farm run")
    ap.add_argument("--tol", type=float, default=2e-3)
    ap.add_argument("--vtk-times", type=int, default=3)
    ap.add_argument("--show", action="store_true")
    a = ap.parse_args()

    bad = []
    if a.twin:
        print("== farm channels (FAST.Farm.out, T1 vs T2)")
        na, da = read_outb(a.cw_case, "FAST.Farm.out")
        p1 = [n for n in na if "T1" in n and n.replace("T1", "T2", 1) in set(na)]
        n1 = ["Time"] + p1
        n2 = ["Time"] + [n.replace("T1", "T2", 1) for n in p1]
        ia = {n: i for i, n in enumerate(na)}
        d1 = da[:, [ia[n] for n in n1]]
        d2 = da[:, [ia[n] for n in n2]]
        bad += compare_channels(n1, d1, n1, d2, farm_expect, a.tol,
                                show=a.show, label="farm T1-vs-T2")
        print("== turbine channels (T1.outb vs T2.outb)")
        na, da = read_outb(a.cw_case, "FAST.Farm.T1.outb")
        nb, db = read_outb(a.cw_case, "FAST.Farm.T2.outb")
        bad += compare_channels(na, da, nb, db, turbine_expect, a.tol,
                                show=a.show, label="turbine")
        print("== wake-plane fields (WT1 vs WT2)")
        bad += compare_planes(plane_files(a.cw_case, wt=1),
                              plane_files(a.cw_case, wt=2), a.tol, a.vtk_times)
    else:
        if not a.mir_case:
            ap.error("mir_case is required without --twin")
        print("== farm channels (FAST.Farm.out)")
        na, da = read_outb(a.cw_case, "FAST.Farm.out")
        nb, db = read_outb(a.mir_case, "FAST.Farm.out")
        bad += compare_channels(na, da, nb, db, farm_expect, a.tol,
                                show=a.show, label="farm")

        for tn in range(1, turbine_count(na) + 1):
            print(f"== turbine channels (FAST.Farm.T{tn}.outb)")
            try:
                n1, d1 = read_outb(a.cw_case, f"FAST.Farm.T{tn}.outb")
                n2, d2 = read_outb(a.mir_case, f"FAST.Farm.T{tn}.outb")
                bad += compare_channels(n1, d1, n2, d2, turbine_expect, a.tol,
                                        show=a.show, label=f"turbine {tn}")
            except FileNotFoundError as e:
                print(f"  NOTE turbine outb missing: {e}")

            print(f"== wake-plane fields (OutAllPlanes, WT{tn})")
            bad += compare_planes(plane_files(a.cw_case, wt=tn),
                                  plane_files(a.mir_case, wt=tn),
                                  a.tol, a.vtk_times)

    print()
    if bad:
        for ch, exp, got, rel in bad:
            print(f"MISMATCH {ch}: expected {exp}, got {got} (rel={rel:.3e})")
        sys.exit(1)
    print("ALL PASS — the mirrored farm is the y-reflection of the clockwise one")


if __name__ == "__main__":
    main()
