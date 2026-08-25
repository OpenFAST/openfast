#!/usr/bin/env python3
"""Phase 1.2/1.3: build and run the OpenFAST clockwise/mirrored pair.

Builds the pair in scratch from the 5MW land baseline so r-test stays clean until
the case is settled. ServoDyn is off because the mirrored ServoDyn interface is
Phase 3 and is still guard-railed, so the rotor runs at a fixed speed and the
comparison isolates the aerodynamic and structural mirror.

By default this reports the *measured* symmetry of every channel rather than
asserting a table, which is how the sign map in CONVENTIONS is meant to be built.
Pass --check to assert the expectations in ED_EXPECT once they are settled.

Usage: openfast_mirror.py [--check] [--tmax 20]
"""
import os
import re
import shutil
import subprocess
import sys

import numpy as np

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REPO = REPO_ROOT
RTEST = f"{REPO}/reg_tests/r-test/glue-codes/openfast"
SRC = f"{RTEST}/5MW_Land_BD_DLL_WTurb" if "--beamdyn" in sys.argv else f"{RTEST}/5MW_Land_DLL_WTurb"
BASE = f"{RTEST}/5MW_Baseline"
EXE = f"{REPO}/build-docker-double-debug/glue-codes/openfast/openfast"
WORK = os.path.join(os.environ.get("TMPDIR", "/tmp"), "of_mirror")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_mirror import classify, blade_permutation  # noqa: E402

RIGID = "--rigid" in sys.argv
# Free drivetrain: exercises the gearbox efficiency factor and shaft torque paths,
# which a fixed rotor speed never reaches.
GENDOF = "--gendof" in sys.argv
# Loose coupling is needed to reach SignLSSTrq, but BeamDyn is unstable with it, so it
# is requested separately rather than being implied by a free drivetrain.
LOOSE = "--loose" in sys.argv

# Vertical shear is symmetric about the mirror plane, so both decks get the same
# exponent. Nacelle yaw is a global-frame orientation rather than a rotor-convention
# input, so the mirrored deck gets the opposite sign.
# BeamDyn blades instead of the ElastoDyn beam model.
BEAMDYN = "--beamdyn" in sys.argv
SHEAR = float(sys.argv[sys.argv.index("--shear") + 1]) if "--shear" in sys.argv else None
YAW = float(sys.argv[sys.argv.index("--yaw") + 1]) if "--yaw" in sys.argv else None
# Free yaw exercises the yaw bearing dynamics and the yaw moment mirror.
YAWDOF = "--yawdof" in sys.argv

# With every structural DOF locked the BEMT iteration converges to slightly
# different round-off in the two runs; the real configurations hold 1e-6.
TOL = 1e-5 if RIGID else 1e-6

# One channel per line: the AeroDyn OutList parser only takes the first quoted
# name on each line and silently drops the rest.
AD_OUTLIST = """"RtAeroFxh"
"RtAeroMxh"
"RtTSR"
"RtSkew"
"B1N1Alpha"
"B1N2Alpha"
"B1N3Alpha"
"B1N1Theta"
"B1N2Theta"
"B1N3Theta"
"B1N1Phi"
"B1N2Phi"
"B1N3Phi"
"B1N1VDisx"
"B1N2VDisx"
"B1N3VDisx"
"B1N1VDisy"
"B1N2VDisy"
"B1N3VDisy"
"B1N1Fn"
"B1N2Fn"
"B1N3Fn"
"B1N1Ft"
"B1N2Ft"
"B1N3Ft"
"B1N1Cl"
"B1N2Cl"
"B1N3Cl"
"B1N1Curve"
"B1N2Curve"
"B1N3Curve"
"""

OUTLIST = """"RotSpeed"    - Rotor speed
"RotAccel"    - Rotor acceleration
"Azimuth"     - Blade 1 azimuth
"GenSpeed"    - Generator speed
"LSShftFxa"   - LSS thrust
"LSShftFya"
"LSShftFza"
"LSShftMxa"   - LSS torque
"LSSTipVxa"
"LSSTipAxa"
"LSSTipPxa"
"LSShftTq"
"HSShftTq"
"HSShftPwr"
"HSSBrTq"
"LSSGagMxa"
"LSSGagPxa"
"GenAccel"
"LSSTipMya"
"LSSTipMza"
"RotPwr"      - Rotor power
"RotThrust"   - Rotor thrust
"RotTorq"     - Rotor torque
"YawBrFxp"
"YawBrFyp"
"YawBrFzp"
"YawBrMxp"
"YawBrMyp"
"YawBrMzp"
"OoPDefl1"    - Blade 1 out-of-plane tip deflection
"IPDefl1"     - Blade 1 in-plane tip deflection
"TipDxc1"
"TipDyc1"
"TipDzc1"
"RootFxc1"
"RootFyc1"
"RootFzc1"
"RootMxc1"
"RootMyc1"
"RootMzc1"
"OoPDefl2"
"IPDefl2"
"RootMxc2"
"RootMyc2"
"OoPDefl3"
"IPDefl3"
"RootMxc3"
"RootMyc3"
"TwrBsFxt"
"TwrBsFyt"
"TwrBsMxt"
"TwrBsMyt"
"TwrBsMzt"
"""


def set_var(text, var, value):
    pat = re.compile(r"^(\s*)(\S+)(\s+" + re.escape(var) + r"\s)", re.MULTILINE)
    new, n = pat.subn(lambda m: f"{m.group(1)}{value}{m.group(3)}", text, count=1)
    if n != 1:
        raise KeyError(f"could not set {var!r}")
    return new


def build(mirrored, tmax):
    name = "mir" if mirrored else "cw"
    # The baseline files are referenced as ../5MW_Baseline/..., so stage each case as a
    # sibling of a symlink to the real one.
    os.makedirs(WORK, exist_ok=True)
    link = os.path.join(WORK, "5MW_Baseline")
    if not os.path.exists(link):
        os.symlink(BASE, link)
    d = os.path.join(WORK, name)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    for f in os.listdir(SRC):
        if f.endswith(".dat"):
            shutil.copy(os.path.join(SRC, f), d)

    # ElastoDyn: fixed rotor speed (no controller), keep the flexible DOFs on so
    # the structural side of the mirror is actually exercised.
    ed_src = ("NRELOffshrBsline5MW_Onshore_ElastoDyn_BDoutputs.dat" if BEAMDYN
              else "NRELOffshrBsline5MW_Onshore_ElastoDyn.dat")
    ed = open(os.path.join(d, ed_src)).read()
    ed = set_var(ed, "GenDOF", "True" if GENDOF else "False")
    if GENDOF:
        # The baseline gearbox is 100% efficient, which makes GBoxEffFac 1.0 whichever
        # branch is taken and hides any error in the power-flow direction test.
        ed = set_var(ed, "GBoxEff", "95.0")
    ed = set_var(ed, "YawDOF", "True" if YAWDOF else "False")
    ed = set_var(ed, "RotSpeed", "9.15")
    if YAW is not None:
        ed = set_var(ed, "NacYaw", str(-YAW if mirrored else YAW))
    ed = set_var(ed, "BlPitch(1)", "1.0")
    ed = set_var(ed, "BlPitch(2)", "1.0")
    ed = set_var(ed, "BlPitch(3)", "1.0")
    if RIGID:
        # Bisection aid: with every structural DOF off, only the geometry handed to
        # AeroDyn can differ between the two runs.
        for dof in ("FlapDOF1", "FlapDOF2", "EdgeDOF", "DrTrDOF",
                    "TwFADOF1", "TwFADOF2", "TwSSDOF1", "TwSSDOF2"):
            ed = set_var(ed, dof, "False")
    # Swap only the channel list, keeping the optional nodal-outputs section intact.
    head, _, rest = ed.partition("OutList")
    _, _, tail = rest.partition("END of OutList section")
    ed = (head + "OutList - Output channels\n" + OUTLIST
          + "END of OutList section" + tail)
    open(os.path.join(d, "ElastoDyn.dat"), "w").write(ed)

    ad_name = "NRELOffshrBsline5MW_Onshore_AeroDyn.dat"
    ad = open(os.path.join(d, ad_name)).read()
    ad = set_var(ad, "NBlOuts", "3")
    ad = set_var(ad, "SumPrint", "True")
    head, _, rest = ad.partition("OutList ")
    _, _, tail = rest.partition("END of OutList section")
    ad = (head + "OutList - Output channels\n" + AD_OUTLIST
          + "END of OutList section" + tail)
    open(os.path.join(d, ad_name), "w").write(ad)

    fst_name = "5MW_Land_BD_DLL_WTurb.fst" if BEAMDYN else "5MW_Land_DLL_WTurb.fst"
    fst = open(os.path.join(SRC, fst_name)).read()
    fst = set_var(fst, "TMax", str(tmax))
    fst = set_var(fst, "CompServo", "0")
    if LOOSE:
        # SignLSSTrq, and so the gearbox efficiency direction, is only reached from
        # the loose-coupling integrators.
        fst = set_var(fst, "ModCoupling", "1")
    fst = set_var(fst, "CompInflow", "1")
    inflow = f"{BASE}/NRELOffshrBsline5MW_InflowWind_Steady8mps.dat"
    if SHEAR is not None:
        iw = open(inflow).read()
        iw = set_var(iw, "PLExp", str(SHEAR))
        inflow = os.path.join(d, "InflowWind.dat")
        open(inflow, "w").write(iw)
        inflow = "InflowWind.dat"
    fst = set_var(fst, "MirrorRotor", "True" if mirrored else "False")
    fst = set_var(fst, "EDFile", '"ElastoDyn.dat"')
    fst = set_var(fst, "InflowFile",
                  f'"{inflow}"')
    fst = set_var(fst, "AeroFile", '"NRELOffshrBsline5MW_Onshore_AeroDyn.dat"')
    fst = set_var(fst, "OutFileFmt", "1")
    fst = set_var(fst, "Echo", "False")
    open(os.path.join(d, f"{name}.fst"), "w").write(fst)

    p = subprocess.run([EXE, f"{name}.fst"], cwd=d, capture_output=True, text=True)
    out = os.path.join(d, f"{name}.out")
    if not os.path.exists(out):
        tail = [l for l in (p.stdout + p.stderr).split("\n") if l.strip()][-8:]
        print("\n".join(tail))
        sys.exit(f"{name} run failed")
    return out


def load(path):
    lines = open(path, errors="replace").read().split("\n")
    h = next(i for i, l in enumerate(lines) if l.strip().startswith("Time"))
    names = lines[h].split()
    rows = [[float(x) for x in l.split()] for l in lines[h + 2:] if l.strip()]
    return names, np.array(rows)


def main():
    tmax = 20
    if "--tmax" in sys.argv:
        tmax = float(sys.argv[sys.argv.index("--tmax") + 1])
    shutil.rmtree(WORK, ignore_errors=True)
    cw, mir = build(False, tmax), build(True, tmax)
    n1, a = load(cw)
    n2, b = load(mir)
    idx = {k: v for v, k in enumerate(n1)}
    sl = slice(int(0.5 * len(a)), None)
    # Solver diagnostics are not physical and can be orders of magnitude larger than
    # any load, which would drag the noise floor up and mask real channels.
    DIAG = {"ConvError", "ConvIter", "NumUJac", "Time"}
    phys = [idx[c] for c in n1 if c not in DIAG]
    floor = 1e-9 * np.abs(a[sl][:, phys]).max()

    groups = {"S": [], "F": [], "A": [], "negligible": [], "?": []}
    if "--dump" in sys.argv:
        want = sys.argv[sys.argv.index("--dump") + 1].split(",")
        rows = list(range(0, len(a), max(1, len(a) // 12)))
        for ch in want:
            print(f"\n{ch}:")
            print(f"{'t':>8s} {'CW':>14s} {'MIRROR':>14s}")
            for r in rows:
                print(f"{a[r,0]:8.3f} {a[r,idx[ch]]:14.6g} {b[r,idx[ch]]:14.6g}")
        return

    print(f"{'channel':14s} {'peak(CW)':>12s} {'class':>10s} {'resid':>10s}")
    for ch in n1:
        if ch == "Time":
            continue
        mate = blade_permutation(ch)
        if mate not in idx:
            mate = ch
        got, rel = classify(a[sl, idx[ch]], b[sl, idx[mate]], TOL, floor)
        if ch in DIAG:
            got = "diagnostic"
            groups.setdefault("diagnostic", []).append(ch)
            continue
        groups[got].append(ch)
        print(f"{ch:14s} {np.abs(a[sl, idx[ch]]).max():12.5g} "
              f"{got:>10s} {rel:10.2e}")

    print("\n--- measured summary (blade permutation applied) ---")
    for k in ("S", "F", "A", "negligible", "?"):
        if groups[k]:
            print(f"{k:11s} ({len(groups[k]):2d}): {' '.join(groups[k])}")
    if groups["?"]:
        print(f"\n{len(groups['?'])} channel(s) neither same, flipped nor mirrored angle")
        sys.exit(1)
    print("\nEvery channel resolves to same, flipped, mirrored angle, or negligible")


if __name__ == "__main__":
    main()
