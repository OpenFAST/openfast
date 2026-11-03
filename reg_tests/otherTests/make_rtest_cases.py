#!/usr/bin/env python3
"""Generate the Phase 1 mirrored-rotor regression cases in r-test.

Creates 5MW_Land_noDLL_Steady_CW and 5MW_Land_noDLL_Steady_MirrorRotor from the
5MW land baseline. The pair is identical apart from the MirrorRotor flag, which
is the whole point: the mirrored deck describes the same clockwise turbine.

Case design, and why each choice is there:
  CompServo = 0        ServoDyn is still guard-railed for mirrored rotors
  GenDOF    = True     free drivetrain, so shaft torque actually matters
  ModCoupling = 1      SignLSSTrq is only reached from the loose-coupling
                       integrators, so tight coupling cannot exercise the
                       gearbox efficiency direction at all
  GBoxEff   = 95       the baseline 100% makes both branches of the efficiency
                       factor identical and hides errors in that direction test
  ShftTilt  = -5       inherited from the baseline; gives azimuthal variation so
                       the blades are not interchangeable
  steady 8 m/s, TMax 20 s, flexible blades and tower

Usage: make_rtest_cases.py [--write]
"""
import os
import re
import shutil
import sys

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

RTEST = os.path.join(REPO_ROOT, "reg_tests", "r-test", "glue-codes", "openfast")
ED_SRC = os.path.join(RTEST, "5MW_Land_DLL_WTurb")
BD_SRC = os.path.join(RTEST, "5MW_Land_BD_DLL_WTurb")

# name -> (mirrored, beamdyn)
CASES = {"5MW_Land_noDLL_Steady_CW":           (False, False),
         "5MW_Land_noDLL_Steady_MirrorRotor":    (True, False),
         "5MW_Land_BD_noDLL_Steady_CW":         (False, True),
         "5MW_Land_BD_noDLL_Steady_MirrorRotor": (True, True)}

OUTLIST = """"RotSpeed"
"RotAccel"
"Azimuth"
"GenSpeed"
"GenAccel"
"LSSTipVxa"
"LSSTipAxa"
"LSSTipPxa"
"LSShftMxa"
"LSSGagMxa"
"RotTorq"
"LSShftTq"
"HSShftTq"
"HSShftPwr"
"RotPwr"
"RotThrust"
"LSShftFxa"
"LSShftFya"
"LSShftFza"
"LSSTipMya"
"LSSTipMza"
"YawBrFxp"
"YawBrFyp"
"YawBrFzp"
"YawBrMxp"
"YawBrMyp"
"YawBrMzp"
"OoPDefl1"
"IPDefl1"
"TipDxc1"
"TipDyc1"
"RootFxc1"
"RootFyc1"
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

AD_OUTLIST = """"RtAeroFxh"
"RtAeroMxh"
"RtAeroPwr"
"RtTSR"
"RtSkew"
"B1N1Alpha"
"B1N2Alpha"
"B1N3Alpha"
"B1N1Theta"
"B1N2Theta"
"B1N3Theta"
"B1N1Fn"
"B1N2Fn"
"B1N3Fn"
"B1N1Ft"
"B1N2Ft"
"B1N3Ft"
"""


def set_var(text, var, value):
    pat = re.compile(r"^(\s*)(\S+)(\s+" + re.escape(var) + r"\s)", re.MULTILINE)
    new, n = pat.subn(lambda m: f"{m.group(1)}{value}{m.group(3)}", text, count=1)
    if n != 1:
        raise KeyError(f"could not set {var!r}")
    return new


def swap_outlist(text, entries, marker="OutList"):
    head, _, rest = text.partition(marker)
    _, _, tail = rest.partition("END of OutList section")
    return head + marker + " - Output channels\n" + entries + "END of OutList section" + tail


def build(case, mirrored, beamdyn, write):
    src = BD_SRC if beamdyn else ED_SRC
    dst = os.path.join(RTEST, case)
    if write:
        # Only the generated files are replaced. Hand-written README.md files and the
        # committed .outb baselines are left alone, so re-running this does not quietly
        # delete them.
        os.makedirs(dst, exist_ok=True)

    ed_name = ("NRELOffshrBsline5MW_Onshore_ElastoDyn_BDoutputs.dat" if beamdyn
               else "NRELOffshrBsline5MW_Onshore_ElastoDyn.dat")
    ed = open(os.path.join(src, ed_name)).read()
    ed = set_var(ed, "GenDOF", "True")
    ed = set_var(ed, "YawDOF", "False")
    ed = set_var(ed, "RotSpeed", "9.15")
    ed = set_var(ed, "GBoxEff", "95.0")
    for k in (1, 2, 3):
        ed = set_var(ed, f"BlPitch({k})", "1.0")
    ed = swap_outlist(ed, OUTLIST)

    ad = open(os.path.join(src, "NRELOffshrBsline5MW_Onshore_AeroDyn.dat")).read()
    ad = set_var(ad, "NBlOuts", "3")
    ad = swap_outlist(ad, AD_OUTLIST, "OutList ")

    fst_name = ("5MW_Land_BD_DLL_WTurb.fst" if beamdyn else "5MW_Land_DLL_WTurb.fst")
    fst = open(os.path.join(src, fst_name)).read()
    fst = set_var(fst, "TMax", "20")
    fst = set_var(fst, "Echo", "False")
    # BeamDyn is unstable with loose coupling, so only the ElastoDyn pair uses it. That
    # costs the BeamDyn pair the gearbox-efficiency path, which the ElastoDyn pair covers.
    if not beamdyn:
        fst = set_var(fst, "ModCoupling", "1")
    fst = set_var(fst, "CompServo", "0")
    fst = set_var(fst, "CompInflow", "1")
    fst = set_var(fst, "MirrorRotor", "True" if mirrored else "False")
    fst = set_var(fst, "EDFile", '"ElastoDyn.dat"')
    fst = set_var(fst, "AeroFile", '"AeroDyn.dat"')
    fst = set_var(fst, "InflowFile",
                  '"../5MW_Baseline/NRELOffshrBsline5MW_InflowWind_Steady8mps.dat"')
    fst = set_var(fst, "OutFileFmt", "2")

    files = {
        "ElastoDyn.dat": ed,
        "AeroDyn.dat": ad,
        f"{case}.fst": fst,
    }
    if write:
        for name, text in files.items():
            open(os.path.join(dst, name), "w").write(text)
        shutil.copy(os.path.join(src, "NRELOffshrBsline5MW_Onshore_ElastoDyn_Tower.dat"),
                    os.path.join(dst, "ElastoDyn_Tower.dat"))
        ed2 = open(os.path.join(dst, "ElastoDyn.dat")).read()
        ed2 = set_var(ed2, "TwrFile", '"ElastoDyn_Tower.dat"')
        open(os.path.join(dst, "ElastoDyn.dat"), "w").write(ed2)
    print(f"  {case:40s} MirrorRotor={'True' if mirrored else 'False':5s} "
          f"{'BeamDyn' if beamdyn else 'ElastoDyn':9s} "
          f"{'written' if write else '(dry run)'}")


def main():
    write = "--write" in sys.argv
    for case, (mirrored, beamdyn) in CASES.items():
        build(case, mirrored, beamdyn, write)
    if not write:
        print("\ndry run; pass --write to create the cases")


if __name__ == "__main__":
    main()
