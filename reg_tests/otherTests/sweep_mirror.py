#!/usr/bin/env python3
"""Phase 1.0a sweep: run the CW/mirrored driver pair across a matrix of conditions.

Each variant edits the base driver deck and AeroDyn primary file, runs the pair,
and requires every channel to resolve to same / flipped / mirrored-angle /
negligible. A single condition proves very little on its own; the point of the
matrix is that each row exercises a different assumption in the sign map.

Quantities defined in the global frame (nacelle yaw) are not rotor-convention
inputs, so the mirrored deck must be given the mirrored value explicitly. That is
what "mirror" entries do. Everything else stays identical between the pair,
which is the whole claim being tested.

Usage: sweep_mirror.py [name ...]      (no args runs the full matrix)
"""
import os
import re
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_mirror import compare  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(HERE))
BASE = os.path.join(os.path.dirname(HERE), "driver_checks")
DRIVER = os.path.join(REPO_ROOT, "build-docker-double-debug",
                       "modules", "aerodyn", "aerodyn_driver")
WORK = os.path.join(os.environ.get("TMPDIR", "/tmp"), "mirror_sweep")

# name -> (dvr edits, primary edits, mirrored-deck overrides, blade permutation[, tol])
# Tolerance defaults to 1e-6; raise it only where a stateful model makes exact
# agreement unreasonable, and say why.
MATRIX = [
    ("base",            {}, {}, {}, False),
    ("pitch_+10",       {"BldPitch(1)": "10.0"}, {}, {}, False),
    ("pitch_-3",        {"BldPitch(1)": "-3.0"}, {}, {}, False),
    ("wind_4",          {"HWindSpeed": "4.0"}, {}, {}, False),
    ("wind_20_stall",   {"HWindSpeed": "20.0"}, {}, {}, False),
    ("highTSR_25rpm",   {"RotSpeed(1)": "25.0"}, {}, {}, False),
    ("lowTSR_3rpm",     {"RotSpeed(1)": "3.0"}, {}, {}, False),
    ("propbrake",       {"HWindSpeed": "25.0", "RotSpeed(1)": "3.0",
                         "BldPitch(1)": "-10.0"}, {}, {}, False),
    # Shaft tilt skews the inflow, so the rotor sees azimuthal variation and the
    # blades no longer pair up index for index.
    ("tilt_5",          {"ShftTilt(1)": "-5.0"}, {}, {}, True),
    ("precone_2.5",     {"Precone(1)": "2.5"}, {}, {}, False),
    # Vertical shear is symmetric about the mirror plane, so the deck is shared;
    # the blades still swap sweep order.
    ("shear_0.2",       {"PLExp": "0.2"}, {}, {}, True),
    # Yaw is a global-frame orientation, so the mirrored run gets -yaw.
    ("yaw_20",          {"NacYaw(1)": "20.0"}, {}, {"NacYaw(1)": "-20.0"}, True),
    ("yaw_20_skew",     {"NacYaw(1)": "20.0"}, {"Skew_Mod": "1"},
                        {"NacYaw(1)": "-20.0"}, True),
    ("yaw_-20",         {"NacYaw(1)": "-20.0"}, {}, {"NacYaw(1)": "20.0"}, True),
    ("wake_off",        {}, {"Wake_Mod": "0"}, {}, False),
    ("bem_polar",       {}, {"BEM_Mod": "2"}, {}, False),
    # The polar BEM path measures skew with its own psiSkewOffset, built from a
    # separate cross product, so it needs covering independently of BEM_Mod=1.
    ("bem_polar_skew",  {"NacYaw(1)": "20.0"}, {"BEM_Mod": "2", "Skew_Mod": "1"},
                        {"NacYaw(1)": "-20.0"}, True),
    ("bem_polar_shear", {"PLExp": "0.2"}, {"BEM_Mod": "2"}, {}, True),
    ("dbemt_2",         {}, {"DBEMT_Mod": "2"}, {}, False),
    ("ua_3",            {}, {"UA_Mod": "3"}, {}, False),
    ("ua_6_oye",        {}, {"UA_Mod": "6"}, {}, False),
    ("ua_4_hgm",        {}, {"UA_Mod": "4"}, {}, False),
    # Shear + yaw + unsteady aero is the stiffest combination in the matrix. UA carries
    # state and is solved iteratively, so the two runs converge to slightly different
    # round-off; the residual here is ~1.5 N-m on a 9e5 N-m moment.
    ("shear_yaw_ua",    {"PLExp": "0.2", "NacYaw(1)": "15.0"},
                        {"UA_Mod": "3", "Skew_Mod": "1"},
                        {"NacYaw(1)": "-15.0"}, True, 1e-5),
]


def set_var(text, var, value):
    """Replace the value column on the line whose comment names `var`."""
    pat = re.compile(r"^(\s*)(\S+)(\s+" + re.escape(var) + r"\s)",
                     re.MULTILINE)
    new, n = pat.subn(lambda m: f"{m.group(1)}{value}{m.group(3)}", text, count=1)
    if n != 1:
        raise KeyError(f"could not set {var!r}")
    return new


def build_case(dst, mirrored, dvr_edits, pri_edits, overrides):
    os.makedirs(dst, exist_ok=True)
    for f in os.listdir(BASE):
        if f.endswith((".dat", ".dvr", ".csv", ".ipt")):
            shutil.copy(os.path.join(BASE, f), dst)

    pri = open(os.path.join(BASE, "ad_primary.dat")).read()
    for k, v in pri_edits.items():
        pri = set_var(pri, k, v)
    open(os.path.join(dst, "ad_primary.dat"), "w").write(pri)

    dvr = open(os.path.join(BASE, "ad_CW.dvr")).read()
    edits = dict(dvr_edits)
    if mirrored:
        edits.update(overrides)
    for k, v in edits.items():
        dvr = set_var(dvr, k, v)
    dvr = set_var(dvr, "MirrorRotor(1)", "True" if mirrored else "False")
    name = "mir" if mirrored else "cw"
    open(os.path.join(dst, f"{name}.dvr"), "w").write(dvr)
    return os.path.join(dst, f"{name}.dvr")


def run(dvr_path):
    d, f = os.path.split(dvr_path)
    p = subprocess.run([DRIVER, f], cwd=d, capture_output=True, text=True)
    out = os.path.join(d, f[:-4] + ".out")
    if p.returncode != 0 or not os.path.exists(out):
        tail = [l for l in (p.stdout + p.stderr).split("\n") if l.strip()][-2:]
        return None, " | ".join(tail)
    return out, ""


def main():
    wanted = sys.argv[1:]
    rows = [r for r in MATRIX if not wanted or r[0] in wanted]
    shutil.rmtree(WORK, ignore_errors=True)

    fails = 0
    print(f"{'variant':18s} {'chk':>4s}  result")
    for row in rows:
        name, dvr_e, pri_e, over, permute = row[:5]
        tol = row[5] if len(row) > 5 else 1e-6
        d = os.path.join(WORK, name)
        try:
            cw = build_case(d, False, dvr_e, pri_e, over)
            mir = build_case(d, True, dvr_e, pri_e, over)
        except KeyError as exc:
            print(f"{name:18s} {'-':>4s}  SETUP  {exc}")
            fails += 1
            continue

        cw_out, err = run(cw)
        if cw_out is None:
            print(f"{name:18s} {'-':>4s}  CW RUN FAILED  {err}")
            fails += 1
            continue
        mir_out, err = run(mir)
        if mir_out is None:
            print(f"{name:18s} {'-':>4s}  MIRROR RUN FAILED  {err}")
            fails += 1
            continue

        checked, bad, unknown, _ = compare(cw_out, mir_out, tol=tol,
                                           permute=permute)
        tag = "perm" if permute else ""
        if tol != 1e-6:
            tag += f" tol={tol:g}"
        if bad or unknown:
            fails += 1
            print(f"{name:18s} {checked:4d}  FAIL {tag}")
            for ch, rel in unknown[:4]:
                print(f"{'':18s}       unclassified {ch} rel={rel:.2e}")
            for ch, exp, got, rel in bad[:4]:
                print(f"{'':18s}       {ch} expected {exp} got {got} rel={rel:.2e}")
            extra = len(bad) + len(unknown) - 8
            if extra > 0:
                print(f"{'':18s}       ... and {extra} more")
        else:
            print(f"{name:18s} {checked:4d}  PASS {tag}")

    print(f"\n{len(rows) - fails}/{len(rows)} variants pass")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
