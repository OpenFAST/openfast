#!/usr/bin/env python3
"""Phase 1.0b canary: the HSS brake must behave identically on a mirrored rotor.

The brake is the one drivetrain load that is signed by the direction the shaft is
actually turning rather than by a convention, so it is the easiest place to get a
mirrored rotor wrong: a blanket sign flip makes the brake drive the rotor instead
of stopping it. This runs the SED HSS-brake case clockwise and mirrored and
requires the rotor-convention channels to be identical, not merely similar.

Usage: sed_canary.py
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
CASE = f"{REPO}/reg_tests/r-test/modules/simple-elastodyn/sed_test_HSSbrk"
DRIVER = f"{REPO}/build-docker-double-debug/modules/simple-elastodyn/sed_driver"
WORK = os.path.join(os.environ.get("TMPDIR", "/tmp"), "sed_canary")

# RotSpeed and Azimuth are reported in the rotor's own convention, so a mirrored
# rotor spinning its design direction must read exactly the same as its CW twin.
IDENTICAL = ["Azimuth", "RotSpeed", "RotAcc", "GenSpeed", "GenAcc"]

# The same two quantities carrying an explicit axis suffix are the physical components
# about the shaft x axis, so they must come out exactly opposite. Requesting them is the
# only way the split is visible at all -- the case ships with the rotor-convention names
# alone, and those read identically whether or not the split exists.
FLIPPED = ["LSSTipVxa", "LSSTipAxa"]


def load(path):
    lines = open(path, errors="replace").read().split("\n")
    h = next(i for i, l in enumerate(lines) if l.strip().startswith("Time"))
    names = lines[h].split()
    rows = [[float(x) for x in l.split()] for l in lines[h + 2:] if l.strip()]
    return names, np.array(rows)


def build(mirrored):
    name = "mir" if mirrored else "cw"
    d = os.path.join(WORK, name)
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(CASE, d)

    inp = os.path.join(d, "sed_primary.inp")
    txt = open(inp).read()
    txt = txt.replace('"RotTorq"\n', '"RotTorq"\n' + "".join(f'"{c}"\n' for c in FLIPPED))
    open(inp, "w").write(txt)

    dvr = os.path.join(d, "sed_driver.dvr")
    txt = open(dvr).read()
    txt = re.sub(r"^(\s*)\S+(\s+MirrorRotor\s)",
                 lambda m: f"{m.group(1)}{'True' if mirrored else 'False'}{m.group(2)}",
                 txt, count=1, flags=re.MULTILINE)
    txt = re.sub(r'^("?)[^"\n]*\1(\s+OutRootName\s)',
                 lambda m: f'"{name}"{m.group(2)}', txt, count=1, flags=re.MULTILINE)
    open(dvr, "w").write(txt)
    p = subprocess.run([DRIVER, "sed_driver.dvr"], cwd=d,
                       capture_output=True, text=True)
    out = os.path.join(d, f"{name}.out")
    if not os.path.exists(out):
        print((p.stdout + p.stderr)[-1500:])
        sys.exit(f"{name} run failed")
    return out


def main():
    cw, mir = build(False), build(True)
    n1, a = load(cw)
    n2, b = load(mir)
    idx = {k: v for v, k in enumerate(n1)}

    print(f"{'channel':10s} {'max|CW|':>12s} {'max abs diff':>14s}  verdict")
    bad = 0
    for ch in IDENTICAL:
        if ch not in idx:
            continue
        x, y = a[:, idx[ch]], b[:, idx[ch]]
        diff = np.abs(y - x).max()
        ok = diff == 0.0
        bad += not ok
        print(f"{ch:10s} {np.abs(x).max():12.5g} {diff:14.3e}  "
              f"{'IDENTICAL' if ok else 'DIFFERS'}")

    for ch in FLIPPED:
        if ch not in idx:
            print(f"{ch:10s} {'':>12s} {'':>14s}  MISSING - split not exercised")
            bad += 1
            continue
        x, y = a[:, idx[ch]], b[:, idx[ch]]
        diff = np.abs(y + x).max()
        ok = diff == 0.0 and np.abs(x).max() > 0.0
        bad += not ok
        print(f"{ch:10s} {np.abs(x).max():12.5g} {diff:14.3e}  "
              f"{'EXACTLY OPPOSITE' if ok else 'NOT OPPOSITE'}")

    # The brake must remove energy in both runs, whichever way the shaft turns.
    sp = a[:, idx["RotSpeed"]]
    spm = b[:, idx["RotSpeed"]]
    print(f"\nCW  RotSpeed {sp[0]:.4f} -> {sp[-1]:.4f} rpm")
    print(f"MIR RotSpeed {spm[0]:.4f} -> {spm[-1]:.4f} rpm")
    if abs(spm[-1]) > abs(spm[0]):
        print("FAIL: mirrored rotor sped up - the brake is driving it")
        bad += 1

    print("\nCANARY PASS" if not bad else "\nCANARY FAIL")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
