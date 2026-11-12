#!/usr/bin/env python3
"""Phase 1.5 acceptance: run the OpenFAST mirror check across the condition matrix.

Each row exercises a different assumption in the sign map. A single condition
proves very little on its own; the point is that they all have to hold at once.

Usage: openfast_mirror_matrix.py [name ...]
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RUNNER = os.path.join(HERE, "openfast_mirror.py")

MATRIX = [
    ("rigid",              ["--rigid"]),
    ("base_flexible",      []),
    ("free_drivetrain",    ["--gendof", "--loose"]),
    ("shear",              ["--shear", "0.2"]),
    ("yaw_+15",            ["--yaw", "15"]),
    ("yaw_-15",            ["--yaw", "-15"]),
    ("shear_yaw",          ["--shear", "0.2", "--yaw", "15"]),
    ("free_yaw",           ["--yawdof", "--yaw", "15"]),
    ("shear_freeyaw_dt",   ["--shear", "0.2", "--yaw", "10", "--yawdof", "--gendof", "--loose"]),
    # BeamDyn blades. Not run with --loose: BeamDyn needs tight coupling.
    ("bd_base",            ["--beamdyn", "--tmax", "10"]),
    ("bd_shear",           ["--beamdyn", "--tmax", "10", "--shear", "0.2"]),
    ("bd_yaw",             ["--beamdyn", "--tmax", "10", "--yaw", "15"]),
    ("bd_free_drivetrain", ["--beamdyn", "--tmax", "10", "--gendof"]),
]


def main():
    wanted = sys.argv[1:]
    rows = [r for r in MATRIX if not wanted or r[0] in wanted]
    fails = 0
    for name, args in rows:
        p = subprocess.run([sys.executable, RUNNER] + args,
                           capture_output=True, text=True)
        tail = [l for l in p.stdout.split("\n") if l.strip()]
        verdict = tail[-1] if tail else "(no output)"
        ok = p.returncode == 0
        fails += not ok
        print(f"  {name:20s} {'PASS' if ok else 'FAIL'}  {verdict[:70]}", flush=True)
        if not ok:
            for l in tail[-6:]:
                print(f"      {l[:150]}")
    print(f"\n{len(rows) - fails}/{len(rows)} conditions pass")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
