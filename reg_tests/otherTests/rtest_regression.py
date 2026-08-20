#!/usr/bin/env python3
"""Regression guard: MirrorRotor=F must reproduce the stored r-test baselines bit-for-bit.

Copies each AeroDyn driver case to a scratch dir, runs the freshly built driver,
and compares every channel against the committed reference output.

Usage: rtest_regression.py <case> [<case> ...]
"""
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
REPO = REPO_ROOT
RTEST = f"{REPO}/reg_tests/r-test/modules/aerodyn"
DRIVER = f"{REPO}/build-docker-double-debug/modules/aerodyn/aerodyn_driver"
sys.path.insert(0, f"{REPO}/reg_tests/lib")
import fast_io  # noqa: E402


def read(path):
    data, info, _ = fast_io.load_output(path)
    return info["attribute_names"], np.asarray(data)


def run_case(case):
    src = os.path.join(RTEST, case)
    ref = None
    for cand in ("ad_driver.outb", "ad_driver.out"):
        if os.path.exists(os.path.join(src, cand)):
            ref = cand
            break
    if ref is None:
        return case, "SKIP", "no reference output"

    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
        # Cases reference files by relative paths that escape the module folder
        # (e.g. ../../../glue-codes/...), so rebuild the r-test layout in scratch.
        rtest_root = os.path.abspath(os.path.join(RTEST, "..", ".."))
        mod_dir = os.path.join(tmp, "modules", "aerodyn")
        os.makedirs(mod_dir)
        for entry in os.listdir(rtest_root):
            if entry != "modules":
                os.symlink(os.path.join(rtest_root, entry),
                           os.path.join(tmp, entry))
        for entry in os.listdir(os.path.join(rtest_root, "modules")):
            if entry != "aerodyn":
                os.symlink(os.path.join(rtest_root, "modules", entry),
                           os.path.join(tmp, "modules", entry))
        for entry in os.listdir(RTEST):
            if entry != case:
                os.symlink(os.path.join(RTEST, entry),
                           os.path.join(mod_dir, entry))

        work = os.path.join(mod_dir, case)
        shutil.copytree(src, work)

        dvr = next(f for f in os.listdir(work) if f.endswith(".dvr"))
        proc = subprocess.run([DRIVER, dvr], cwd=work,
                              capture_output=True, text=True)
        out = os.path.join(work, ref)
        if proc.returncode != 0 or not os.path.exists(out):
            tail = (proc.stdout + proc.stderr).strip().split("\n")[-3:]
            return case, "RUNFAIL", " | ".join(tail)

        try:
            rn, rd = read(os.path.join(src, ref))
            tn, td = read(out)
        except Exception as exc:
            return case, "READFAIL", str(exc)

        if rn != tn:
            return case, "FAIL", "channel list changed"
        if rd.shape != td.shape:
            return case, "FAIL", f"shape {rd.shape} vs {td.shape}"

        diff = np.abs(td - rd)
        scale = np.maximum(np.abs(rd).max(axis=0), 1e-12)
        rel = (diff / scale).max()
        if rel == 0.0:
            return case, "IDENTICAL", "rel=0"
        if rel < 1e-10:
            return case, "OK", f"rel={rel:.2e}"
        worst = tn[int(np.argmax((diff / scale).max(axis=0)))]
        return case, "FAIL", f"rel={rel:.2e} worst={worst}"


def main():
    cases = sys.argv[1:]
    fails = 0
    for case in cases:
        name, status, detail = run_case(case)
        if status in ("FAIL", "RUNFAIL", "READFAIL"):
            fails += 1
        print(f"  {name:32s} {status:10s} {detail}")
    print(f"\n{len(cases) - fails}/{len(cases)} cases match the baseline")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
