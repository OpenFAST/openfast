#!/usr/bin/env python3
"""Run OpenFAST regression cases and compare against their stored baselines.

Runs each case in a scratch copy of the r-test tree so the committed reference
outputs are never overwritten.

Usage: openfast_regression.py <case> [<case> ...]
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
RTEST = f"{REPO}/reg_tests/r-test/glue-codes/openfast"
EXE = f"{REPO}/build-docker-double-debug/glue-codes/openfast/openfast"
sys.path.insert(0, f"{REPO}/reg_tests/lib")
import fast_io  # noqa: E402

DETAIL = "--detail" in sys.argv


def read(path):
    data, info, _ = fast_io.load_output(path)
    return info["attribute_names"], np.asarray(data)


def run_case(case):
    src = os.path.join(RTEST, case)
    fst = [f for f in os.listdir(src) if f.endswith(".fst")]
    if not fst:
        return case, "SKIP", "no .fst"
    fst = fst[0]
    stem = fst[:-4]
    ref = None
    for cand in (f"{stem}.outb", f"{stem}.out"):
        if os.path.exists(os.path.join(src, cand)):
            ref = cand
            break
    if ref is None:
        return case, "SKIP", "no reference output"

    with tempfile.TemporaryDirectory(dir=os.environ.get("TMPDIR")) as tmp:
        for entry in os.listdir(RTEST):
            if entry != case:
                os.symlink(os.path.join(RTEST, entry), os.path.join(tmp, entry))
        work = os.path.join(tmp, case)
        shutil.copytree(src, work)

        p = subprocess.run([EXE, fst], cwd=work, capture_output=True, text=True)
        out = os.path.join(work, ref)
        if not os.path.exists(out):
            tail = [l for l in (p.stdout + p.stderr).split("\n") if l.strip()][-3:]
            return case, "RUNFAIL", " | ".join(tail)

        rn, rd = read(os.path.join(src, ref))
        tn, td = read(out)
        if rn != tn:
            return case, "FAIL", "channel list changed"
        if rd.shape != td.shape:
            return case, "FAIL", f"shape {rd.shape} vs {td.shape}"

        scale = np.maximum(np.abs(rd).max(axis=0), 1e-12)
        col = (np.abs(td - rd) / scale).max(axis=0)
        rel = col.max()
        if DETAIL and rel > 0.0:
            order = np.argsort(col)[::-1]
            print(f"    {case}: {int((col > 0).sum())} of {len(col)} channels differ")
            for k in order[:10]:
                if col[k] > 0:
                    print(f"      {tn[k]:16s} rel={col[k]:.3e}")
        if rel == 0.0:
            return case, "IDENTICAL", "rel=0"
        if rel < 1e-10:
            return case, "OK", f"rel={rel:.2e}"
        return case, "FAIL", f"rel={rel:.2e} worst={tn[int(np.argmax(col))]}"


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    fails = 0
    for case in args:
        name, status, detail = run_case(case)
        if status in ("FAIL", "RUNFAIL"):
            fails += 1
        print(f"  {name:34s} {status:10s} {detail}", flush=True)
    print(f"\n{len(args) - fails}/{len(args)} cases match")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
