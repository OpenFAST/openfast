#!/usr/bin/env python3
"""Run `<exe> -CheckInput <deck>` on a copy of an r-test case (optionally corrupted first)
and assert exit code + .verify.yaml contents. No baseline comparison needed.

Fixture handling: r-test glue-codes/openfast cases are not self-contained -- they
reference sibling case directories by relative path (e.g. AOC_WSt's ElastoDyn deck
points at "../AOC/AOC_Blade.dat"; the 5MW cases point at "../5MW_Baseline/...").
Copying only the case directory therefore breaks every module that follows such a
reference. After copying the case, this script scans the copied input files for
"../<name>/" references and copies each referenced sibling directory straight from
the pristine r-test tree, placing it as a sibling of the copied case directory --
the same relative layout r-test itself has -- so the relative paths resolve exactly
as they do "at home". Nothing is ever written back into the r-test submodule.
"""
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Matches a relative reference to a sibling case directory, e.g. "../AOC/AOC_Blade.dat"
# or "../5MW_Baseline/ServoData/DISCON.dll" -> captures "AOC" / "5MW_Baseline".
SIBLING_RE = re.compile(r'\.\./([A-Za-z0-9_.\-]+)/')
# Extensions worth scanning for sibling-directory references.
SIBLING_SCAN_GLOBS = ("*.fst", "*.dat", "*.inp", "*.txt")


def copy_case_with_siblings(src: Path, container: Path) -> Path:
    """Copy the r-test case dir `src` into container/case, then discover and copy
    every sibling case directory it references (via '../<name>/' relative paths)
    from src's parent into container/<name> -- i.e. as a sibling of container/case,
    matching the layout the case expects in r-test itself. Returns the case dir."""
    work = container / "case"
    shutil.copytree(src, work)

    siblings = set()
    for pattern in SIBLING_SCAN_GLOBS:
        for f in work.rglob(pattern):
            try:
                text = f.read_text(errors="replace")
            except OSError:
                continue
            siblings.update(SIBLING_RE.findall(text))

    for name in sorted(siblings):
        sdir = src.parent / name
        dst = container / name
        if sdir.is_dir() and not dst.exists():
            shutil.copytree(sdir, dst)

    return work


def corrupt(case_dir: Path, spec: str):
    """spec: '<fileglob>::<regex>::<replacement>' applied once, first matching file."""
    fileglob, pattern, repl = spec.split("::", 2)
    for f in sorted(case_dir.rglob(fileglob)):
        text = f.read_text(errors="replace")
        new, n = re.subn(pattern, repl, text, count=1, flags=re.M)
        if n:
            f.write_text(new)
            return f
    sys.exit(f"corruption spec matched no file: {spec}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("executable"); ap.add_argument("source_case"); ap.add_argument("build_dir")
    ap.add_argument("--corrupt", action="append", default=[])
    ap.add_argument("--expect-exit", type=int, required=True)
    ap.add_argument("--expect-status", choices=["passed", "failed"], required=True)
    ap.add_argument("--expect-min-fatals", type=int, default=0)
    ap.add_argument("--expect-component-failed", action="append", default=[],
                    help="component name that must have status: failed")
    a = ap.parse_args()

    src = Path(a.source_case); container = Path(a.build_dir)
    if container.exists(): shutil.rmtree(container)
    container.mkdir(parents=True)
    work = copy_case_with_siblings(src, container)
    for spec in a.corrupt: corrupt(work, spec)

    fst = next(work.glob("*.fst"))
    r = subprocess.run([a.executable, "-CheckInput", fst.name], cwd=work,
                       capture_output=True, text=True, timeout=600)
    print(r.stdout[-4000:]); print(r.stderr[-2000:], file=sys.stderr)

    ok = True
    if r.returncode != a.expect_exit:
        print(f"FAIL: exit {r.returncode} != {a.expect_exit}"); ok = False
    if "INPUT CHECK" not in r.stdout:
        print("FAIL: stdout has no INPUT CHECK summary"); ok = False

    reports = list(work.glob("*.verify.yaml"))
    if not reports:
        print("FAIL: no .verify.yaml written"); sys.exit(1)
    y = reports[0].read_text()

    m = re.search(r"^overall_status:\s*(\w+)", y, re.M)
    if not m:
        print("FAIL: no trailing overall_status block (crash-liveness contract)"); ok = False
    elif m.group(1) != a.expect_status:
        print(f"FAIL: overall_status {m.group(1)} != {a.expect_status}"); ok = False

    nfatal = len(re.findall(r"severity: (?:fatal|error)", y))
    if nfatal < a.expect_min_fatals:
        print(f"FAIL: {nfatal} fatal/error messages < {a.expect_min_fatals}"); ok = False

    for comp in a.expect_component_failed:
        block = re.search(rf"- name: {re.escape(comp)}\n(?:    .*\n)*?    status: (\w+)", y)
        if not block or block.group(1) != "failed":
            got = block.group(1) if block else "absent"
            print(f"FAIL: component {comp} status {got} != failed"); ok = False

    # console/file parity: every yaml message must appear on stdout.
    # The console writer hard-wraps long messages (sometimes mid-word), so compare
    # against a whitespace-normalized stdout, but keep summary-entry boundaries:
    # entries start with a "[error]/[warn]/[info]" tag, so join wrapped lines onto
    # their preceding tagged line before comparing.
    entries = []
    for line in r.stdout.splitlines():
        s = line.strip()
        if re.match(r"\[(error|warn|info)\]", s) or not entries:
            entries.append(s)
        else:
            entries[-1] += " " + s
    def squash(t): return re.sub(r"\s+", "", t)
    squashed_entries = [squash(e) for e in entries]
    for text in re.findall(r'text: "(.*)"', y):
        if text and not any(squash(text[:60]) in e for e in squashed_entries):
            print(f"FAIL: yaml message missing from stdout summary: {text[:60]}"); ok = False

    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
