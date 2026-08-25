#!/usr/bin/env python3
"""Clone an r-test OpenFAST case and run it clockwise and mirrored.

Unlike openfast_mirror.py this changes *nothing* except the MirrorRotor flag, so
it is the right tool once ServoDyn is in the loop: the controller, its DLL and
every input file stay exactly as the clockwise case ships them.

Usage: openfast_case_mirror.py <CaseName> [--tmax N] [--dump CH1,CH2] [--tol X]
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
# Cases with a Bladed-style controller reference a DISCON library that is built, not
# stored, so it exists only under the build tree. Prefer that copy of a referenced
# sibling directory when it is there; ctest stages it from the r-test source.
STAGED = f"{REPO}/build-docker-double/reg_tests/glue-codes/openfast"
EXE = f"{REPO}/build-docker-double-debug/glue-codes/openfast/openfast"
WORK = os.path.join(os.environ.get("TMPDIR", "/tmp"), "of_case_mirror")

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, f"{REPO}/reg_tests/lib")
from compare_mirror import classify, blade_permutation  # noqa: E402
import fast_io  # noqa: E402

DIAG = {"Time", "ConvError", "ConvIter", "NumUJac"}


def set_var(text, var, value):
    pat = re.compile(r"^(\s*)(\S+)(\s+" + re.escape(var) + r"\s)", re.MULTILINE)
    new, n = pat.subn(lambda m: f"{m.group(1)}{value}{m.group(3)}", text, count=1)
    if n != 1:
        raise KeyError(f"could not set {var!r}")
    return new


def add_partner_channels(path, nblades):
    """Ensure every blade-indexed output has its mirror partner in the OutList.

    A case that only prints blades 1 and 3 gives the mirrored run's blade 3 nothing
    to be compared against, since it corresponds to the clockwise run's blade 2.
    """
    txt = open(path, errors="replace").read()
    lines = txt.split("\n")
    cand = [i for i, l in enumerate(lines) if "outlist" in l.lower() and '"' not in l]
    if not cand:
        return
    # These files carry a second, usually empty, nodal-output OutList at the end.
    start = cand[0]
    end = next((i for i in range(start, len(lines)) if lines[i].strip().upper().startswith("END")), None)
    if end is None:
        return

    have, order = set(), []
    for i in range(start + 1, end):
        for ch in re.findall(r'"([^"]+)"', lines[i]):
            for c in ch.split(","):
                c = c.strip()
                if c:
                    have.add(c)
                    order.append(c)

    add = []
    for c in order:
        mate = blade_permutation(c, nblades)
        if mate != c and mate not in have:
            have.add(mate)
            add.append(mate)
    if not add:
        return

    lines[end:end] = [f'"{c}"' for c in add]
    open(path, "w").write("\n".join(lines))


def negate_vars(d, names):
    """Negate asymmetric scalar inputs in the mirrored clone.

    Quantities such as an initial nacelle yaw or a wind direction are ordinary
    asymmetric inputs, not rotor-convention ones, so mirroring the turbine does not
    mirror them. A symmetry test has to mirror them by hand.
    """
    for f in sorted(os.listdir(d)):
        if not f.lower().endswith((".dat", ".fst")):
            continue
        p = os.path.join(d, f)
        txt = open(p, errors="replace").read()
        orig = txt
        for var in names:
            pat = re.compile(r"^(\s*)(-?[\d.]+(?:[eEdD][-+]?\d+)?)(\s+" + re.escape(var) + r"\s)",
                             re.MULTILINE)
            txt = pat.sub(lambda m: f"{m.group(1)}{-float(m.group(2)):g}{m.group(3)}", txt)
        if txt != orig:
            open(p, "w").write(txt)


def mirror_series(d, names):
    """Apply the mirror to a prescribed force/moment series (t, FX..MZ).

    Under S = diag(1, -1, 1) the true vector components FY flips, and for the moment
    pseudovector MX and MZ flip instead.
    """
    sign = [1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
    for f in names:
        p = os.path.join(d, f)
        if not os.path.exists(p):
            sys.exit(f"no such series file: {f}")
        out = []
        for line in open(p, errors="replace").read().split("\n"):
            if line.lstrip().startswith(("#", "!")):
                out.append(line)
                continue
            # Rows may carry a trailing inline comment.
            cut = min((i for i in (line.find("#"), line.find("!")) if i >= 0), default=len(line))
            body, tail = line[:cut], line[cut:]
            parts = body.split()
            if len(parts) != len(sign):
                out.append(line)
                continue
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                out.append(line)
                continue
            out.append("  ".join(f"{s*v:.8E}" for s, v in zip(sign, vals)) + ("  " + tail if tail else ""))
        open(p, "w").write("\n".join(out))


def build(case, mirrored, tmax, negate=(), series=(), reuse=False):
    name = "mir" if mirrored else "cw"
    d = os.path.join(WORK, f"{case}_{name}")
    if reuse and os.path.exists(os.path.join(d, f"{case}.out")):
        return os.path.join(d, f"{case}.out")
    os.makedirs(WORK, exist_ok=True)
    # Sibling directories are referenced by relative path from the case folder.
    for entry in os.listdir(RTEST):
        link = os.path.join(WORK, entry)
        if not os.path.exists(link):
            staged = os.path.join(STAGED, entry)
            os.symlink(staged if os.path.isdir(staged) else os.path.join(RTEST, entry), link)

    d = os.path.join(WORK, f"{case}_{name}")
    shutil.rmtree(d, ignore_errors=True)
    shutil.copytree(os.path.join(RTEST, case), d)

    fst = os.path.join(d, f"{case}.fst")
    txt = open(fst).read()
    txt = set_var(txt, "MirrorRotor", "True" if mirrored else "False")
    txt = set_var(txt, "OutFileFmt", "1")
    if tmax is not None:
        txt = set_var(txt, "TMax", str(tmax))
    open(fst, "w").write(txt)

    # Give every blade-indexed channel a partner to be compared against.
    for f in os.listdir(d):
        if not f.lower().endswith(".dat"):
            continue
        p = os.path.join(d, f)
        body = open(p, errors="replace").read()
        head = body[:400].upper()
        if "ELASTODYN" in head or "AERODYN" in head:
            m = re.search(r"^\s*(\d+)\s+NumBl\b", body, re.MULTILINE)
            add_partner_channels(p, int(m.group(1)) if m else 3)

    if mirrored and negate:
        negate_vars(d, negate)
    if mirrored and series:
        mirror_series(d, series)

    p = subprocess.run([EXE, f"{case}.fst"], cwd=d, capture_output=True, text=True)
    out = os.path.join(d, f"{case}.out")
    if not os.path.exists(out):
        tail = [l for l in (p.stdout + p.stderr).split("\n") if l.strip()][-10:]
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
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    case = args[0]
    tmax = float(sys.argv[sys.argv.index("--tmax") + 1]) if "--tmax" in sys.argv else None
    tol = float(sys.argv[sys.argv.index("--tol") + 1]) if "--tol" in sys.argv else 1e-6
    reuse = "--reuse" in sys.argv
    negate = sys.argv[sys.argv.index("--negate") + 1].split(",") if "--negate" in sys.argv else ()
    series = sys.argv[sys.argv.index("--mirror-series") + 1].split(",") if "--mirror-series" in sys.argv else ()

    # Components that are not blades also permute under the mirror, but which ones pair up
    # depends on the layout rather than on blade order -- the OC4 semi's mooring line 2
    # lies on the mirror plane, so lines 1 and 3 swap. State those pairings explicitly.
    swap = {}
    if "--swap" in sys.argv:
        for pair in sys.argv[sys.argv.index("--swap") + 1].split(","):
            a, _, b = pair.partition("=")
            swap[a.strip()], swap[b.strip()] = b.strip(), a.strip()

    cw = build(case, False, tmax, negate, series, reuse)
    mir = build(case, True, tmax, negate, series, reuse)
    n1, a = load(cw)
    n2, b = load(mir)
    if n1 != n2:
        sys.exit("channel lists differ")
    idx = {k: v for v, k in enumerate(n1)}
    sl = slice(int(0.25 * len(a)), None)
    phys = [idx[c] for c in n1 if c not in DIAG]
    floor = 1e-9 * np.abs(a[sl][:, phys]).max()
    nb = int(sys.argv[sys.argv.index("--nblades") + 1]) if "--nblades" in sys.argv else 3

    if "--dump" in sys.argv:
        rows = list(range(0, len(a), max(1, len(a) // 14)))
        for ch in sys.argv[sys.argv.index("--dump") + 1].split(","):
            print(f"\n{ch}:")
            print(f"{'t':>8s} {'CW':>14s} {'MIRROR':>14s}")
            for r in rows:
                print(f"{a[r,0]:8.3f} {a[r,idx[ch]]:14.6g} {b[r,idx[ch]]:14.6g}")
        return

    groups = {}
    for ch in n1:
        if ch in DIAG:
            continue
        mate = swap.get(ch) or blade_permutation(ch, nb)
        if mate not in idx:
            mate = ch
        got, rel = classify(a[sl, idx[ch]], b[sl, idx[mate]], tol, floor)
        groups.setdefault(got, []).append((ch, rel))

    for k in ("S", "F", "A", "negligible", "?"):
        if groups.get(k):
            print(f"{k:11s} ({len(groups[k]):3d}): {' '.join(c for c, _ in groups[k])}")
    bad = groups.get("?", [])
    if bad:
        print(f"\nFAIL: {len(bad)} channel(s) unresolved")
        for c, r in bad[:10]:
            print(f"   {c:14s} rel={r:.3e}")
        sys.exit(1)
    print("\nPASS: every channel resolves")


if __name__ == "__main__":
    main()
