"""V2: one clockwise rotor and one mirrored rotor inside a single driver run.

Turbine 2 carries MirrorRotor, turbine 1 does not, and the two are otherwise the
same machine in the same uniform inflow. Their loads must therefore be exact
mirror images of each other, with no run-to-run difference of any kind to hide
behind: same solver, same time steps, same wind.

Mirroring reverses the order in which the blades sweep, so blade 1 keeps its
place and blades 2 and 3 exchange.
"""
import os
import re
import sys

import numpy as np

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

sys.path.insert(0, os.path.join(REPO_ROOT, "reg_tests", "lib"))
import fast_io  # noqa: E402

WORK = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.environ.get("TMPDIR", "/tmp"), "v2test", "mir")
TOL = float(sys.argv[sys.argv.index("--tol") + 1]) if "--tol" in sys.argv else 1e-6

# Driver channels are AB<blade>N<node><quantity>; the module form is B<blade>N<node>...
AB = re.compile(r"^(A?B)(\d+)(N.*)$")


def swap(name, nblades=3):
    m = AB.match(name)
    if not m:
        return name
    k = int(m.group(2))
    if not 1 <= k <= nblades:
        return name
    return f"{m.group(1)}{1 if k == 1 else nblades - k + 2}{m.group(3)}"


def load(path):
    d, info, _ = fast_io.load_output(path)
    return info["attribute_names"], np.asarray(d)


def find(turbine):
    for f in sorted(os.listdir(WORK)):
        if f".T{turbine}." in f and f.endswith((".outb", ".out")):
            return os.path.join(WORK, f)
    sys.exit(f"no output file for turbine {turbine} in {WORK}")


n1, a = load(find(1))
n2, b = load(find(2))
if n1 != n2:
    sys.exit("turbine 1 and turbine 2 report different channels")

idx = {k: v for v, k in enumerate(n1)}
sl = slice(int(0.25 * len(a)), None)
phys = [idx[c] for c in n1 if c != "Time"]
floor = 1e-9 * np.abs(a[sl][:, phys]).max()

groups = {}
for ch in n1:
    if ch == "Time":
        continue
    mate = swap(ch)
    if mate not in idx:
        mate = ch
    x, y = a[sl, idx[ch]], b[sl, idx[mate]]
    peak = np.abs(x).max()
    if peak < floor:
        groups.setdefault("negligible", []).append((ch, 0.0))
        continue
    rs = np.abs(y - x).max() / peak
    rf = np.abs(y + x).max() / peak
    if rs < TOL:
        groups.setdefault("S", []).append((ch, rs))
    elif rf < TOL:
        groups.setdefault("F", []).append((ch, rf))
    else:
        groups.setdefault("?", []).append((ch, min(rs, rf)))

for k, label in (("S", "same"), ("F", "flipped"), ("negligible", "negligible")):
    if groups.get(k):
        print(f"{label:11s} {len(groups[k]):4d}")
bad = groups.get("?", [])
print(f"unresolved  {len(bad):4d}")
for c, r in sorted(bad, key=lambda t: -t[1])[:25]:
    print(f"   {c:16s} {r:.3e}")
sys.exit(1 if bad else 0)
