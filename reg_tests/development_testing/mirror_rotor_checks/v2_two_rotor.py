"""V2: one clockwise rotor and one mirrored rotor inside a single driver run.

Turbine 2 carries MirrorRotor, turbine 1 does not, and the two are otherwise the
same machine in the same uniform inflow. Their loads must therefore be exact
mirror images of each other, with no run-to-run difference of any kind to hide
behind: same solver, same time steps, same wind.

Blade spacing follows the rotation sense, so a mirrored blade compares directly
against the clockwise blade of the same number.
"""
import os
import sys

import numpy as np

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

sys.path.insert(0, os.path.join(REPO_ROOT, "reg_tests", "lib"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fast_io  # noqa: E402
from compare_mirror import expected_class  # noqa: E402

_args = [a for i, a in enumerate(sys.argv[1:], 1)
         if not a.startswith("--") and sys.argv[i - 1] != "--tol"]
WORK = _args[0] if _args else os.path.join(
    os.environ.get("TMPDIR", "/tmp"), "v2test", "mir")
TOL = float(sys.argv[sys.argv.index("--tol") + 1]) if "--tol" in sys.argv else 1e-6


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
wrong = []
for ch in n1:
    if ch == "Time":
        continue
    x, y = a[sl, idx[ch]], b[sl, idx[ch]]
    # Both turbines must sit below the floor before the channel is skipped: a
    # channel quiet on turbine 1 but live on turbine 2 is exactly the one-sided
    # asymmetry this comparison exists to catch.
    peak = max(np.abs(x).max(), np.abs(y).max())
    if peak < floor:
        groups.setdefault("negligible", []).append((ch, 0.0))
        continue
    rs = np.abs(y - x).max() / peak
    rf = np.abs(y + x).max() / peak
    if rs < TOL:
        got, r = "S", rs
    elif rf < TOL:
        got, r = "F", rf
    else:
        got, r = "?", min(rs, rf)
    groups.setdefault(got, []).append((ch, r))
    # Channels with a known expected class must land in it: two clockwise
    # turbines agree on every channel and must not pass as a mirror pair.
    exp = expected_class(ch)
    if exp is not None and got not in ("negligible", exp):
        wrong.append((ch, exp, got))

for k, label in (("S", "same"), ("F", "flipped"), ("negligible", "negligible")):
    if groups.get(k):
        print(f"{label:11s} {len(groups[k]):4d}")
bad = groups.get("?", [])
print(f"unresolved  {len(bad):4d}")
for c, r in sorted(bad, key=lambda t: -t[1])[:25]:
    print(f"   {c:16s} {r:.3e}")
if wrong:
    print("wrong class: "
          + " ".join(f"{c}(expected {e}, got {g})" for c, e, g in wrong))
sys.exit(1 if bad or wrong else 0)
