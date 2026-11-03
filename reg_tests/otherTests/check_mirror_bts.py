"""Check that the mirrored .bts really is the y-reflection of the original.

Reads both boxes back through the same reader and requires u and w to be the
reverse of each other along y, and v to be the reverse and sign-flipped. Reading
the written file back matters: the format stores int16 with a per-component
scale and offset, so the round trip is where a mistake would show up.
"""
import os
import sys

import numpy as np

import bts_io

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

BASE = os.path.join(REPO_ROOT, "reg_tests", "r-test", "glue-codes",
                     "openfast", "5MW_Baseline", "Wind")
orig = bts_io.read_bts(f"{BASE}/90m_12mps_twr.bts")
mir = bts_io.read_bts(f"{BASE}/90m_12mps_twr_MirrorY.bts")

for name in ("FileID", "NYGrids", "NZGrids", "NTGrids", "NSteps",
             "dy", "dz", "dt", "MeanWS", "RefHeight", "GridBase"):
    a, b = getattr(orig, name), getattr(mir, name)
    if a != b:
        sys.exit(f"{name} differs: {a} vs {b}")

o = orig.decode_grid()   # (NSteps, NZGrids, NYGrids, 3)
m = mir.decode_grid()
print(f"grid {o.shape}  (time, z, y, component)")

sign = (1.0, -1.0, 1.0)
worst = 0.0
for c, s in enumerate(sign):
    expect = s * o[:, :, ::-1, c]
    err = np.abs(m[:, :, :, c] - expect).max() / max(np.abs(o[..., c]).max(), 1e-12)
    worst = max(worst, err)
    print(f"  component {'uvw'[c]}: max rel error {err:.3e}")

# The tower line sits at y = 0, so it maps onto itself and only v changes sign.
ot, mt = orig.decode_tower(), mir.decode_tower()
if ot is not None:
    for c, s in enumerate(sign):
        err = np.abs(mt[..., c] - s * ot[..., c]).max() / max(np.abs(ot[..., c]).max(), 1e-12)
        worst = max(worst, err)
        print(f"  tower {'uvw'[c]}: max rel error {err:.3e}")

# One int16 quantum is 1/32768 of the stored range, and the two files are quantised
# independently, so a couple of quanta is the floor.
ok = worst < 1e-3
print("PASS" if ok else "FAIL", f"worst={worst:.3e}")
sys.exit(0 if ok else 1)
