"""Fix the BasicHAWTFormat comment, and add the missing MirrorRotor input.

The comment claimed the basic inputs were the next 7 lines. There are 8, from
BaseOriginInit to Twr2Shft, and MirrorRotor now sits between the flag and them,
so a positional description is doubly wrong. Describe the range by name instead.

The two example files under docs/ never gained the MirrorRotor input that the
driver now requires, so anyone copying one would have written a file the driver
rejects.
"""
import os
import re
import pathlib

ROOT = pathlib.Path(os.environ.get("OPENFAST_REPO") or
                     pathlib.Path(__file__).resolve().parents[2])

OLD = "True: next 7 lines are basic inputs, False:"
NEW = "True: basic HAWT inputs BaseOriginInit to Twr2Shft, False:"

MIRROR = ("False           MirrorRotor({n}) - Flag indicating the rotor rotation "
          "direction is mirrored (counter-clockwise viewed from upwind)")

EXAMPLES = [ROOT / "docs/source/user/aerodyn/examples/ad_driver_example.dvr",
            ROOT / "docs/source/user/aerodyn/examples/ad_driver_multiple.dvr"]

added = 0
for path in EXAMPLES:
    lines = path.read_text().split("\n")
    out = []
    for line in lines:
        out.append(line)
        m = re.match(r"\s*\S+\s+BasicHAWTFormat\((\d+)\)", line)
        if m:
            out.append(MIRROR.format(n=m.group(1)))
            added += 1
    path.write_text("\n".join(out))
print(f"added {added} MirrorRotor lines")

fixed = 0
for path in ROOT.rglob("*"):
    if not path.is_file() or "/build" in str(path):
        continue
    if path.suffix not in (".dvr", ".rst"):
        continue
    text = path.read_text(errors="replace")
    if OLD not in text:
        continue
    path.write_text(text.replace(OLD, NEW))
    fixed += 1
    print(f"  {path.relative_to(ROOT)}")
print(f"reworded {fixed} files")
