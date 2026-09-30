"""Check no output channel name appears twice in a tab.

A name listed both as a primary and as somebody else's alias is exactly the
condition the alias split exists to remove: the workbook would claim two
channels are the same while the code gives them opposite signs.
"""
import collections
import subprocess
import sys

TABS = ["ElastoDyn", "SimpleElastoDyn"]

for tab in TABS:
    out = subprocess.run([sys.executable, "dump_outlist_xlsx.py", tab],
                         capture_output=True, text=True).stdout
    seen = collections.defaultdict(list)
    for line in out.split("\n"):
        if not line.strip():
            continue
        row, _, rest = line.strip().partition(" ")
        for field in rest.split("|"):
            field = field.strip()
            if field.startswith("B:"):
                seen[field[2:].strip()].append(row)
            elif field.startswith("C:"):
                for name in field[2:].split(","):
                    if name.strip():
                        seen[name.strip()].append(row)

    dupes = {k: v for k, v in seen.items() if len(v) > 1 and k not in ("Name", "Other Name(s)")}
    print(f"{tab}: {len(seen)} names")
    for name, rows in sorted(dupes.items()):
        print(f"  DUPLICATE {name!r} in rows {', '.join(rows)}")
    if not dupes:
        print("  no duplicates")
