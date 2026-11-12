"""Check no output channel name appears twice in a tab.

A name listed both as a primary and as somebody else's alias is exactly the
condition the alias split exists to remove: the workbook would claim two
channels are the same while the code gives them opposite signs.
"""
import collections
import os
import subprocess
import sys

TABS = ["ElastoDyn", "SimpleElastoDyn"]
HERE = os.path.dirname(os.path.abspath(__file__))

failed = False
for tab in TABS:
    # The helper is resolved next to this script, not against the caller's working
    # directory, and its failure is fatal: a silently empty dump used to report
    # "0 names / no duplicates" from any other directory.
    proc = subprocess.run([sys.executable, os.path.join(HERE, "dump_outlist_xlsx.py"), tab],
                          capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(proc.stderr)
        sys.exit(f"dump_outlist_xlsx.py failed for tab {tab!r}")
    out = proc.stdout
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
    if dupes:
        failed = True
    else:
        print("  no duplicates")

if failed:
    sys.exit(1)
