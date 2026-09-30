#!/usr/bin/env python3
"""Check the ElastoDyn output-parameter tables stay consistent.

ValidParamAry, ParamIndxAry and ParamUnitsAry are three parallel Fortran array
literals indexed by position. Nothing in the build checks they line up, so a
misalignment silently reports the wrong channel. This parses all three straight
out of the source and verifies:

  - all three have the declared length, and the same length as each other
  - ValidParamAry is sorted alphabetically (SetOutParam binary-searches it)
  - no duplicate channel names
  - every ParamIndxAry symbol is a declared AllOuts index within MaxOutPts

With channel names as arguments it also prints what each currently maps to,
which is the safe way to retarget an alias.

Usage: check_ed_outparams.py [CHANNEL ...]
"""
import os
import re
import sys

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ED = os.path.join(REPO_ROOT, "modules", "elastodyn", "src", "ElastoDyn.f90")
IO = os.path.join(REPO_ROOT, "modules", "elastodyn", "src", "ElastoDyn_IO.f90")


def grab_array(text, name):
    """Return (declared_length, [entries]) for a Fortran array literal."""
    m = re.search(rf"{name}\((\d+)\)\s*=\s*\(/", text)
    if not m:
        sys.exit(f"could not find {name}")
    declared = int(m.group(1))
    start = m.end()
    end = text.index("/)", start)
    body = text[start:end]
    # Strip Fortran trailing comments and continuations before splitting, or the
    # comment on the declaration line gets parsed as the first entry.
    body = "\n".join(l.split("!")[0] for l in body.split("\n"))
    body = re.sub(r"character\(ChanLen\)\s*::", "", body)
    body = body.replace("&", " ")
    items = [t.strip() for t in body.split(",")]
    items = [t for t in items if t]
    items = [t[1:-1].strip() if t.startswith('"') else t for t in items]
    return declared, items


def main():
    ed = open(ED).read()
    io = open(IO).read()

    indices = dict(re.findall(
        r"INTEGER\(IntKi\),\s*PARAMETER\s*::\s*(\w+)\s*=\s*(\d+)", io))
    indices = {k: int(v) for k, v in indices.items()}
    maxout = indices.get("MaxOutPts")

    n_names, names = grab_array(ed, "ValidParamAry")
    n_indx, indx = grab_array(ed, "ParamIndxAry")
    n_unit, units = grab_array(ed, "ParamUnitsAry")

    bad = 0
    print(f"MaxOutPts = {maxout}")
    for label, declared, got in (("ValidParamAry", n_names, len(names)),
                                 ("ParamIndxAry", n_indx, len(indx)),
                                 ("ParamUnitsAry", n_unit, len(units))):
        ok = declared == got
        bad += not ok
        print(f"{label:15s} declared {declared:5d}  parsed {got:5d}  "
              f"{'ok' if ok else 'MISMATCH'}")

    if not (len(names) == len(indx) == len(units)):
        print("ARRAYS ARE NOT THE SAME LENGTH - they are indexed in parallel")
        bad += 1

    unsorted = [(a, b) for a, b in zip(names, names[1:]) if a.strip() > b.strip()]
    if unsorted:
        bad += 1
        print(f"NOT ALPHABETICAL ({len(unsorted)}): first at {unsorted[0]}")
    else:
        print("ValidParamAry    alphabetical  ok")

    seen, dupes = set(), []
    for nm in names:
        s = nm.strip()
        if s in seen:
            dupes.append(s)
        seen.add(s)
    if dupes:
        bad += 1
        print(f"DUPLICATE NAMES: {dupes}")
    else:
        print("ValidParamAry    no duplicates ok")

    unknown = sorted({s for s in indx if s not in indices})
    if unknown:
        bad += 1
        print(f"UNKNOWN INDEX SYMBOLS: {unknown[:10]}")
    else:
        over = sorted({s for s in indx if indices[s] > maxout})
        if over:
            bad += 1
            print(f"INDEX SYMBOLS BEYOND MaxOutPts: {over}")
        else:
            print("ParamIndxAry     all symbols known and within MaxOutPts  ok")

    wanted = [a.upper() for a in sys.argv[1:]]
    if wanted:
        print("\nrequested channels:")
        lookup = {nm.strip(): i for i, nm in enumerate(names)}
        for w in wanted:
            if w not in lookup:
                print(f"  {w:12s} NOT FOUND")
                continue
            i = lookup[w]
            sym = indx[i]
            print(f"  {w:12s} pos {i:5d} -> {sym:12s} "
                  f"(AllOuts {indices.get(sym, '?')})  units {units[i].strip()}")

    print("\nOK" if not bad else f"\n{bad} problem(s)")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
