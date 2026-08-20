#!/usr/bin/env python3
"""Dump rows from OutListParameters.xlsx without needing openpyxl.

An xlsx is a zip of XML, so the sheet can be read directly. Used to give exact
row numbers for the OutListParameters edits.

Usage: dump_outlist_xlsx.py <TabName> [CHANNEL ...]
       dump_outlist_xlsx.py --tabs
"""
import os
import re
import sys
import zipfile
from xml.etree import ElementTree as ET

REPO_ROOT = os.environ.get("OPENFAST_REPO") or os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

XLSX = os.path.join(REPO_ROOT, "docs", "OtherSupporting", "OutListParameters.xlsx")
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
RNS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def load():
    z = zipfile.ZipFile(XLSX)
    strings = []
    if "xl/sharedStrings.xml" in z.namelist():
        root = ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall(f"{NS}si"):
            strings.append("".join(t.text or "" for t in si.iter(f"{NS}t")))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = {r.get("Id"): r.get("Target") for r in rels}
    sheets = {}
    for sh in wb.find(f"{NS}sheets"):
        rid = sh.get(f"{RNS}id")
        path = target[rid]
        if not path.startswith("xl/"):
            path = "xl/" + path.lstrip("/")
        sheets[sh.get("name")] = path
    return z, strings, sheets


def cells(z, strings, path):
    root = ET.fromstring(z.read(path))
    out = {}
    for row in root.iter(f"{NS}row"):
        r = int(row.get("r"))
        vals = {}
        for c in row.findall(f"{NS}c"):
            ref = c.get("r")
            col = re.match(r"([A-Z]+)", ref).group(1)
            v = c.find(f"{NS}v")
            isel = c.find(f"{NS}is")
            if c.get("t") == "s" and v is not None:
                vals[col] = strings[int(v.text)]
            elif isel is not None:
                vals[col] = "".join(t.text or "" for t in isel.iter(f"{NS}t"))
            elif v is not None:
                vals[col] = v.text
        if vals:
            out[r] = vals
    return out


def main():
    z, strings, sheets = load()
    if "--tabs" in sys.argv:
        for name in sheets:
            print(name)
        return

    tab = sys.argv[1]
    if tab not in sheets:
        sys.exit(f"tab {tab!r} not found; have: {', '.join(sheets)}")
    rows = cells(z, strings, sheets[tab])
    wanted = [a.upper() for a in sys.argv[2:]]

    for r in sorted(rows):
        vals = rows[r]
        joined = " | ".join(f"{k}:{v}" for k, v in sorted(vals.items()))
        if not wanted:
            print(f"{r:5d}  {joined[:200]}")
            continue
        for w in wanted:
            if any(str(v).strip().upper() == w for v in vals.values()):
                print(f"{r:5d}  {joined[:220]}")
                break


if __name__ == "__main__":
    main()
