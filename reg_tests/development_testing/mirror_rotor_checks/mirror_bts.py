#!/usr/bin/env python3
"""Create a mirror-image (across the XZ plane, i.e. Y -> -Y) copy of a TurbSim
.bts full-field wind file.

Usage:
    python3 mirror_bts.py <input.bts> <output.bts>
"""
import sys

import bts_io


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    print(f"Reading {in_path} ...")
    bts = bts_io.read_bts(in_path)
    print(f"  FileID={bts.FileID} NYGrids={bts.NYGrids} NZGrids={bts.NZGrids} "
          f"NTGrids={bts.NTGrids} NSteps={bts.NSteps} dy={bts.dy} dz={bts.dz} dt={bts.dt}")
    print(f"  DescStr: {bts.DescStr!r}")

    print("Mirroring across XZ plane (Y -> -Y; reversing Y grid index, flipping V sign)...")
    mirrored = bts_io.mirror_xz(bts, note="Mirrored across XZ plane (Y -> -Y; V-component sign flipped) for comparison testing.")

    print(f"Writing {out_path} ...")
    bts_io.write_bts(out_path, mirrored)
    print("Done.")


if __name__ == "__main__":
    main()
