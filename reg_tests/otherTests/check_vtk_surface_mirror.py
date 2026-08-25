#!/usr/bin/env python3
"""Check that the VTK blade surfaces mirror.

The blade surface written for ``VTK_type = 1`` is built from the airfoil
coordinate files, not from the mesh nodes, and it is the one part of the
mirrored geometry that no output channel and no node comparison can police.  A
mirrored rotor whose nodes are all in the right place can still be drawn with
its aerofoils the wrong way round, which is what happened: the coordinates were
carried across without ``RotDir`` and the surfaces sat 2.33 m from where they
belonged, while every channel and every node position agreed exactly.

The vertices are placed by ``MeshWrVTK_Ln2Surface`` as
``matmul(xyz, Orientation)``, so ``AirfoilCoords`` component 1 rides row 1 of the
node's direction cosine matrix and component 2 rides row 2.  Under ``R' = S R S``
the rows do not transform alike -- row 1 becomes ``S*row1`` and row 2 becomes
``-S*row2`` -- so the coordinates need ``RotDir`` on component 2 alone for the
two negations to cancel.  See ``AD_SetInitOut`` in ``modules/aerodyn/src/AeroDyn.f90``.

Three configurations are checked, because the surface is drawn three different
ways and each carries the sign separately:

``coords``
    The airfoil files supply coordinates, and ``AD_SetInitOut`` builds the
    section.  This is the normal case.
``nocoords``
    ``NumCoords = 0``, so ``BladeShape`` is never allocated and the glue code
    synthesises a generic S809 section instead.  The same fallback is reached
    when the airfoil files merely disagree on how many coordinates they have,
    which the code itself calls an overly restrictive check.
``noaero``
    ``CompAero = 0``, so the blades come from BeamDyn and the fallback draws a
    rectangle.  Nearly symmetric, so nearly immune -- but only nearly: the
    closing vertex of the rectangle sits at the midpoint of one long edge and
    its mirror lies on the opposite edge, which displaces one vertex per node.

Run it directly, or through ``run_guards.sh``::

    python3 check_vtk_surface_mirror.py
    python3 check_vtk_surface_mirror.py --only coords
    python3 check_vtk_surface_mirror.py --keep   # leave the runs in place

It clones the registered pair `5MW_Land_BD_DLL_WTurb` and
`5MW_Land_BD_DLL_WTurb_MirrorRotor` from the build tree, so the DISCON library
resolves, shortens them, turns on surface output and compares.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from render_vtk_pair import read_vtp, S  # noqa: E402


CW = "5MW_Land_BD_DLL_WTurb"
MR = "5MW_Land_BD_DLL_WTurb_MirrorRotor"

# Blade 1 lies on the mirror plane; blades 2 and 3 exchange.
PAIRS = {"1": "1", "2": "3", "3": "2"}

# The ASCII .vtp files carry five decimal places, so agreement cannot be
# demonstrated below 1e-5 m no matter how exact the arithmetic is.
TOL = 2.0e-5

# Each draws the surface a different way, and each carries the sign separately.
MODES = ("coords", "nocoords", "noaero")


def repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.environ.get("OPENFAST_REPO", os.path.abspath(os.path.join(here, "..", "..")))


def strip_coords(work, airfoil_src, aero_deck):
    """Copy the airfoil files with NumCoords switched off and point the deck at them.

    Setting NumCoords to 0 is what a user does when they have no profile geometry
    to hand, and it sends the surface down the generic-shape fallback.
    """
    dst = os.path.join(work, "Airfoils_nocoords")
    os.makedirs(dst, exist_ok=True)
    for fn in os.listdir(airfoil_src):
        if not fn.endswith(".dat"):
            continue
        with open(os.path.join(airfoil_src, fn)) as fh:
            text = fh.read()
        text = re.sub(r'^@"[^"]*"(\s+NumCoords)', r"0\1", text, flags=re.M)
        with open(os.path.join(dst, fn), "w") as fh:
            fh.write(text)
    with open(aero_deck) as fh:
        text = fh.read()
    text = text.replace('"../5MW_Baseline/Airfoils/', '"../Airfoils_nocoords/')
    with open(aero_deck, "w") as fh:
        fh.write(text)


def stage(repo, work, tmax, mode):
    """Clone both cases from the staged build tree and turn on surface VTK."""
    built = os.path.join(repo, "build-docker-double", "reg_tests", "glue-codes", "openfast")
    baseline = os.path.join(built, "5MW_Baseline")
    if not os.path.isdir(baseline):
        raise SystemExit(f"no staged 5MW_Baseline in {built}; run ctest once first")
    # The decks reach for ../5MW_Baseline, and the mirrored deck also reaches for
    # ../5MW_Land_BD_DLL_WTurb, so both must sit beside each other.
    os.symlink(baseline, os.path.join(work, "5MW_Baseline"))

    for case in (CW, MR):
        src = os.path.join(built, case)
        if not os.path.isdir(src):
            raise SystemExit(f"case not staged: {src}; run ctest -R {case} once first")
        dst = os.path.join(work, case)
        shutil.copytree(src, dst)
        shutil.rmtree(os.path.join(dst, "vtk"), ignore_errors=True)

        fst = os.path.join(dst, case + ".fst")
        with open(fst) as fh:
            text = fh.read()
        text = re.sub(r"^(\s*)\d+(\s+WrVTK\s)", r"\g<1>2\g<2>", text, flags=re.M)
        text = re.sub(r"^(\s*)\d+(\s+VTK_type\s)", r"\g<1>1\g<2>", text, flags=re.M)
        text = re.sub(r"^(\s*)[\d.]+(\s+TMax\s)", rf"\g<1>{tmax}\g<2>", text, flags=re.M)
        if mode == "noaero":
            # No AeroDyn at all, so the blades come from BeamDyn and the fallback
            # draws a rectangle.  ServoDyn goes too, since its controller expects
            # aerodynamic torque.
            text = re.sub(r"^(\s*)\d+(\s+CompAero\s)", r"\g<1>0\g<2>", text, flags=re.M)
            text = re.sub(r"^(\s*)\d+(\s+CompServo\s)", r"\g<1>0\g<2>", text, flags=re.M)
        with open(fst, "w") as fh:
            fh.write(text)

    if mode == "nocoords":
        strip_coords(work,
                     os.path.join(baseline, "Airfoils"),
                     os.path.join(work, CW, "NRELOffshrBsline5MW_Onshore_AeroDyn.dat"))


def run(repo, work, case):
    exe = os.path.join(repo, "build-docker-double", "glue-codes", "openfast", "openfast")
    if not os.path.isfile(exe):
        raise SystemExit(f"no openfast binary at {exe}")
    d = os.path.join(work, case)
    with open(os.path.join(d, "run.log"), "w") as log:
        rc = subprocess.run([exe, case + ".fst"], cwd=d, stdout=log, stderr=log).returncode
    if rc != 0:
        raise SystemExit(f"{case} exited {rc}; see {d}/run.log")


def surfaces(work, case, blade):
    """Every frame of one blade's surface, in file order."""
    vtk = os.path.join(work, case, "vtk")
    pat = re.compile(rf"^{re.escape(case)}\.Blade_R1B{blade}Surface\.(\d+)\.vtp$")
    found = []
    for fn in sorted(os.listdir(vtk)):
        m = pat.match(fn)
        if m:
            found.append((int(m.group(1)), os.path.join(vtk, fn)))
    return [p for _, p in sorted(found)]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tmax", type=float, default=2.0, help="run length (s), default 2")
    ap.add_argument("--tol", type=float, default=TOL, help=f"tolerance (m), default {TOL}")
    ap.add_argument("--keep", action="store_true", help="keep the run directories")
    ap.add_argument("--only", choices=MODES, help="check just one configuration")
    args = ap.parse_args(argv)

    repo = repo_root()
    modes = (args.only,) if args.only else MODES
    ok = True

    for mode in modes:
        work = tempfile.mkdtemp(prefix=f"vtk_surface_{mode}_")
        try:
            stage(repo, work, args.tmax, mode)
            for case in (CW, MR):
                run(repo, work, case)

            for i, j in PAIRS.items():
                a_files = surfaces(work, CW, i)
                b_files = surfaces(work, MR, j)
                if not a_files:
                    print(f"FAIL  {mode:9s} blade {i}: no surface files written")
                    ok = False
                    continue
                if len(a_files) != len(b_files):
                    print(f"FAIL  {mode:9s} blade {i}->{j}: "
                          f"{len(a_files)} frames vs {len(b_files)}")
                    ok = False
                    continue

                worst = 0.0
                bad = False
                for pa, pb in zip(a_files, b_files):
                    a = read_vtp(pa)["points"] @ S.T
                    b = read_vtp(pb)["points"]
                    if a.shape != b.shape:
                        print(f"FAIL  {mode:9s} blade {i}->{j}: vertex counts differ, "
                              f"{a.shape[0]} vs {b.shape[0]}")
                        ok = False
                        bad = True
                        break
                    worst = max(worst, float(np.abs(a - b).max()))
                if bad:
                    continue
                good = worst <= args.tol
                ok &= good
                print(f"{'PASS' if good else 'FAIL'}  {mode:9s} blade {i}->{j}  "
                      f"{len(a_files)} frames  {a_files and len(read_vtp(a_files[0])['points'])} verts  "
                      f"max separation {worst:.3e} m  (tol {args.tol:.1e})")
        finally:
            if args.keep:
                print(f"  runs kept in {work}")
            else:
                shutil.rmtree(work, ignore_errors=True)

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
