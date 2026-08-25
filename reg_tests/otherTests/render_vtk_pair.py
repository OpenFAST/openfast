#!/usr/bin/env python3
"""Render an OpenFAST VTK mesh pair for visual inspection of the mirror.

Every other check in the MirrorRotor work is numerical, and a tolerance cannot
see a blade drawn on the wrong side, a rotor turning the right way with its
geometry built the wrong way round, or a mesh that mirrors correctly in its own
frame but is attached to the hub backwards.  This renders the geometry so that a
human can look at it.

The container this was written in has no OpenGL -- GLX, EGL and OSMesa are all
absent, so VTK cannot open a render window -- and the VTK python module is
therefore not used at all.  The ``.vtp`` files OpenFAST writes are ASCII XML
PolyData, so they are parsed directly and drawn with matplotlib.

Three figures are produced for a pair of case directories:

``geometry.png``
    The clockwise and mirrored geometry side by side at one instant, in three
    views, on identical axes.  Look at these first, and independently: they are
    the only output that can show an error the mirror operator itself shares.

``overlay.png``
    The mirrored geometry with ``y -> -y`` re-applied, drawn over the clockwise
    geometry.  Where the implementation is right the two coincide.  This is much
    the easier judgement to make, but it is a numerical test wearing a visual
    costume, so it does not replace ``geometry.png``.

``rotation.png``
    Blade tip paths over the run, viewed looking downwind, coloured by time.
    This is what shows the sense of rotation.

Usage::

    python3 render_vtk_pair.py --cw <case-dir> --mirror <case-dir> --out <dir>

Each case directory is the one holding the ``vtk/`` subdirectory that OpenFAST
wrote.  ``--frame`` picks the instant for the first two figures.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection


# The mirror operator.  See docs/source/user/glue-code/mirror_rotor.rst.
S = np.diag([1.0, -1.0, 1.0])


# ---------------------------------------------------------------------------
# reading


def _array(el):
    """Text of a DataArray as a flat float array."""
    return np.fromstring(el.text.replace("\n", " "), sep=" ")


def _cells(piece, tag):
    """(connectivity, offsets) for a Lines or Polys block, or None."""
    block = piece.find(tag)
    if block is None:
        return None
    named = {da.get("Name"): da for da in block.findall("DataArray")}
    if "connectivity" not in named or "offsets" not in named:
        return None
    conn = _array(named["connectivity"]).astype(int)
    offs = _array(named["offsets"]).astype(int)
    if conn.size == 0:
        return None
    return conn, offs


def _split(conn, offs):
    """Split a connectivity array at its offsets into per-cell index arrays."""
    return [conn[a:b] for a, b in zip(np.concatenate(([0], offs[:-1])), offs)]


def read_vtp(path):
    """Parse one OpenFAST .vtp file.

    Returns a dict with ``points`` (n,3), ``lines`` and ``polys`` (lists of
    index arrays) and ``fields`` (name -> (n,3) array) for whatever point data
    the file carries.  OpenFAST writes orientation and, with VTK_fields, motion
    and load fields; all are read the same way.
    """
    piece = ET.parse(path).getroot().find(".//Piece")
    if piece is None:
        raise ValueError(f"{path}: no Piece element")

    pts_el = piece.find("Points/DataArray")
    points = _array(pts_el).reshape(-1, 3)

    fields = {}
    pd = piece.find("PointData")
    if pd is not None:
        for da in pd.findall("DataArray"):
            name = da.get("Name")
            ncomp = int(da.get("NumberOfComponents", "1"))
            if name:
                fields[name] = _array(da).reshape(-1, ncomp)

    lines = _cells(piece, "Lines")
    polys = _cells(piece, "Polys")
    return {
        "points": points,
        "lines": _split(*lines) if lines else [],
        "polys": _split(*polys) if polys else [],
        "fields": fields,
    }


# VTK_tWidth scales with the run length, so the frame field is not always four
# digits; a greedy root takes the last dotted number group, which is the frame.
_FRAME_RE = re.compile(r"^(?P<root>.+)\.(?P<frame>\d+)\.vtp$")


def scan(case_dir):
    """Index a case's vtk directory.

    Returns (meshes, frames) where meshes maps a mesh name -- the file root with
    the case prefix stripped -- to {frame_number: path}, and frames is the
    sorted list of frames present in every mesh.  Reference and static files are
    indexed under frame None.
    """
    vtk_dir = os.path.join(case_dir, "vtk")
    if not os.path.isdir(vtk_dir):
        raise SystemExit(f"no vtk/ directory in {case_dir}")

    meshes = defaultdict(dict)
    for fn in sorted(os.listdir(vtk_dir)):
        if not fn.endswith(".vtp"):
            continue
        path = os.path.join(vtk_dir, fn)
        m = _FRAME_RE.match(fn)
        if m:
            meshes[m.group("root")][int(m.group("frame"))] = path
        else:
            meshes[fn[: -len(".vtp")]][None] = path

    # Strip the leading case name, which differs between the two members of a
    # pair and would otherwise stop the meshes from lining up.
    prefix = os.path.basename(os.path.normpath(case_dir)) + "."
    stripped = {}
    for name, byframe in meshes.items():
        stripped[name[len(prefix):] if name.startswith(prefix) else name] = byframe

    animated = [f for m in stripped.values() for f in m if f is not None]
    frames = sorted(set(animated))
    return stripped, frames


# ---------------------------------------------------------------------------
# what each mesh is, and how to draw it


def vtk_dt(case_dir):
    """Seconds between VTK frames, from the case summary file.

    OpenFAST rounds 1/VTK_fps to an integer multiple of DT, so the frame
    spacing is generally not 1/VTK_fps and cannot be taken from the deck.  The
    summary reports the effective rate directly.  Note that the time-step table
    just above it truncates the exponent off every row that carries a subcycle
    count -- see FAST_Subs.f90, the T37 in the format string -- so the rate line
    is the one to read, not the table.
    """
    for fn in sorted(os.listdir(case_dir)):
        if not fn.endswith(".sum"):
            continue
        with open(os.path.join(case_dir, fn), errors="ignore") as fh:
            for line in fh:
                if "Frame rate" in line:
                    m = re.search(r"([0-9]+\.[0-9]+)\s*fps", line)
                    if m and float(m.group(1)) > 0:
                        return 1.0 / float(m.group(1))
    return None


def when(dt, frame):
    """Label for a frame: its time if we know the rate, else the index."""
    if dt is None:
        return f"frame {frame}"
    return f"t = {dt * frame:.2f} s (frame {frame})"


def style(name):
    """Colour, line width and z-order for a mesh, by family."""
    if "GroundSurface" in name:
        return "0.80", 0.8, 0
    # VTK_type = 1 writes surfaces rather than the debug line meshes.
    if "BladeSurface" in name or ("Blade" in name and name.endswith("Surface")):
        return "#1f77b4", 0.4, 3
    if "TowerSurface" in name:
        return "0.35", 0.4, 1
    if "NacelleSurface" in name or "HubSurface" in name:
        return "#ff7f0e", 0.4, 5
    if name.startswith("BD_Blade"):
        return "#1f77b4", 2.2, 3          # the structural blade
    if name.startswith("AD_Blade_"):
        return "#d62728", 1.4, 2          # the aerodynamic blade
    if "BladeRootMotion" in name:
        return "#2ca02c", 3.0, 4
    if "ReactionForce" in name:
        return "#9467bd", 3.0, 4
    if "HubMotion" in name:
        return "#ff7f0e", 3.0, 5
    if "Tower" in name:
        return "0.35", 2.0, 1
    return "0.5", 1.0, 1


def collect(meshes, frame, want=None):
    """Read every mesh at a frame.  Falls back to the static file if a mesh has
    no animated frames, so the ground plane and any reference-only mesh still
    appear."""
    out = {}
    for name, byframe in sorted(meshes.items()):
        if want is not None and not want(name):
            continue
        path = byframe.get(frame)
        if path is None:
            path = byframe.get(None)
        if path is None:
            continue
        out[name] = read_vtp(path)
    return out


def draw(ax, data, transform=None, colour=None, alpha=1.0, points=True):
    """Draw a set of meshes into a 3-D axis.

    ``transform`` is applied to every point, which is how the overlay re-applies
    the mirror.  ``colour`` overrides the per-family colour, which is how the
    overlay tells the two turbines apart.
    """
    for name, mesh in data.items():
        p = mesh["points"]
        if transform is not None:
            p = p @ transform.T
        base, lw, zo = style(name)
        c = colour or base

        segs = [p[idx] for idx in mesh["lines"] if len(idx) >= 2]
        if segs:
            ax.add_collection3d(
                Line3DCollection(segs, colors=c, linewidths=lw, alpha=alpha, zorder=zo)
            )
        faces = [p[idx] for idx in mesh["polys"] if len(idx) >= 3]
        if faces:
            ax.add_collection3d(
                Poly3DCollection(
                    faces, facecolors=c, edgecolors="none", alpha=alpha * 0.35, zorder=zo
                )
            )
        # A mesh with no cells is a single node -- the hub, a blade root -- and
        # is invisible unless its points are drawn.
        if points and not segs and not faces:
            ax.scatter(p[:, 0], p[:, 1], p[:, 2], c=c, s=28, alpha=alpha, zorder=zo,
                       depthshade=False)


def draw2d(ax, data, axes, transform=None, colour=None, alpha=1.0):
    """Draw a set of meshes as a flat orthographic projection.

    ``axes`` picks the two coordinate indices to plot.  True coordinates are
    plotted and the axis is reversed by ``setup2d`` where a view needs it, so
    that the tick labels keep meaning what they say; negating the data instead
    puts a node at y = +63 under a tick reading -60.

    These views are planar, so drawing them in a 3-D axis only buys
    matplotlib's pane and tick furniture over the top of the geometry.
    """
    i, j = axes
    si = sj = 1
    for name, mesh in data.items():
        p = mesh["points"]
        if transform is not None:
            p = p @ transform.T
        base, lw, zo = style(name)
        c = colour or base
        u, v = si * p[:, i], sj * p[:, j]

        segs = [np.column_stack((u[idx], v[idx])) for idx in mesh["lines"] if len(idx) >= 2]
        if segs:
            ax.add_collection(
                LineCollection(segs, colors=c, linewidths=lw, alpha=alpha, zorder=zo)
            )
        faces = [np.column_stack((u[idx], v[idx])) for idx in mesh["polys"] if len(idx) >= 3]
        if faces:
            ax.add_collection(
                PolyCollection(faces, facecolors=c, edgecolors="none",
                               alpha=alpha * 0.25, zorder=zo)
            )
        if not segs and not faces:
            ax.scatter(u, v, c=c, s=30, alpha=alpha, zorder=zo)


# Two planar views plus one isometric.  (label, axis indices, which axes to
# reverse, axis labels).  Looking downwind means +x into the page; the frame is
# right-handed with z up, so +y runs to the left and the horizontal axis of
# that view is reversed.
FLAT_VIEWS = [
    ("looking downwind, +x into page", (1, 2), (True, False),
     "y (m), increasing to the left", "z (m)"),
    ("plan, from above", (0, 1), (False, False), "x (m), downwind", "y (m)"),
]


def frame_limits(*datasets):
    """A single cubic bounding box covering everything, so that the two members
    of a pair are drawn to exactly the same scale and cannot be compared
    misleadingly."""
    pts = []
    for data in datasets:
        for mesh in data.values():
            pts.append(mesh["points"])
    allp = np.vstack(pts)
    # The ground plane is huge and would shrink the turbine to nothing.
    turbine = allp[allp[:, 2] > 1.0]
    if len(turbine) > 8:
        allp = turbine
    lo, hi = allp.min(axis=0), allp.max(axis=0)
    ctr = 0.5 * (lo + hi)
    half = 0.55 * (hi - lo).max()
    return np.array([ctr - half, ctr + half])


def setup(ax, lim, title, elev=22, azim=-125):
    ax.set_xlim(lim[0, 0], lim[1, 0])
    ax.set_ylim(lim[0, 1], lim[1, 1])
    ax.set_zlim(lim[0, 2], lim[1, 2])
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel("x", fontsize=7)
    ax.set_ylabel("y", fontsize=7)
    ax.set_zlabel("z", fontsize=7)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=5)


def setup2d(ax, lim, view, title):
    _, (i, j), (rx, ry), xl, yl = view
    ax.set_xlim(lim[0, i], lim[1, i])
    ax.set_ylim(lim[0, j], lim[1, j])
    if rx:
        ax.invert_xaxis()
    if ry:
        ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xlabel(xl, fontsize=8)
    ax.set_ylabel(yl, fontsize=8)
    ax.set_title(title, fontsize=9)
    ax.tick_params(labelsize=7)
    ax.grid(alpha=0.25, lw=0.5)


# ---------------------------------------------------------------------------
# figures


LEGEND = ("blue BD blade, red AD blade, green blade root, "
          "purple BD reaction, orange hub, grey tower")


def fig_geometry(cw, mr, lim, label, out):
    fig = plt.figure(figsize=(15, 9.5))
    for row, (case, data) in enumerate((("clockwise", cw), ("mirrored", mr))):
        for col, view in enumerate(FLAT_VIEWS):
            ax = fig.add_subplot(2, 3, row * 3 + col + 1)
            draw2d(ax, data, view[1])
            setup2d(ax, lim, view, f"{case} -- {view[0]}")
        ax = fig.add_subplot(2, 3, row * 3 + 3, projection="3d")
        draw(ax, data)
        setup(ax, lim, f"{case} -- isometric")
    fig.suptitle(
        f"Geometry at {label}.  Identical axes throughout.  {LEGEND}",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out, dpi=130)
    plt.close(fig)


def fig_overlay(cw, mr, lim, label, out, deviations):
    fig = plt.figure(figsize=(15, 5.5))
    for col, view in enumerate(FLAT_VIEWS):
        ax = fig.add_subplot(1, 3, col + 1)
        draw2d(ax, cw, view[1], colour="#1f77b4", alpha=1.0)
        draw2d(ax, mr, view[1], transform=S, colour="#d62728", alpha=0.6)
        setup2d(ax, lim, view, view[0])
    ax = fig.add_subplot(1, 3, 3, projection="3d")
    draw(ax, cw, colour="#1f77b4", alpha=1.0)
    draw(ax, mr, transform=S, colour="#d62728", alpha=0.6)
    setup(ax, lim, "isometric")
    worst = max((v for v, _ in deviations.values()), default=float("nan"))
    fig.suptitle(
        f"Mirrored geometry with y -> -y re-applied (red) over clockwise (blue), {label}.  "
        f"They should coincide -- no blue should be visible.  "
        f"Largest node separation {worst:.3e} m",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=130)
    plt.close(fig)


def tip_paths(meshes, frames, pattern="AD_Blade_"):
    """Tip position of each blade over the run."""
    out = {}
    for name, byframe in sorted(meshes.items()):
        if not name.startswith(pattern):
            continue
        pts = []
        for f in frames:
            path = byframe.get(f)
            if path is None:
                continue
            # The last point of the blade line mesh is the tip.
            pts.append(read_vtp(path)["points"][-1])
        if pts:
            out[name] = np.array(pts)
    return out


def azimuth(path, hub):
    """Unwrapped blade azimuth about the shaft, radians.

    Measured in the plane normal to the wind, from vertical, positive in the
    sense a clockwise rotor turns when viewed from upwind.  Looking downwind
    with z up the frame is right-handed and +y runs to the left, so that sense
    carries -y into +z, giving atan2(-y, z).
    """
    d = path - hub
    return np.unwrap(np.arctan2(-d[:, 1], d[:, 2]))


def fig_rotation(cw_paths, mr_paths, cw_hub, mr_hub, out):
    """Sense of rotation.

    The earlier version of this figure drew every frame of the run.  The rotor
    turns about four times in twenty seconds, so the tip paths wrapped over
    themselves and the colour gradient aliased into noise that showed nothing at
    all.  It is restricted to a single revolution here, and the direction of
    travel is drawn as arrows taken from consecutive positions rather than left
    to be inferred from colour.  The right-hand panel is the unambiguous one:
    the sign of the slope is the sense of rotation, and no question of how the
    view is oriented enters into it.
    """
    fig = plt.figure(figsize=(15, 6.0))
    cases = (("clockwise", cw_paths, cw_hub, "#1f77b4"),
             ("mirrored", mr_paths, mr_hub, "#d62728"))

    for col, (label, paths, hub, _) in enumerate(cases):
        ax = fig.add_subplot(1, 3, col + 1)
        # One blade only.  Drawing all three puts three colour ramps on one
        # circle, which reads as a ramp that cycles and hides the very thing
        # the panel exists to show.
        for name, p in sorted(paths.items())[:1]:
            psi = azimuth(p, hub)
            # One revolution's worth of frames, or the whole run if shorter.
            adv = np.abs(psi - psi[0])
            n = int(np.argmax(adv >= 2 * np.pi)) + 1 if (adv >= 2 * np.pi).any() else len(p)
            q = p[:n]
            u, v = q[:, 1], q[:, 2]
            ax.scatter(u, v, c=np.arange(n), cmap="viridis", s=14, zorder=3)
            step = max(1, n // 16)
            for k in range(0, n - step, step):
                ax.annotate(
                    "", xy=(u[k + step], v[k + step]), xytext=(u[k], v[k]),
                    arrowprops=dict(arrowstyle="-|>", lw=1.6, color="k",
                                    shrinkA=0, shrinkB=0),
                    zorder=4,
                )
            ax.annotate(name.replace("AD_Blade_R1", ""), (u[0], v[0]),
                        fontsize=9, weight="bold", zorder=5)
        ax.plot([hub[1]], [hub[2]], "k+", ms=14)
        ax.invert_xaxis()   # looking downwind: +y to the left
        ax.set_aspect("equal")
        ax.set_xlabel("y (m), increasing to the left", fontsize=8)
        ax.set_ylabel("z (m)", fontsize=8)
        ax.set_title(f"{label}: blade 1 tip path, one revolution", fontsize=10)
        ax.grid(alpha=0.25, lw=0.5)

    ax = fig.add_subplot(1, 3, 3)
    for label, paths, hub, colour in cases:
        for i, (name, p) in enumerate(sorted(paths.items())):
            psi = np.degrees(azimuth(p, hub) - azimuth(p, hub)[0])
            ax.plot(np.arange(len(psi)), psi, color=colour, lw=1.4,
                    label=label if i == 0 else None)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("frame", fontsize=8)
    ax.set_ylabel("azimuth advance from start (deg)", fontsize=8)
    ax.set_title("unwrapped azimuth: the slope is the sense", fontsize=10)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(fontsize=8)

    fig.suptitle(
        "Sense of rotation, viewed from upwind looking downwind (+x into page).  "
        "The two rotors must turn in opposite senses.",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=130)
    plt.close(fig)


def hub_position(meshes, frame):
    """Hub node position at a frame, for the azimuth origin."""
    for name, byframe in meshes.items():
        if "HubMotion" in name and "_Reference" not in name:
            path = byframe.get(frame) or byframe.get(None)
            if path:
                return read_vtp(path)["points"][0]
    return np.array([0.0, 0.0, 90.0])


# ---------------------------------------------------------------------------



# ---------------------------------------------------------------------------
# animation


class MeshCache:
    """Parsed meshes, kept so each .vtp is read once.

    A run of this length holds four thousand files per case and every figure
    wants the same ones, so without this the encode spends most of its time in
    the XML parser.  The meshes are small -- a blade is nineteen points -- so
    holding all of them costs little.
    """

    def __init__(self, meshes, want=None):
        self.meshes = meshes
        self.want = want
        self._cache = {}

    def at(self, frame):
        if frame not in self._cache:
            self._cache[frame] = collect(self.meshes, frame, self.want)
        return self._cache[frame]


def animation_limits(cw, mr, frames):
    """Axis limits fixed for the whole animation.

    Taken over a spread of frames rather than one, so that nothing wanders out
    of frame part way through -- an axis that rescales mid-clip would make two
    cases look different when they are not.
    """
    sample = frames[:: max(1, len(frames) // 12)] + [frames[-1]]
    return frame_limits(*[cw.at(f) for f in sample], *[mr.at(f) for f in sample])


def rotation_bounds(cw_paths, mr_paths, cw_hub, mr_hub):
    """Fixed axis limits for the rotation animation.

    Every axis in an animation has to be pinned to its final extent.  Left to
    autoscale, the first frames hold one point, so the tip circle appears to
    grow out of nothing and the azimuth trace always touches the top of its
    box -- both of which show motion that is not there.
    """
    tips = np.vstack([p for p in list(cw_paths.values()) + list(mr_paths.values())])
    y, z = tips[:, 1], tips[:, 2]
    pad = 0.08 * max(np.ptp(y), np.ptp(z))
    yz = (y.min() - pad, y.max() + pad, z.min() - pad, z.max() + pad)

    psis = []
    for paths, hub in ((cw_paths, cw_hub), (mr_paths, mr_hub)):
        for p in sorted(paths.values(), key=lambda a: a[0, 1])[:1]:
            a = azimuth(p, hub)
            psis.append(np.degrees(a - a[0]))
    # Unwrapping assumes successive samples are less than half a revolution
    # apart.  A large --stride breaks that: at 0.85 rev per sample the unwrap
    # resolves the wrong way and the azimuth trace understates the turning by
    # a factor of several, quietly.
    for q in psis:
        step = np.abs(np.diff(q)).max() if len(q) > 1 else 0.0
        if step > 150.0:
            print(f"warning: {step:.0f} deg between sampled frames -- the stride "
                  "aliases the azimuth unwrap; use a smaller --stride for a "
                  "trustworthy rotation panel", file=sys.stderr)
            break

    lo = min(float(q.min()) for q in psis)
    hi = max(float(q.max()) for q in psis)
    m = 0.05 * (hi - lo)
    return yz, (lo - m, hi + m)


def anim_rotation(fig, axes, cw_paths, mr_paths, cw_hub, mr_hub, upto, nframes,
                  bounds):
    """One frame of the rotation animation: the trace accumulated so far."""
    (y0, y1, z0, z1), (psi0, psi1) = bounds
    ax_cw, ax_mr, ax_psi = axes
    cases = (("clockwise", cw_paths, cw_hub, ax_cw, "#1f77b4"),
             ("mirrored", mr_paths, mr_hub, ax_mr, "#d62728"))

    for label, paths, hub, ax, _ in cases:
        ax.cla()
        for name, p in sorted(paths.items())[:1]:
            q = p[: upto + 1]
            ax.plot(q[:, 1], q[:, 2], "-", color="0.7", lw=1.0, zorder=2)
            ax.scatter(q[:, 1], q[:, 2], c=np.arange(len(q)), cmap="viridis",
                       s=10, vmin=0, vmax=nframes - 1, zorder=3)
            ax.plot(q[-1, 1], q[-1, 2], "o", color="k", ms=9, zorder=5)
        ax.plot([hub[1]], [hub[2]], "k+", ms=14)
        ax.set_xlim(y0, y1)
        ax.set_ylim(z0, z1)
        ax.set_aspect("equal")
        ax.invert_xaxis()
        ax.set_xlabel("y (m), increasing to the left", fontsize=8)
        ax.set_ylabel("z (m)", fontsize=8)
        ax.set_title(f"{label}: blade 1 tip", fontsize=10)
        ax.grid(alpha=0.25, lw=0.5)

    ax_psi.cla()
    for label, paths, hub, _, colour in cases:
        for i, (name, p) in enumerate(sorted(paths.items())[:1]):
            psi = np.degrees(azimuth(p, hub) - azimuth(p, hub)[0])
            ax_psi.plot(np.arange(upto + 1), psi[: upto + 1], color=colour, lw=1.8,
                        label=label if i == 0 else None)
    ax_psi.axhline(0, color="k", lw=0.6)
    ax_psi.set_xlim(0, nframes - 1)
    ax_psi.set_ylim(psi0, psi1)
    ax_psi.set_xlabel("frame", fontsize=8)
    ax_psi.set_ylabel("azimuth advance from start (deg)", fontsize=8)
    ax_psi.set_title("unwrapped azimuth: the slope is the sense", fontsize=10)
    ax_psi.grid(alpha=0.25, lw=0.5)
    ax_psi.legend(fontsize=8, loc="upper left")


def encode(pattern, out, fps):
    """PNG sequence to H.264 mp4, the format a GitHub PR comment accepts."""
    try:
        import imageio_ffmpeg
        exe = imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        exe = shutil.which("ffmpeg")
    if not exe:
        raise SystemExit("no ffmpeg available; install imageio-ffmpeg or ffmpeg")

    cmd = [
        exe, "-y", "-loglevel", "error",
        "-framerate", f"{fps:.6f}", "-i", pattern,
        # H.264 with yuv420p needs even dimensions, and a figure saved at an
        # arbitrary dpi rarely has them.
        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2:0:0:white",
        "-c:v", "libx264", "-preset", "slow", "-crf", "23",
        "-pix_fmt", "yuv420p", "-movflags", "+faststart",
        out,
    ]
    subprocess.run(cmd, check=True)
    return out


def animate(cw_meshes, mr_meshes, frames, out_dir, want, fps, dpi, dt):
    """Write geometry.mp4, overlay.mp4 and rotation.mp4."""
    cw = MeshCache(cw_meshes, want)
    mr = MeshCache(mr_meshes, want)
    lim = animation_limits(cw, mr, frames)

    tmp = os.path.join(out_dir, "_frames")
    if os.path.isdir(tmp):
        shutil.rmtree(tmp)
    os.makedirs(tmp)

    # Figures are built once and their axes cleared each frame; rebuilding a
    # figure per frame is most of the cost otherwise.
    fig_g = plt.figure(figsize=(15, 9.5))
    ax_g = [fig_g.add_subplot(2, 3, k + 1) for k in (0, 1)] + \
           [fig_g.add_subplot(2, 3, 3, projection="3d")] + \
           [fig_g.add_subplot(2, 3, k + 1) for k in (3, 4)] + \
           [fig_g.add_subplot(2, 3, 6, projection="3d")]
    sup_g = fig_g.suptitle("", fontsize=10)

    fig_o = plt.figure(figsize=(15, 5.5))
    ax_o = [fig_o.add_subplot(1, 3, 1), fig_o.add_subplot(1, 3, 2),
            fig_o.add_subplot(1, 3, 3, projection="3d")]
    sup_o = fig_o.suptitle("", fontsize=10)

    fig_r = plt.figure(figsize=(15, 6.0))
    ax_r = [fig_r.add_subplot(1, 3, k + 1) for k in range(3)]
    sup_r = fig_r.suptitle(
        "Sense of rotation, viewed from upwind looking downwind (+x into page).  "
        "The two rotors must turn in opposite senses.", fontsize=11)

    cw_paths = tip_paths(cw_meshes, frames)
    mr_paths = tip_paths(mr_meshes, frames)
    cw_hub = hub_position(cw_meshes, frames[0])
    mr_hub = hub_position(mr_meshes, frames[0])
    rbounds = rotation_bounds(cw_paths, mr_paths, cw_hub, mr_hub)

    for k, f in enumerate(frames):
        dcw, dmr = cw.at(f), mr.at(f)
        label = when(dt, f)

        for a in ax_g:
            a.cla()
        for col, view in enumerate(FLAT_VIEWS):
            draw2d(ax_g[col], dcw, view[1])
            setup2d(ax_g[col], lim, view, f"clockwise -- {view[0]}")
            draw2d(ax_g[3 + col], dmr, view[1])
            setup2d(ax_g[3 + col], lim, view, f"mirrored -- {view[0]}")
        draw(ax_g[2], dcw)
        setup(ax_g[2], lim, "clockwise -- isometric")
        draw(ax_g[5], dmr)
        setup(ax_g[5], lim, "mirrored -- isometric")
        sup_g.set_text(f"Geometry at {label}.  Identical axes throughout.  {LEGEND}")
        fig_g.savefig(os.path.join(tmp, f"geometry_{k:05d}.png"), dpi=dpi)

        for a in ax_o:
            a.cla()
        for col, view in enumerate(FLAT_VIEWS):
            draw2d(ax_o[col], dcw, view[1], colour="#1f77b4", alpha=1.0)
            draw2d(ax_o[col], dmr, view[1], transform=S, colour="#d62728", alpha=0.6)
            setup2d(ax_o[col], lim, view, view[0])
        draw(ax_o[2], dcw, colour="#1f77b4", alpha=1.0)
        draw(ax_o[2], dmr, transform=S, colour="#d62728", alpha=0.6)
        setup(ax_o[2], lim, "isometric")
        worst = max((v for v, _ in deviation(dcw, dmr).values()),
                    default=float("nan"))
        sup_o.set_text(
            f"Mirrored geometry with y -> -y re-applied (red) over clockwise (blue), "
            f"{label}.  No blue should be visible.  "
            f"Largest node separation {worst:.3e} m")
        fig_o.savefig(os.path.join(tmp, f"overlay_{k:05d}.png"), dpi=dpi)

        anim_rotation(fig_r, ax_r, cw_paths, mr_paths, cw_hub, mr_hub, k,
                      len(frames), rbounds)
        fig_r.savefig(os.path.join(tmp, f"rotation_{k:05d}.png"), dpi=dpi)

        if k % 20 == 0 or k == len(frames) - 1:
            print(f"  frame {k + 1}/{len(frames)}  ({label})", flush=True)

    for fig in (fig_g, fig_o, fig_r):
        plt.close(fig)

    outs = []
    for stem in ("geometry", "overlay", "rotation"):
        out = os.path.join(out_dir, f"{stem}.mp4")
        encode(os.path.join(tmp, f"{stem}_%05d.png"), out, fps)
        outs.append(out)
        print(f"  {out}  {os.path.getsize(out) / 1e6:.2f} MB")

    shutil.rmtree(tmp)
    return outs

# Meshes built by sweeping an angle -- the ground quad, the tower, nacelle and
# hub surfaces -- have their vertex order reversed by the reflection, so a
# node-by-node comparison reports a large separation for two identical shapes.
# They are compared as point clouds instead.
_SWEPT = ("GroundSurface", "TowerSurface", "NacelleSurface", "HubSurface")

# The blade index is not reliably the end of a mesh name: it is followed by
# "_Reference" on reference meshes and by "Surface" on surface meshes.  Matching
# the end of the name silently compared a blade against itself, which shows up as
# a clean 120 degree azimuth error and looks convincingly like a real defect.
_BLADE_RE = re.compile(r"B(\d+)")


def partner(name):
    """The mirrored mesh a clockwise mesh should be compared against.

    Blade 1 lies on the mirror plane and blades 2 and 3 exchange.
    """
    swap = {"2": "3", "3": "2"}
    return _BLADE_RE.sub(lambda m: "B" + swap.get(m.group(1), m.group(1)), name, count=1)


def _cloud_sep(a, b):
    """Largest distance from a point of b to the nearest point of a."""
    step = max(1, len(a) // 1500)
    a, b = a[::step], b[::step]
    d = np.linalg.norm(b[:, None, :] - a[None, :, :], axis=2)
    return float(d.min(axis=1).max())


def deviation(cw, mr):
    """Largest node separation between each clockwise mesh and its mirrored
    partner with the mirror re-applied.

    Returns {name: (separation, ordered)} where ordered is False for the swept
    shapes measured as point clouds.
    """
    out = {}
    for name, mesh in cw.items():
        p = partner(name)
        if p not in mr:
            continue
        a = mesh["points"]
        b = mr[p]["points"] @ S.T
        if a.shape != b.shape:
            continue
        if any(k in name for k in _SWEPT):
            out[name] = (_cloud_sep(a, b), False)
        else:
            out[name] = (float(np.abs(a - b).max()), True)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cw", required=True, help="clockwise case directory")
    ap.add_argument("--mirror", required=True, help="mirrored case directory")
    ap.add_argument("--out", required=True, help="directory to write PNGs into")
    ap.add_argument("--frame", type=int, default=None,
                    help="frame for the geometry and overlay figures "
                         "(default: the last one)")
    ap.add_argument("--no-ground", action="store_true",
                    help="omit the ground plane")
    ap.add_argument("--with-reference", action="store_true",
                    help="also draw the _Reference meshes.  They sit at the "
                         "reference azimuth rather than the frame's, so by "
                         "default they are measured but not drawn")
    ap.add_argument("--animate", action="store_true",
                    help="also write geometry.mp4, overlay.mp4 and rotation.mp4 "
                         "over every frame")
    ap.add_argument("--fps", type=float, default=None,
                    help="playback frame rate (default: the VTK output rate, "
                         "so the clip runs in real time)")
    ap.add_argument("--dpi", type=int, default=100,
                    help="dpi for the animation frames (default 100)")
    ap.add_argument("--stride", type=int, default=1,
                    help="use every Nth frame in the animation (default 1)")
    args = ap.parse_args(argv)

    os.makedirs(args.out, exist_ok=True)

    cw_meshes, cw_frames = scan(args.cw)
    mr_meshes, mr_frames = scan(args.mirror)

    common = sorted(set(cw_frames) & set(mr_frames))
    if not common:
        raise SystemExit("the two cases share no animated frames")
    frame = args.frame if args.frame is not None else common[-1]
    if frame not in common:
        raise SystemExit(f"frame {frame} not present in both cases")

    want = (lambda n: "GroundSurface" not in n) if args.no_ground else None
    cw = collect(cw_meshes, frame, want)
    mr = collect(mr_meshes, frame, want)

    # Every mesh is measured, including the reference ones, but drawing a
    # reference mesh beside the current geometry puts a second copy of the rotor
    # at a different azimuth in the same picture and reads as a defect.
    devs = deviation(cw, mr)
    if not args.with_reference:
        cw = {k: v for k, v in cw.items() if "_Reference" not in k}
        mr = {k: v for k, v in mr.items() if "_Reference" not in k}

    missing = set(cw) ^ set(mr)
    if missing:
        print(f"warning: meshes present in only one case: {sorted(missing)}",
              file=sys.stderr)

    lim = frame_limits(cw, mr)

    dt = vtk_dt(args.cw)
    label = when(dt, frame)

    fig_geometry(cw, mr, lim, label, os.path.join(args.out, "geometry.png"))
    fig_overlay(cw, mr, lim, label, os.path.join(args.out, "overlay.png"), devs)
    fig_rotation(
        tip_paths(cw_meshes, common),
        tip_paths(mr_meshes, common),
        hub_position(cw_meshes, common[0]),
        hub_position(mr_meshes, common[0]),
        os.path.join(args.out, "rotation.png"),
    )

    if dt is not None:
        print(f"frames {common[0]}..{common[-1]} at {dt:.4f} s spacing, "
              f"t = 0 to {dt * common[-1]:.2f} s")
    print(f"rendered {label}")
    print(f"{len(cw)} meshes\n")
    print("largest node separation after re-applying the mirror, per mesh")
    print("(blades 2 and 3 compared against their permuted partners)")
    for name in sorted(devs, key=lambda n: devs[n][0], reverse=True):
        val, ordered = devs[name]
        note = "" if ordered else "   (point cloud; vertex order reverses)"
        print(f"  {name:<44s} {val:.6e} m{note}")
    print(f"\nPNGs in {args.out}")

    if args.animate:
        fps = args.fps if args.fps else (1.0 / dt if dt else 15.0)
        seq = common[:: max(1, args.stride)]
        # A stride changes how much simulated time each played frame covers, so
        # the rate has to follow it or the clip no longer runs in real time.
        if args.stride > 1 and args.fps is None:
            fps = fps / args.stride
        print(f"\nanimating {len(seq)} frames at {fps:.4f} fps "
              f"({len(seq) / fps:.1f} s of video)")
        animate(cw_meshes, mr_meshes, seq, args.out, want, fps, args.dpi, dt)

    return 0


if __name__ == "__main__":
    sys.exit(main())
