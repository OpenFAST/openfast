"""
Reader/writer for TurbSim binary full-field (.bts) files, matching the format
implemented in OpenFAST:
  - Writer: modules/turbsim/src/TS_FileIO.f90  (SUBROUTINE WrBinTURBSIM)
  - Reader: modules/inflowwind/src/InflowWind_IO.f90  (SUBROUTINE IfW_TurbSim_Init)

Header layout (stream/unformatted, little-endian, 70 bytes total):
    int16   FileID            (7 = non-periodic, 8 = periodic)
    int32   NZGrids
    int32   NYGrids
    int32   NTGrids           (number of tower points, may be 0)
    int32   NSteps
    float32 dz
    float32 dy
    float32 dt
    float32 MeanWS            (mws, hub-height mean wind speed)
    float32 RefHeight         (hub height)
    float32 GridBase          (height of bottom of grid)
    float32 VslopeX, VoffsetX (U-component scale/offset)
    float32 VslopeY, VoffsetY (V-component scale/offset)
    float32 VslopeZ, VoffsetZ (W-component scale/offset)
    int32   DescLen           (number of ASCII bytes in description string)

Followed by:
    DescLen bytes  -> ASCII description string (no null terminator)

Then, for each of NSteps time steps:
    int16 grid data, shape (NZGrids, NYGrids, 3) in C order
        -> for iz in 0..NZGrids-1: for iy in 0..NYGrids-1: for ic in (U,V,W)
        -> IY varies fastest, then IZ; component (U,V,W) is fastest of all
    if NTGrids > 0:
        int16 tower data, shape (NTGrids, 3) in C order

Grid Y-index convention (confirmed from InflowWind_IO.f90, IfW_FlowField.f90):
    Y_grid_index (1-based) = (Y_position + YHWid) / dy + 1
    => iy = 0 (0-based) corresponds to Y = -YHWid (most negative Y)
    => iy = NYGrids-1    corresponds to Y = +YHWid (most positive Y)
    i.e. increasing iy <=> increasing Y coordinate.

Wind components are ordered (U, V, W) = (longitudinal-X, lateral-Y, vertical-Z).
Tower points all lie on the Y=0 centerline (X,Y=0 column, various Z below grid base).
"""

import struct
from dataclasses import dataclass
from typing import Optional

import numpy as np

HEADER_STRUCT = struct.Struct("<h iiii fff fff ff ff ff i".replace(" ", ""))
HEADER_SIZE = HEADER_STRUCT.size  # 70 bytes

_HEADER_FIELDS = [
    "FileID",
    "NZGrids", "NYGrids", "NTGrids", "NSteps",
    "dz", "dy", "dt",
    "MeanWS", "RefHeight", "GridBase",
    "VslopeX", "VoffsetX",
    "VslopeY", "VoffsetY",
    "VslopeZ", "VoffsetZ",
    "DescLen",
]

assert HEADER_SIZE == 70, f"unexpected header size {HEADER_SIZE}"


@dataclass
class BTSFile:
    FileID: int
    NZGrids: int
    NYGrids: int
    NTGrids: int
    NSteps: int
    dz: float
    dy: float
    dt: float
    MeanWS: float
    RefHeight: float
    GridBase: float
    VslopeX: float
    VoffsetX: float
    VslopeY: float
    VoffsetY: float
    VslopeZ: float
    VoffsetZ: float
    DescStr: str
    # grid velocity, raw int16, shape (NSteps, NZGrids, NYGrids, 3)
    grid: np.ndarray
    # tower velocity, raw int16, shape (NSteps, NTGrids, 3), or None if NTGrids == 0
    tower: Optional[np.ndarray]

    @property
    def Vslope(self):
        return (self.VslopeX, self.VslopeY, self.VslopeZ)

    @property
    def Voffset(self):
        return (self.VoffsetX, self.VoffsetY, self.VoffsetZ)

    def decode_grid(self) -> np.ndarray:
        """Return de-normalized (physical, m/s) grid velocities, shape (NSteps, NZGrids, NYGrids, 3), float32."""
        out = np.empty(self.grid.shape, dtype=np.float32)
        for ic in range(3):
            out[..., ic] = (self.grid[..., ic].astype(np.float32) - self.Voffset[ic]) / self.Vslope[ic]
        return out

    def decode_tower(self) -> Optional[np.ndarray]:
        if self.tower is None:
            return None
        out = np.empty(self.tower.shape, dtype=np.float32)
        for ic in range(3):
            out[..., ic] = (self.tower[..., ic].astype(np.float32) - self.Voffset[ic]) / self.Vslope[ic]
        return out


def read_bts(path: str) -> BTSFile:
    with open(path, "rb") as f:
        raw_header = f.read(HEADER_SIZE)
        if len(raw_header) != HEADER_SIZE:
            raise IOError(f"Could not read {HEADER_SIZE}-byte header from {path}")
        vals = HEADER_STRUCT.unpack(raw_header)
        h = dict(zip(_HEADER_FIELDS, vals))

        desc_bytes = f.read(h["DescLen"])
        if len(desc_bytes) != h["DescLen"]:
            raise IOError(f"Could not read {h['DescLen']}-byte description string from {path}")
        desc_str = desc_bytes.decode("ascii", errors="replace")

        n_grid_vals_per_step = 3 * h["NYGrids"] * h["NZGrids"]
        n_twr_vals_per_step = 3 * h["NTGrids"]

        grid = np.empty((h["NSteps"], h["NZGrids"], h["NYGrids"], 3), dtype=np.int16)
        tower = None
        if h["NTGrids"] > 0:
            tower = np.empty((h["NSteps"], h["NTGrids"], 3), dtype=np.int16)

        for it in range(h["NSteps"]):
            buf = f.read(2 * n_grid_vals_per_step)
            if len(buf) != 2 * n_grid_vals_per_step:
                raise IOError(f"Truncated grid data at time step {it} in {path}")
            grid[it] = np.frombuffer(buf, dtype="<i2").reshape(h["NZGrids"], h["NYGrids"], 3)

            if h["NTGrids"] > 0:
                buf = f.read(2 * n_twr_vals_per_step)
                if len(buf) != 2 * n_twr_vals_per_step:
                    raise IOError(f"Truncated tower data at time step {it} in {path}")
                tower[it] = np.frombuffer(buf, dtype="<i2").reshape(h["NTGrids"], 3)

        trailing = f.read()
        if len(trailing) != 0:
            raise IOError(f"{len(trailing)} unexpected trailing bytes in {path}")

    del h["DescLen"]  # derived from len(DescStr) on write, not stored as a field
    return BTSFile(DescStr=desc_str, grid=grid, tower=tower, **h)


def write_bts(path: str, bts: BTSFile) -> None:
    desc_bytes = bts.DescStr.encode("ascii")
    desc_len = len(desc_bytes)

    header_vals = [
        bts.FileID,
        bts.NZGrids, bts.NYGrids, bts.NTGrids, bts.NSteps,
        bts.dz, bts.dy, bts.dt,
        bts.MeanWS, bts.RefHeight, bts.GridBase,
        bts.VslopeX, bts.VoffsetX,
        bts.VslopeY, bts.VoffsetY,
        bts.VslopeZ, bts.VoffsetZ,
        desc_len,
    ]

    with open(path, "wb") as f:
        f.write(HEADER_STRUCT.pack(*header_vals))
        f.write(desc_bytes)

        for it in range(bts.NSteps):
            f.write(np.ascontiguousarray(bts.grid[it], dtype="<i2").tobytes())
            if bts.NTGrids > 0:
                f.write(np.ascontiguousarray(bts.tower[it], dtype="<i2").tobytes())


def _renormalize(v_phys: np.ndarray, slope: float, offset: float) -> np.ndarray:
    """Encode physical (m/s) values back to clamped, rounded int16 raw values."""
    raw = np.round(v_phys * slope + offset)
    raw = np.clip(raw, -32768, 32767)
    return raw.astype(np.int16)


def mirror_xz(bts: BTSFile, note: Optional[str] = None) -> BTSFile:
    """Return a new BTSFile mirrored across the XZ plane (Y -> -Y).

    - Reverses the Y-index ordering of the grid data (grid point at +Y swaps
      with the corresponding point at -Y).
    - Flips the sign of the V (lateral, Y) velocity component for both the
      grid and tower data (physical velocity, decoded/re-encoded using the
      file's own slope/offset so the int16 raw values remain internally
      consistent).
    - Tower points lie on the Y=0 centerline, so only the V-component sign
      flip applies there (no reordering needed).

    Note on V scaling: the original file's VslopeY/VoffsetY were calibrated
    (by TurbSim) to the min/max of the *original* (generally asymmetric) V
    data. After negating V, those extremes land on the opposite side of the
    range and would clip against the original int16 scaling. So VslopeY and
    VoffsetY are recomputed here from the mirrored V data's own min/max,
    following TurbSim's own scaling formula (WrBinTURBSIM in TS_FileIO.f90).
    U and W values are unchanged by mirroring (only their Y location moves),
    so their slope/offset are left untouched.
    """
    IntMin, IntMax = -32768.0, 32767.0
    IntRng = IntMax - IntMin

    grid_phys = bts.decode_grid()
    grid_phys = grid_phys[:, :, ::-1, :]          # reverse Y index
    grid_phys[..., 1] = -grid_phys[..., 1]        # flip V component

    twr_phys = None
    if bts.tower is not None:
        twr_phys = bts.decode_tower()
        twr_phys[..., 1] = -twr_phys[..., 1]

    # Recompute V-component slope/offset from the mirrored data's min/max
    # (grid + tower combined), matching TurbSim's own scaling algorithm.
    v_vals = [grid_phys[..., 1]]
    if twr_phys is not None:
        v_vals.append(twr_phys[..., 1])
    v_min = min(float(v.min()) for v in v_vals)
    v_max = max(float(v.max()) for v in v_vals)
    if v_max == v_min:
        new_VslopeY = 1.0
    else:
        new_VslopeY = IntRng / (v_max - v_min)
    new_VoffsetY = IntMin - new_VslopeY * v_min

    new_slope = (bts.VslopeX, new_VslopeY, bts.VslopeZ)
    new_offset = (bts.VoffsetX, new_VoffsetY, bts.VoffsetZ)

    new_grid = np.empty_like(bts.grid)
    for ic in range(3):
        new_grid[..., ic] = _renormalize(grid_phys[..., ic], new_slope[ic], new_offset[ic])

    new_tower = None
    if twr_phys is not None:
        new_tower = np.empty_like(bts.tower)
        for ic in range(3):
            new_tower[..., ic] = _renormalize(twr_phys[..., ic], new_slope[ic], new_offset[ic])

    new_desc = bts.DescStr
    if note:
        new_desc = (new_desc + " " + note) if new_desc else note

    return BTSFile(
        FileID=bts.FileID,
        NZGrids=bts.NZGrids, NYGrids=bts.NYGrids, NTGrids=bts.NTGrids, NSteps=bts.NSteps,
        dz=bts.dz, dy=bts.dy, dt=bts.dt,
        MeanWS=bts.MeanWS, RefHeight=bts.RefHeight, GridBase=bts.GridBase,
        VslopeX=bts.VslopeX, VoffsetX=bts.VoffsetX,
        VslopeY=new_VslopeY, VoffsetY=new_VoffsetY,
        VslopeZ=bts.VslopeZ, VoffsetZ=bts.VoffsetZ,
        DescStr=new_desc,
        grid=new_grid,
        tower=new_tower,
    )
