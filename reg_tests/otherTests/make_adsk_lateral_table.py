"""Derive an AeroDisk coefficient table that carries lateral coefficients.

The shipped 5MW table has C_Fy, C_Fz, C_My and C_Mz identically zero, so the
mirror-related sign factors on Force(2), Force(3), Moment(2) and Moment(3) -
and the z_hat pseudovector flip that multiplies them - never execute.  This
scales the lateral columns off the axial ones so every component is non-zero
while keeping the TSR and pitch dependence of the original.

The ratios are arbitrary but fixed; the table is a regression fixture, not a
physical rotor.
"""
import io, sys

RATIO = {6: -0.08, 7: 0.05, 9: 0.06, 10: -0.04}   # col idx -> factor
SRC_COL = {6: 5, 7: 5, 9: 8, 10: 8}               # lateral col -> driving col

def main(src, dst):
    out = []
    for n, line in enumerate(io.open(src, encoding="utf-8")):
        s = line.rstrip("\n")
        if n < 3 or not s.strip() or s.lstrip().startswith("#"):
            out.append(s); continue
        f = [c.strip() for c in s.split(",")]
        if len(f) < 11:
            out.append(s); continue
        for col, fac in RATIO.items():
            f[col] = "%10.4f" % (fac * float(f[SRC_COL[col]]))
        out.append(", ".join("%10s" % c for c in f))
    io.open(dst, "w", encoding="utf-8").write("\n".join(out) + "\n")
    print("wrote", dst, len(out), "lines")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
