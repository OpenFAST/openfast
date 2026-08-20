# otherTests

Verification tooling that is not part of the ctest regression suite but is
needed to justify what the regression suite asserts. These scripts measure
behaviour that a stored baseline cannot express on its own — most of them exist
to support the counter-clockwise rotor work described in
`docs/source/user/glue-code/mirror_rotor.rst`.

Nothing here runs automatically. The scripts are run by hand when a claim needs
re-establishing, or when a change touches the rotor convention.

Each script locates the repository from its own position on disk, so they work
from any checkout. Set `OPENFAST_REPO` to override that.

## The mirror comparison

A mirrored rotor is verified by running a model twice, once clockwise and once
with `MirrorRotor = T`, and requiring **every** output channel to resolve to one
of: identical, exactly sign-flipped, a mirrored angle, or below the numerical
noise floor. Anything else means two quantities have been combined while
expressed in different frames.

| Script | What it does |
|---|---|
| `compare_mirror.py` | the shared classifier. `classify(x, y, tol, floor)` returns `S`, `F`, `A`, `negligible` or `?`. Also holds `blade_permutation`, since blade 1 lies on the mirror plane and blades 2 and 3 exchange. Everything else imports this |
| `emit_sign_table.py` | regenerates the measured sign table in `mirror_rotor.rst` by running every registered mirrored pair and classifying every channel |
| `check_rtest_mirror_pair.py` | checks the registered clockwise/mirrored pairs, with per-pair blade count, tolerance and start time |
| `openfast_case_mirror.py` | clones a regression case verbatim and toggles only `MirrorRotor` |
| `mr_two_rotor.py` | compares rotor 1 against rotor 2 inside one multi-rotor glue-code run |
| `v2_two_rotor.py` | the same for the AeroDyn driver, handling the `A?B<k>N...` nodal naming |
| `sweep_mirror.py`, `openfast_mirror.py`, `openfast_mirror_matrix.py` | sweep the comparison over pitch, wind speed, tip-speed ratio, tilt, precone and the wake and unsteady-aerodynamics models |
| `sed_canary.py` | the SimplifiedElastoDyn brake check |
| `run_guards.sh` | proves the `MirrorRotor` guard rails fire, and fire only for a mirrored rotor |
| `make_adsk_lateral_table.py` | derives an AeroDisk coefficient table carrying lateral coefficients, since the shipped 5MW table has them identically zero |
| `bts_io.py`, `mirror_bts.py`, `check_mirror_bts.py` | reflect a TurbSim box in `y` and verify the result. `bts_io.py` also documents the `.bts` format |

## Reading output files

`compare_mirror.main()` reads **text** `.out` files. For `.outb`, import
`classify` and use `reg_tests/lib/pass_fail.readFASTOut` directly, as the other
scripts do.

`.outb` written with `OutFileFmt = 2` packs each channel into int16, giving
about `1/65536`, or `1.5e-5`, relative resolution. **Do not claim agreement
tighter than that from such a file.** For precision work set `OutFileFmt = 4`.
Text output is worse, not better, at the default `OutFmt`.

## Other checks

| Script | What it does |
|---|---|
| `check_ed_outparams.py`, `check_bd_vs_ed.py` | cross-check output parameter tables against the source |
| `check_outlist_dupes.py`, `dump_outlist_xlsx.py` | check `OutListParameters.xlsx` against the registered channels |
| `make_rtest_cases.py`, `rtest_regression.py`, `openfast_regression.py`, `fix_dvr_comment.py` | case generation and bulk regression helpers |
