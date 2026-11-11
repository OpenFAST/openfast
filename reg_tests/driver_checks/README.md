# driver_checks

Base decks for `reg_tests/otherTests/sweep_mirror.py`: a steady, uniform-inflow
AeroDyn driver case on the NREL 5MW rotor, swept over pitch, wind speed, yaw and
rotor speed as clockwise/mirrored pairs.

- `ad_CW.dvr` — the clockwise driver deck; the sweep writes the mirrored twin by
  editing `MirrorRotor(1)`.
- `ad_primary.dat` — the AeroDyn primary file. Airfoil and blade references are
  written as `@RTEST@/...`, which the sweep replaces with the checkout's
  `reg_tests/r-test/glue-codes/openfast` path when it stages a case, so the decks
  work from any checkout.

Nothing here runs automatically; see `reg_tests/otherTests/README.md`.
