2025.11.10

This is a work in progress. Some things are not complete or functional yet:


Froude scaling is not complete, nor is it tested!!!!

At present, only some inputs are scaled, but equations
have not been verified yet. This has been disabled by
removing the reading of the `*Fact` input lines in the
input file parsing and input file.

TODO:
  - verify equations in FroudeScaling* functions
  - scale resulting forces / moments
  - add scaled time, pos, vel, acc, frc, mom to output
    channels and add subscripting to differentiate
