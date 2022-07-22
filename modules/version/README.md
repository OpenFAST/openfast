# Version Module

## Overview
The Version module provides all driver and glue codes with the version based
on the git status. OpenFAST follows [semantic versioning](https://semver.org).
In summary, this means that with a version number as MAJOR.MINOR.PATCH, the
components will be incremented as follows:

- MAJOR version when introducing incompatible API changes,
- MINOR version when adding functionality in a backwards-compatible manner, and
- PATCH version when making backwards-compatible bug fixes.

For example, ``OpenFAST-v1.0.0-123-gabcd1234-dirty`` describes OpenFAST as:

- v1.0.0 is the MAJOR.MINOR.PATCH numbering system and corresponds to a tagged
  commit made by NREL on GitHub
- 123-g is the number of additional commits after the most recent tag for a
  build [the ``-g`` is for ``git``]
- abcd1234 is the first 8 characters of the current commit hash
- dirty denotes that local changes have been made but not committed

To make use of version information when checking command-line arguments, 
the `CheckArgs` routine that was previously in NWTC_Library has been moved
to this module. This ensures that the version information does not have to
be one of the first modules to be compiled. (The goal is that NWTC Library 
does not depend on the VersionInfo module.)
