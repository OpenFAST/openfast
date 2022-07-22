.. _versioning:

Versioning
~~~~~~~~~~
OpenFAST follows `semantic versioning <https://semver.org>`_. In summary, this
means that with a version number as MAJOR.MINOR.PATCH, the components will be
incremented as follows:

- MAJOR version when introducing incompatible API changes,
- MINOR version when adding functionality in a backwards-compatible manner, and
- PATCH version when making backwards-compatible bug fixes.

For example, ``OpenFAST-v1.0.0-123-gabcd1234-dirty`` describes OpenFAST as:

=================== =============
 Version Component   Explanation
=================== =============
 v1.0.0              MAJOR.MINOR.PATCH numbering system; corresponds to a tagged commit made by NREL on GitHub
 123-g               Number of additional commits after the most recent tag for a build (the ``-g`` is for ``git``)
 abcd1234            First 8 characters of the current commit hash
 dirty               Denotes that local changes have been made but not committed; omitted if there are no local changes
=================== =============
