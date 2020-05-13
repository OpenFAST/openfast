.. OpenFAST documentation master file, created by
   sphinx-quickstart on Wed Jan 25 13:52:07 2017.

OpenFAST Documentation
======================

.. only:: html

   :Version: |release|
   :Date: |today|

OpenFAST is a multi-physics, multi-fidelity tool for simulating the coupled
dynamic response of wind turbines. Practically speaking, OpenFAST is the
framework (or "glue code") that couples computational modules for
aerodynamics, hydrodynamics for offshore structures, control and electrical
system (servo) dynamics, and structural dynamics to enable coupled nonlinear
aero-hydro-servo-elastic simulation in the time domain. OpenFAST enables the
analysis of a range of wind turbine configurations, including two- or
three-blade horizontal-axis rotor, pitch or stall regulation, rigid or
teetering hub, upwind or downwind rotor, and lattice or tubular tower.
The wind turbine can be modeled on land or offshore on fixed-bottom or floating
substructures.

Established in 2017, OpenFAST is an open-source software package
that builds on FAST v8 (see :ref:`fast_to_openfast`). The glue code
and underlying modules are mostly written in Fortran
(adhering to the 2003 standard), and modules can also be written in C or
C++. It was created with the goal of being a community model developed
and used by research laboratories, academia, and industry. It is
managed by a dedicated team at the National Renewable Energy Lab.
Our objective is to ensure that OpenFAST is well tested, well
documented, and self-sustaining software. To that end, we are continually
improving the documentation and test coverage for existing code, and we
expect that new capabilities will include adequate testing and
documentation. If you'd like to contribute, see the
:ref:`dev_guide` and any open GitHub issues with the
`Help Wanted <https://github.com/OpenFAST/openfast/issues?q=is%3Aopen+is%3Aissue+label%3A"Help+wanted">`_
tag.

The following links provide more insight into OpenFAST as a software
package:

- `OpenFAST Github Organization <https://github.com/OpenFAST>`_
- `Github Repository <https://github.com/OpenFAST/OpenFAST>`_

**Documentation Directory**

.. toctree::
   :numbered:
   :maxdepth: 2

   source/this_doc.rst
   source/install/index.rst
   source/testing/index.rst
   source/user/index.rst
   source/dev/index.rst
   source/license.rst
   source/help.rst
   source/acknowledgements.rst
