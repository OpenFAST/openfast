.. OpenFAST documentation master file, created by
   sphinx-quickstart on Wed Jan 25 13:52:07 2017.

OpenFAST Documentation
======================

.. only:: html

   :Version: |release|
   :Date: |today|

OpenFAST is an open-source wind turbine simulation tool that was established in
2017 with the FAST v8 code as its starting point (see :ref:`fast_to_openfast`).
OpenFAST is a multi-physics, multi-fidelity tool for simulating the coupled
dynamic response of wind turbines. Practically speaking, OpenFAST is the
framework (or glue code) that couples computational modules for
aerodynamics, hydrodynamics for offshore structures, control and electrical
system (servo) dynamics, and structural dynamics to enable coupled nonlinear
aero-hydro-servo-elastic simulation in the time domain. OpenFAST enables the
analysis of a range of wind turbine configurations, including two- or
three-blade horizontal-axis rotor, pitch or stall regulation, rigid or
teetering hub, upwind or downwind rotor, and lattice or tubular tower.
The wind turbine can be modeled on land or offshore on fixed-bottom or floating
substructures.

OpenFAST and its underlying modules are mostly written in Fortran (adhering to
the 2003 standard), but modules can be written in C/C++. OpenFAST was created
with the goal of being a community model, with developers and users from
research laboratories, academia, and industry. Our goal is also to ensure that
OpenFAST is sustainable software that is well tested and well documented. To
that end, we are continually improving the documentation and test coverage for
existing code, and we expect that new capabilities will include adequate
testing and documentation.

Here are some important links:
- `Github Organization Page <https://github.com/OpenFAST>`_
- `Github Repository <https://github.com/OpenFAST/OpenFAST>`_
- `Nightly Testing Results <http://my.cdash.org/index.php?project=OpenFAST>`_

.. toctree::
   :numbered:
   :maxdepth: 2

   source/this_doc.rst
   source/get_started.rst
   source/install/index.rst
   source/testing/index.rst
   source/user/index.rst
   source/dev/index.rst
   source/license.rst
   source/help.rst
   source/acknowledgements.rst

   Nightly Testing Results <http://my.cdash.org/index.php?project=OpenFAST>
   github.com Repository <https://github.com/openfast/openfast>
