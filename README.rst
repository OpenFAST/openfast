OpenFAST
========

|actions| |nbsp| |rtfd|

.. |actions| image:: https://github.com/openfast/openfast/actions/workflows/automated-dev-tests.yml/badge.svg?branch=dev
   :target: https://github.com/OpenFAST/openfast/actions/workflows/automated-dev-tests.yml?query=workflow%3A%22Development+Pipeline%22
   :alt: Build Status
.. |rtfd| image:: https://readthedocs.org/projects/openfast/badge/?version=dev
   :target: https://openfast.readthedocs.io/en/dev
   :alt: Documentation Status
.. |nbsp| unicode:: 0xA0
   :trim:

OpenFAST is a wind turbine simulation tool which builds on FAST v8. FAST.Farm
extends the capability of OpenFAST to simulate multi-turbine wind farms. They were
created with the goal of being community models developed and used by research
laboratories, academia, and industry. They are managed by a dedicated team at the
National Renewable Energy Lab. Our objective is to ensure that OpenFAST and FAST.Farm
are sustainable software that are well tested and well documented. If you'd like
to contribute, see the `Developer Documentation <https://openfast.readthedocs.io/en/dev/source/dev/index.html>`_
and any open GitHub issues with the
`Help Wanted <https://github.com/OpenFAST/openfast/issues?q=is%3Aopen+is%3Aissue+label%3A"Help+wanted">`_
tag.

**OpenFAST is under active development**.

FAST v8 - OpenFAST
------------------
The transition from FAST v8 to OpenFAST represents the effort to better
support an open-source developer community around FAST-based aero-hydro-servo-
elastic engineering models of wind-turbines and wind-plants. OpenFAST is the
next generation of FAST analysis tools. More information is available in the
`transition notes <http://openfast.readthedocs.io/en/latest/source/user/fast_to_openfast.html>`_.

FAST v8, now OpenFAST, is a physics-based engineering tool for simulating the coupled dynamic
response of wind turbines. OpenFAST joins aerodynamics models, hydrodynamics models
for offshore structures, control and electrical system (servo) dynamics models,
and structural (elastic) dynamics models to enable coupled nonlinear aero-
hydro-servo-elastic simulation in the time domain. The OpenFAST tool enables the
analysis of a range of wind turbine configurations, including two- or
three-blade horizontal-axis rotor, pitch or stall regulation, rigid or
teetering hub, upwind or downwind rotor, and lattice or tubular tower. The wind
turbine can be modeled on land or offshore on fixed-bottom or floating
substructures. OpenFAST is based on advanced engineering models derived from
fundamental laws, but with appropriate simplifications and assumptions, and
supplemented where applicable with computational solutions and test data.

With OpenFAST, you can run large numbers of nonlinear time-domain simulations
in approximately real time to enable standards-based loads analysis for predicting
wind system ultimate and fatigue loads. You can also linearize the underlying
nonlinear model about an operating point to understand the system response
and enable the calculation of natural frequencies, damping, and mode shapes;
the design of controllers, and analysis of aero-elastic instabilities.

The aerodynamic models use wind-inflow data and solve for the rotor-wake
effects and blade-element aerodynamic loads, including dynamic stall. The
hydrodynamics models simulate the regular or irregular incident waves and
currents and solve for the hydrostatic, radiation, diffraction, and viscous
loads on the offshore substructure. The control and electrical system models
simulate the controller logic, sensors, and actuators of the blade-pitch,
generator-torque, nacelle-yaw, and other control devices, as well as the
generator and power-converter components of the electrical drive. The
structural-dynamics models apply the control and electrical system
reactions, apply the aerodynamic and hydrodynamic loads, adds gravitational
loads, and simulate the elasticity of the rotor, drivetrain, and support
structure. Coupling between all models is achieved through a modular
interface and coupler (glue code).

FAST.Farm extends the capabilities of OpenFAST to provide physics-based
engineering simulation of multi-turbine land-based, fixed-bottom offshore,
and floating offshore wind farms. With FAST.Farm, you can simulate each wind
turbine in the farm with an OpenFAST model and capture the relevant
physics for prediction of wind farm power performance and structural loads,
including wind farm-wide ambient wind, super controller, and wake advection,
meandering, and merging. FAST.Farm maintains computational efficiency
through parallelization to enable loads analysis for predicting the ultimate
and fatigue loads of each wind turbine in the farm.


Documentation
-------------
The full documentation is available at http://openfast.readthedocs.io/.

This documentation is stored and maintained alongside the source code.
It is compiled into HTML with Sphinx and is tied to a particular version
of OpenFAST. `Readthedocs <http://openfast.readthedocs.io>`_ hosts the following
versions of the documentation:

* ``latest`` - The latest commit on the ``main`` branch
* ``stable`` - Corresponds to the last tagged release
* ``dev`` - The latest commit on the ``dev`` branch

These can be toggled with the ``v: latest`` button in the lower left corner of
the docs site.

Obtaining OpenFAST and FAST.Farm
--------------------------------
OpenFAST and FAST.Farm are hosted entirely on GitHub so you are in the
`right place <https://github.com/OpenFAST/OpenFAST>`_!
The repository is structured with two branches following the
"git-flow" convention:

* ``main``
* ``dev``

The ``main`` branch is stable, well tested, and represents the most up to
date released versions of OpenFAST and FAST.Farm. The latest commit on ``main``
contains a tag with version info and brief release notes. The tag history can be
obtained with the ``git tag`` command and viewed in more detail on
`GitHub Releases <https://github.com/OpenFAST/openfast/releases>`_. For general
use, the ``main`` branch is highly recommended.

The ``dev`` branch is generally stable and tested, but not static. It contains
new features, bug fixes, and documentation updates that have not been compiled
into a production release. Before proceeding with new development, it is
recommended to explore the ``dev`` branch. This branch is updated regularly
through pull requests, so be sure to ``git fetch`` often and check
`outstanding pull requests <https://github.com/OpenFAST/openfast/pulls>`_.

For those not familiar with git and GitHub, there are many resources:

* https://guides.github.com
* https://try.github.io
* https://help.github.com/categories/bootcamp/
* https://desktop.github.com/
* http://nvie.com/posts/a-successful-git-branching-model/

Compilation, Usage, and Development
-----------------------------------
Details for compiling
`compiling <http://openfast.readthedocs.io/en/latest/source/install/index.html>`_,
`using <http://openfast.readthedocs.io/en/latest/source/user/index.html>`_, and
`developing <http://openfast.readthedocs.io/en/latest/source/dev/index.html>`_
OpenFAST and FAST.Farm on Unix-based and Windows machines are available at
`readthedocs <http://openfast.readthedocs.io>`_.

Help
----
Please use `GitHub Issues <https://github.com/OpenFAST/OpenFAST/issues>`_ to:

* ask usage questions
* report bugs
* request code enhancements

Users and developers may also be interested in the NREL National Wind
Technology Center (NWTC) `phpBB Forum <https://wind.nrel.gov/forum/wind/>`_,
which is still maintained and has a long history of FAST-related questions
and answers.

Acknowledgments
---------------

OpenFAST and FAST.Farm are maintained and developed by researchers and software
engineers at the `National Renewable Energy Laboratory <http://www.nrel.gov/>`_
(NREL), with support from the US Department of Energy's Wind Energy Technology
Office. NREL gratefully acknowledges development contributions from the following
organizations:

* Envision Energy USA, Ltd
* Brigham Young University
* The University of Massachusetts
* `Intel® Parallel Computing Center (IPCC) <https://software.intel.com/en-us/ipcc>`_
