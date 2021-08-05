.. _hd-getting-started:

Installation and Getting Started
================================
HydroDyn is included in the OpenFAST software repository and consists
of two major components:

* `hydrodyn_driver` is the standalone HydroDyn executable
* `hydrodynlib` is the OpenFAST module library; it is most commonly
  used when driven through the HydroDyn driver or the OpenFAST glue code  

For installation instructions, see :ref:`installation`. In sections where
an installation target can be specific, use `hydrodyn_driver`.

Running the HydroDyn Driver
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The HydroDyn Driver has a simple command line interface:

.. code-block:: bash

    hydrodyn_driver <input_file>

where `input_file` is the file described in :ref:`hd-driver-input`.
Additional input files are required, including the :ref:`hd-primary-input`.
The time-series output as well as other output from HydroDyn are
described in :ref:`hd-output`.

Running HydroDyn coupled to OpenFAST
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To run an OpenFAST simulation with the HydroDyn module enabled, the
`CompHydro` flag must be switched on and the :ref:`hd-primary-input`
path supplied in the OpenFAST primary input file:

.. code-block::

    # In the "Feature switches" section
    1               CompHydro   - Compute hydrodynamic loads (switch) {0=None; 1=HydroDyn}

    # In the "Input files" section
    "HydroDyn.dat"  HydroFile   - Name of file containing hydrodynamic input parameters (quoted string)

The time-series output as well as other output from HydroDyn are
described in :ref:`hd-output`.
