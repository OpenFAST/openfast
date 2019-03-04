C++ API Users Guide
===================

.. only:: html

   This document offers a quick reference guide for the C++ API and glue 
   code. It is intended to be used by the general user in combination
   with other OpenFAST manuals. The manual will be updated as new releases are
   issued and as needed to provide further information on advancements or
   modifications to the software.

The C++ API provides a high level API to run OpenFAST through a C++ gluecode. The primary purpose of the C++ API is to help interface OpenFAST to external programs like CFD solvers that are typically written in C++. The installation of C++ API is enabled via CMake by turning on the :cmakeval:`BUILD_OPENFAST_CPP_API` flag.

A sample glue-code `FAST_Prog.cpp <https://github.com/OpenFAST/openfast/blob/dev/glue-codes/openfast-cpp/src/FAST_Prog.cpp>`_ is provided as a demonstration of the usage of the C++ API. The glue-code allows for the simulation of multiple turbines using OpenFAST in parallel over multiple processors. The glue-code takes a single input file named ``cDriver.i`` (:download:`download <files/cDriver.i>`).

.. literalinclude:: files/cDriver.i
   :language: yaml

Command line invocation
-----------------------

.. code-block:: bash

   mpiexec -np <N> openfastcpp 
   
Common input file options
-------------------------

.. confval:: nTurbinesGlob

   Total number of turbines in the simulation. The input file must contain a number of turbine specific sections (`Turbine0`, `Turbine1`, ..., `Turbine(n-1)`) that is consistent with `nTurbinesGlob`.

.. confval:: debug
   
   Enable debug outputs if set to true

.. confval:: dryRun

   The simulation will not run if dryRun is set to true. However, the simulation will read the input files, allocate turbines to processors and prepare to run the individual turbine instances. This flag is useful to test the setup of the simulation before running it.
   
.. confval:: simStart

   Flag indicating whether the simulation starts from scratch or restart. ``simStart`` takes on one of three values:

   * ``init`` - Use this option when starting a simulation from `t=0s`.
   * ``trueRestart`` - While OpenFAST allows for restart of a turbine simulation, external components like the Bladed style controller may not. Use this option when all components of the simulation are known to restart.
   * ``restartDriverInitFAST`` - When the ``restartDriverInitFAST`` option is selected, the individual turbine models start from `t=0s` and run up to the specified restart time using the inflow data stored at the actuator nodes from a hdf5 file. The C++ API stores the inflow data at the actuator nodes in a hdf5 file at every OpenFAST time step and then reads it back when using this restart option. This restart option is especially useful when the glue code is a CFD solver.
   
.. confval:: tStart
   
   Start time of the simulation

.. confval:: tEnd

   End time of the simulation. tEnd <= tMax

.. confval:: tMax

   Max time of the simulation   

.. confval:: dtFAST

   Time step for FAST. All turbines should have the same time step.

.. confval:: nEveryCheckPoint

   Restart files will be written every so many time steps   

Turbine specific input options
------------------------------

.. confval:: turbine_base_pos

   The position of the turbine base for actuator-line simulations

.. confval:: num_force_pts_blade
   
   The number of actuator points along each blade for actuator-line simulations   
   
.. confval:: num_force_pts_tower

   The number of actuator points along the tower for actuator-line simulations.
  
.. confval:: restart_filename

   The checkpoint file for this turbine when restarting a simulation
   
.. confval:: FAST_input_filename

   The FAST input file for this turbine
  
.. confval:: turb_id

   A unique turbine id for each turbine
