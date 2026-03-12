.. _glue-code-overview:

Glue Code Overview
==================

The OpenFAST glue code is the software layer that connects individual physics modules
(ElastoDyn, AeroDyn, HydroDyn, ServoDyn, etc.) into a coupled simulation.  It lives
primarily in ``modules/openfast-library/src/`` and relies on the ``ModVar`` module in
``modules/nwtc-library/src/ModVar.f90`` for an abstract description of every
variable exchanged between modules.

High-level responsibilities include:

* Initialising each module and registering its variables with the glue code
* Managing multi-rate sub-stepping (modules whose time step is a divisor of the
  global time step)
* Mapping outputs of one module to inputs of another (motion meshes, load meshes,
  and scalar variables)
* Running the time-stepping loop under either loose coupling or tight
  generalized-alpha coupling
* Linearising the assembled system to produce state-space matrices

Source files
------------

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - File
     - Purpose
   * - ``FAST_Subs.f90``
     - Top-level initialisation: reads input file, calls each module's ``_Init``,
       and calls ``MV_AddModule`` to register every module with the glue code.
   * - ``FAST_ModGlue.f90``
     - Combines per-module variable descriptions into a monolithic ``ModGlueType``
       structure via ``ModGlue_CombineModules``; performs linearization
       (``ModGlue_Linearize_OP``) and steady-state trimming
       (``ModGlue_CalcSteady``).
   * - ``FAST_Solver.f90``
     - Implements the generalized-alpha tight-coupling solver
       (``FAST_SolverStep``), the input-output convergence loop
       (``FAST_CalcOutputsAndSolveForInputs``), and Jacobian assembly.
   * - ``FAST_Mapping.f90``
     - Mesh-to-mesh and variable-to-variable transfer mappings.
   * - ``FAST_Funcs.f90``
     - Wrappers around module-level ``CalcOutput``, ``UpdateStates``,
       ``CalcContStateDeriv``, ``GetOperatingPoint``, and ``SetOperatingPoint``
       that dispatch to the correct module instance.
   * - ``ModVar.f90`` (nwtc-library)
     - The ``ModVar`` module: data structures (``ModVarType``, ``ModVarsType``,
       ``ModDataType``, ``DatLoc``) and all ``MV_*`` subroutines.

Module coupling categories
--------------------------

Each module is assigned to exactly one coupling category during initialisation in
``FAST_SolverInit``:

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - Category
     - Flag
     - Description
   * - Tight Coupling (TC)
     - ``MC_Tight``
     - States and accelerations are solved simultaneously via the generalized-alpha
       Newton iteration.  ElastoDyn, BeamDyn, and SubDyn are tight-coupling modules
       when ``ModCoupling`` ≥ 2.
   * - Option 1
     - ``MC_Option1``
     - Modules whose inputs depend on TC outputs and are converged in the same
       Newton loop (e.g. HydroDyn, MoorDyn, ServoDyn with structural controllers).
   * - Option 2
     - ``MC_Option2``
     - Loosely coupled modules that are called once per step before the convergence
       loop (InflowWind, SeaState, AeroDyn, etc.).
   * - Post
     - ``MC_Post``
     - Modules whose input solve is deferred until after the convergence loop
       (ServoDyn, ExternalInflow).

Time-stepping loop (overview)
------------------------------

Each call to ``FAST_SolverStep`` follows this sequence:

1. **Correction iterations** (outer loop) – at most ``p%MaxConvIter`` iterations.
2. **Option 2** – input solve + state update + ``CalcOutput`` for loosely coupled modules.
3. **Option 1** – input solve + state update for semi-implicit modules.
4. **TC input solve** – gather inputs for the tight-coupling modules.
5. **Convergence iterations** (inner loop) – Newton-Raphson updates of TC
   states and inputs until the update norm falls below ``ConvTol`` or the
   iteration limit is reached.
6. **Post-solve input solves** – ServoDyn, ExternalInflow.

Module registration and variable ordering are described in detail in
:ref:`glue-code-modvar`.  How the per-module variables are assembled into
global arrays and Jacobian matrices is covered in :ref:`glue-code-modglue`.
The solver algorithm and Jacobian construction are covered in
:ref:`glue-code-solver`.
