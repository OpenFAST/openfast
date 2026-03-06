.. _glue-code-solver:

Solver
======

The OpenFAST tight-coupling solver is implemented in
``modules/openfast-library/src/FAST_Solver.f90``.  It integrates the continuous
states and resolves the input-output coupling between modules using a
generalized-alpha scheme with Newton-Raphson convergence iterations.

.. contents::
   :local:
   :depth: 2

.. _glue-code-solver-inputs:

User input parameters
---------------------

All solver parameters are set in the main OpenFAST input file
(``*.fst``) under the **Feature Switches and Flags** and
**Tight-Coupling / Solver** sections.

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Parameter
     - Type
     - Description
   * - ``DT``
     - real
     - Global (solver) time step in seconds.  All module time steps must be
       equal to or an integer sub-divisor of ``DT``.
   * - ``ModCoupling``
     - integer
     - Coupling method.

       * ``1`` – Loose coupling: structural modules (ED/BD/SD) are treated as
         Option 1 and do **not** participate in the tight Newton loop.
       * ``2`` – Tight coupling with fixed Jacobian updates (``DT_UJac``
         controls update frequency).
       * ``3`` – Tight coupling with adaptive Jacobian updates (the Jacobian is
         rebuilt whenever the Newton loop fails to converge within the
         iteration budget).
   * - ``RhoInf``
     - real
     - Numerical damping parameter ρ∞ for the generalized-alpha integrator.
       Range [0, 1]; 1 = no numerical damping (second-order accurate), 0 =
       maximum damping (first-order accurate).  Typical value: **0.9**.
       Reducing ``RhoInf`` below 1 damps high-frequency numerical noise at the
       cost of slightly reduced accuracy.
   * - ``MaxConvIter``
     - integer
     - Maximum number of Newton convergence iterations per time step before the
       solver declares convergence failure.  Typical value: **20**.
       With ``ModCoupling=2`` or ``1``, a fatal error is issued on failure;
       with ``ModCoupling=3`` the Jacobian is rebuilt first and the step is
       retried before a warning is emitted.
   * - ``ConvTol``
     - real
     - Convergence tolerance.  The iteration stops when the average
       `L2`-norm of the Newton update vector falls below this value.
       Typical value: ``1.0e-4``.  Tighter tolerances increase
       computational cost but may be needed for stiff problems.
   * - ``DT_UJac``
     - real
     - Time interval (seconds) between Jacobian rebuilds when
       ``ModCoupling=2``.

       * If ``DT_UJac < DT``: the Jacobian is rebuilt at a fraction of the
         convergence-iteration budget.
       * If ``DT_UJac ≥ DT``: the Jacobian is rebuilt every
         ``CEILING(DT_UJac/DT)`` time steps.
       * Setting ``DT_UJac`` very large (e.g. ``9999``) freezes the Jacobian
         for the entire simulation; useful for profiling or when the system
         is nearly linear and the Jacobian is expensive.
   * - ``UJacSclFact``
     - real
     - Conditioning scale factor applied to load rows and columns of the
       Jacobian.  Force and moment variables are divided by this factor before
       the linear solve and multiplied back afterwards, equalising the magnitude
       of load entries relative to displacement/velocity entries.  Typical
       value: **1.0e5** for offshore systems; may need adjustment for very
       large or very small turbines.
   * - ``CompElast``
     - integer
     - Select the structural dynamics module: ``1`` = ElastoDyn,
       ``2`` = BeamDyn (blades only, ElastoDyn still handles the tower/platform),
       ``3`` = Simplified ElastoDyn.  The chosen modules become TC members when
       ``ModCoupling ≥ 2``.
   * - ``CompSub``
     - integer
     - Sub-structural module: ``0`` = none, ``1`` = SubDyn, ``2`` = ExtPtfm,
       ``3`` = SlD (SoilDyn).  SubDyn joins the TC set when
       ``ModCoupling ≥ 2``.
   * - ``CompHydro``
     - integer
     - ``0`` = none, ``1`` = HydroDyn.  HydroDyn is always Option 1.
   * - ``CompMooring``
     - integer
     - ``0`` = none, ``1`` = MAP++, ``2`` = FEAMooring, ``3`` = MoorDyn,
       ``4`` = OrcaFlex.  Mooring modules are always Option 1.
   * - ``CompAero``
     - integer
     - Aerodynamics module: ``0`` = none, ``1`` = AeroDisk, ``2`` = AeroDyn.
       AeroDyn is Option 2 for land-based turbines and Option 1 for MHK.
   * - ``CompServo``
     - integer
     - Controller module: ``0`` = none, ``1`` = ServoDyn.  ServoDyn is
       Post-solve by default but becomes Option 1 when structural controllers
       (tower, blade, nacelle StC) are active.

Generalized-alpha integration
------------------------------

The tight-coupling solver integrates second-order ODEs of the form

.. math::

   \mathbf{M}\,\ddot{\mathbf{q}} + \mathbf{f}(\mathbf{q}, \dot{\mathbf{q}}, t) = 0

using the **generalized-alpha method** (Chung & Hulbert, 1993).  Given the
spectral radius ρ∞ specified by ``RhoInf``, the method parameters are:

.. math::

   \alpha_m &= \frac{2\rho_\infty - 1}{\rho_\infty + 1} \\
   \alpha_f &= \frac{\rho_\infty}{\rho_\infty + 1} \\
   \gamma   &= \tfrac{1}{2} - \alpha_m + \alpha_f \\
   \beta    &= \tfrac{1}{4}(1 - \alpha_m + \alpha_f)^2

Two derived coefficients used throughout the convergence loop are:

.. math::

   \beta'   &= h^2 \beta \frac{1 - \alpha_f}{1 - \alpha_m} \\
   \gamma'  &= h \gamma  \frac{1 - \alpha_f}{1 - \alpha_m}

where *h* = ``DT``.

**State vector layout** – the solver maintains a per-module *generalized
coordinate* (q) vector with four columns:

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Column
     - Meaning
   * - ``q``
     - Displacement / orientation states (``DerivOrder = 0``)
   * - ``v``
     - Velocity states (``DerivOrder = 1``)
   * - ``vd``
     - Acceleration (physical, from module ``CalcContStateDeriv``)
   * - ``a``
     - Algorithmic acceleration (generalized-alpha internal variable)

State prediction at the start of each step:

.. math::

   q_{n+1}^{\rm pred}  &= q_n + h v_n + h^2[(\tfrac{1}{2} - \beta)a_n + \beta\, a_{n+1}] \\
   v_{n+1}^{\rm pred}  &= v_n + h[(1-\gamma)a_n + \gamma\, a_{n+1}]

Module ordering
---------------

During ``FAST_SolverInit`` each module is categorised based on ``ModCoupling``
and its own physics type, and assigned to one of the ordered index arrays
in the ``Glue_TCParam`` structure:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Array
     - Modules (in order)
   * - ``iModTC``
     - ElastoDyn, BeamDyn, SubDyn (when ``ModCoupling ≥ 2``)
   * - ``iModOpt1``
     - ServoDyn (when StC active), SED, AD (MHK), ExtPtfm, HydroDyn, OrcaFlex,
       MoorDyn; ED/BD/SD also appear here when ``ModCoupling = 1``
   * - ``iModOpt2``
     - ServoDyn, SED, ED, BD, SD, InflowWind, SeaState, AeroDyn (land),
       AeroDisk, ExtLoads, MAP++, FEAMooring, IceDyn, IceFloe, SoilDyn
   * - ``iModPost``
     - ServoDyn, ExternalInflow
   * - ``iModInit``
     - SED, ED, BD, SD, InflowWind, ExtLoads (Step 0 initialisation only)

Jacobian construction
---------------------

Two separate Jacobians are assembled:

1. **TC/Option-1 Jacobian** (``BuildJacobianTC``) — for the main time-stepping
   convergence loop.
2. **IO Jacobian** (``BuildJacobianIO``) — for the initial and linearization
   input-output solve.

Variable selection (``VF_Solve`` flag)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

During ``FAST_SolverInit → SetVarSolveFlags``, the ``VF_Solve`` flag is set on
the variables that must appear in the Jacobian:

* **Continuous states** of all TC modules (automatically).
* **Motion mesh** inputs/outputs of TC-to-TC mappings (all fields).
* **Motion mesh** input accelerations of TC-to-Option1 or
  Option1-to-TC mappings.
* **Load mesh** inputs and outputs involved in any TC/Option1 mapping.
* **Load mesh** displacement outputs of the destination module when the mapping
  carries moments (needed for moment-arm Jacobian terms).
* **Variable-to-variable** mapped inputs/outputs of TC/Option1 modules.
* Any variable with ``VF_NoLin`` is excluded from ``VF_Solve``.

Jacobian structure (TC Jacobian)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The assembled TC Jacobian **J** has size ``NumJ × NumJ``, where:

.. math::

   N_J = \underbrace{N_Q}_{\text{TC states}} +
         \underbrace{N_{U_T}}_{\text{TC inputs}} +
         \underbrace{N_{U_1}}_{\text{Option-1 inputs}}

The columns and rows are partitioned as:

.. math::

   \mathbf{J} = \begin{bmatrix}
     J_{11} & J_{12} \\
     J_{21} & J_{22}
   \end{bmatrix}

where

* **J₁₁** (``NumQ × NumQ``) — derivative of the acceleration residual with
  respect to TC displacement/velocity states (formed from the module
  ``dXdx`` sub-Jacobians plus the generalized-alpha tangent).
* **J₁₂** (``NumQ × NumU_T``) — derivative of the acceleration residual with
  respect to TC inputs (from ``dXdu``).
* **J₂₁** (``NumU_T × NumQ``) — derivative of the input residual with
  respect to TC states (from ``dUdx = dUdy · dydx``).
* **J₂₂** (``NumU × NumU``) — derivative of the input residual with
  respect to inputs, including load conditioning rows/columns.

The right-hand side (XB) contains the residuals:

* **State residual** (rows ``iJX``): difference between the predicted
  velocity derivative and the module-computed accelerations.
* **Input residual** (rows ``iJU``): difference between the inputs computed
  from mesh mappings (``FAST_InputSolve``) and the current iterate.

The loads portion (rows ``iJL``) is pre-divided by ``UJacSclFact`` before the
factorisation to improve conditioning.

Jacobian update strategy
~~~~~~~~~~~~~~~~~~~~~~~~

``ModCoupling = 2`` (fixed updates)
   The Jacobian is rebuilt if either of these counters reaches zero:

   * ``UJacStepsRemain`` — steps remaining; initialised to
     ``CEILING(DT_UJac/DT)`` each time the Jacobian is rebuilt.
   * ``UJacIterRemain`` — iteration budget; initialised to
     ``CEILING(DT_UJac/DT · MaxConvIter)`` when ``DT_UJac < DT``.

   On convergence failure the solver returns a fatal error immediately.

``ModCoupling = 3`` (adaptive updates)
   The Jacobian is rebuilt the first time the convergence loop fails.  If the
   step still does not converge after the forced rebuild, a non-fatal warning
   is issued and the simulation proceeds.

Per-module Jacobian contributions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The module-level Jacobian sub-matrices are computed by finite differencing
inside ``BuildJacobianTC`` and ``BuildJacobianIO`` using the ``MV_Perturb`` /
``MV_ComputeDiff`` / ``MV_ComputeCentralDiff`` utilities from ``ModVar``.  For
each variable flagged ``VF_Solve``:

1. Apply a positive perturbation of magnitude ``Var%Perturb`` to the working
   state/input array.
2. Call ``FAST_CalcOutput`` (or ``FAST_GetContStateDeriv``).
3. Apply an equal negative perturbation.
4. Call again.
5. Compute the central difference: ``(y_plus - y_minus) / (2·Perturb)``.

For orientation variables (``FieldOrientation``), perturbations are applied by
quaternion composition rather than direct addition (``MV_Perturb``), and
differences are extracted as rotation vectors (``MV_ComputeDiff``).

Linear solve
~~~~~~~~~~~~

The LU factorisation of **J** is computed with ``LAPACK_getrf`` and the system
is solved with ``LAPACK_getrs`` (packed in ``NWTC_LAPACK``).  The same
factored matrix is reused across convergence iterations until the update
strategy triggers a rebuild.

Convergence check
~~~~~~~~~~~~~~~~~

After each Newton step the convergence error is the average `L2`-norm of the
update vector:

.. math::

   e = \frac{\|\Delta \mathbf{z}\|_2}{N_J}

where :math:`\Delta \mathbf{z}` combines state and input updates.  The loop
exits if ``e < ConvTol`` (``ErrID_None``) or the iteration count reaches
``MaxConvIter`` (``ErrID_Fatal`` or ``ErrID_Warn`` depending on
``ModCoupling``).

Output channels from the solver
--------------------------------

Three output channels are written to ``DriverWriteOutput`` each step and
appear in the output file when enabled:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Index
     - Content
   * - 1
     - Total convergence iterations in the step (``TotalIter``)
   * - 2
     - Final convergence error (``ConvError``)
   * - 3
     - Number of Jacobian rebuilds in the step (``NumUJac``)
