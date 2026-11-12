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
   * - ``ConvTol``
     - real
     - Convergence tolerance.  The iteration stops when the average
       `L2`-norm of the Newton update vector falls below this value.
       Typical value: ``1.0e-4``.  Tighter tolerances increase
       computational cost but may be needed for stiff problems.
   * - ``MaxConvIter``
     - integer
     - Maximum number of Newton convergence iterations per time step before the
       solver declares convergence failure.  Typical value: **20**.
       With ``ModCoupling=2`` or ``1``, a fatal error is issued on failure;
       with ``ModCoupling=3`` the Jacobian is rebuilt first and the step is
       retried before a warning is emitted.
   * - ``AutoRelax``
     - bool
     - Adaptive under-relaxation for the tight-coupling iterative solver.
       If set to ``true``, OpenFAST will only use the user provided
       ``RelaxFactor`` in the first iteration of each predictor or corrector
       step. In subsequent iterations, OpenFAST will increase the relaxation
       factor from the previous iteration by a factor of 1.2 if converging
       (residual decreasing) or halve the relaxation factor if diverging
       (residual increasing). The relaxation factor is bounded between a
       minimum of 0.01 and a maximum of 0.8. If set to ``false``, OpenFAST
       will use a constant relaxation factor given by ``RelaxFactor``.
       Default is true.
   * - ``RelaxFactor``
     - real
     - Constant (if ``AutoRelax`` is ``False``) or initial (if ``AutoRelax``
       is ``True``) under-relaxation factor for the tight-coupling iterative
       solver. Must be a positive number less than or equal to **1.0**. If
       ``AutoRelax`` is ``False``, default value is **0.7**. If ``AutoRelax``
       is ``True``, default value is **0.3**.
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
factorization to improve conditioning.

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

Each convergence iteration solves the linear system

.. math::

   \mathbf{J}\,\Delta\mathbf{z} = \mathbf{r}

for the update :math:`\Delta\mathbf{z}` of the states and inputs, where
:math:`\mathbf{r}` is the residual (``XB``).  **J** is factored once, when it is
built in ``BuildJacobianTC``.  The factors are reused in every convergence
iteration until the update strategy triggers a rebuild, so the cost of each
iteration is the cost of applying the factors.

The factorization and solve are performed by the ``BlockTriSolve`` module of the
NWTC Library (``modules/nwtc-library/src/BlockTriSolve.f90``):

* ``BTS_Factor`` is called in ``BuildJacobianTC``;
* ``BTS_Solve`` is called in ``FAST_SolverStep``.

The method is a direct solver that exploits the structure of **J**: most solver
variables depend directly on only a few others, and much of the coupling runs in
one direction.  The solution agrees with a dense LU solve to round-off.  The
linear solver has no user inputs.  The IO Jacobian (``BuildJacobianIO``), used
for the initial and linearization input-output solve, is still factored and
solved with dense ``LAPACK_getrf``/``LAPACK_getrs``.

Block triangular ordering
^^^^^^^^^^^^^^^^^^^^^^^^^

Unknown *j* is a dependency of equation *i* when :math:`J_{ij} \neq 0` for
:math:`i \neq j`.  ``BTS_Factor`` builds this dependency graph from the exact
nonzero pattern of **J** and finds its strongly connected components with
Tarjan's algorithm.  Ordering the components so that each one comes after the
components it depends on gives a symmetric permutation **P** such that
:math:`\mathbf{P}\mathbf{J}\mathbf{P}^T` is **block lower triangular**.  Here
:math:`n = N_J` is the number of unknowns and *K* is the number of diagonal
blocks, one per strongly connected component:

.. code-block:: none

               block 1   block 2   block 3            block K
             +---------+---------+---------+-------+---------+
     block 1 |  J_11   |    0    |    0    |  ...  |    0    |
             +---------+---------+---------+-------+---------+
     block 2 |  J_21   |  J_22   |    0    |  ...  |    0    |
             +---------+---------+---------+-------+---------+
     block 3 |  J_31   |  J_32   |  J_33   |  ...  |    0    |
             +---------+---------+---------+-------+---------+
       ...   |  ...    |  ...    |  ...    |  ...  |   ...   |
             +---------+---------+---------+-------+---------+
     block K |  J_K1   |  J_K2   |  J_K3   |  ...  |  J_KK   |
             +---------+---------+---------+-------+---------+

* **Diagonal blocks** :math:`J_{bb}` are either:

  * a single unknown, whose diagonal entry can have any nonzero value; or
  * a group of unknowns that are mutually coupled, i.e. each depends directly or
    indirectly on all the others.  Such a group cannot be split further by
    reordering.

* **Blocks below the diagonal**, :math:`J_{bc}` with :math:`c < b`, contain
  arbitrary coupling to earlier blocks.  They are often sparse and may be
  empty.
* **Blocks above the diagonal** are exactly zero.  No entries are dropped.

Every square matrix has such an ordering, because the strongly connected
components always form an acyclic dependency graph.  The two extremes are:

* *K* = *n*: every block is a single unknown, the reordered matrix is lower
  triangular, and the solve is pure forward substitution;
* *K* = 1: all unknowns are coupled, and the method reduces to a dense LU
  factorization.

The load scaling by ``UJacSclFact`` is applied to matching rows and columns and
does not change the nonzero pattern, so it does not affect the block structure.

Structure of the TC Jacobian
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reordered TC Jacobian typically has three regions:

.. code-block:: none

            upstream U    core C                   downstream D
          +------------+------------------------+--------------+
          | \          |                        |              |
      U   |   \  J_UU  |           0            |       0      |
          |     \      |                        |              |
          +------------+------------------------+--------------+
          |            | [C1]                   |              |
          |            |  *   [C2]              |              |
      C   |    J_CU    |  *    *   ...          |       0      |
          |            |  *    *    *   [Cm]    |              |
          +------------+------------------------+--------------+
          |                                     | \            |
      D   |          [ J_DU    J_DC ]           |   \  J_DD    |
          |                                     |     \        |
          +-------------------------------------+--------------+

The diagonal ``\`` marks the lower-triangular regions, ``[Ck]`` the coupled
diagonal blocks, and ``*`` the coupling between them.  The regions are:

* **Upstream (U):** unknowns that do not depend on any core or downstream
  unknown.  Each is a single-unknown diagonal block, so :math:`J_{UU}` is lower
  triangular; it may be diagonal, or even the identity.
* **Core (C):** one or more coupled diagonal blocks :math:`[C_1], \dots, [C_m]`,
  separated by any single unknowns that fall between them, with coupling
  ``*`` between them.  The core typically contains the structural states of the
  TC modules together with the interface loads and motions that couple them.
  It is the only part of **J** that needs an LU factorization.
* **Downstream (D):** unknowns that depend on upstream or core unknowns, but on
  which no upstream or core unknown depends.  Each is a single-unknown diagonal
  block, so :math:`J_{DD}` is lower triangular; it is diagonal when the
  downstream unknowns do not depend on each other.

The sizes of the regions depend on the model.  The structure is detected
automatically from the nonzero pattern; nothing in the solver is specific to a
module or configuration.

Factorization
^^^^^^^^^^^^^

Once per Jacobian update, ``BTS_Factor``:

* factors each coupled block (:math:`[C_1], \dots, [C_m]`) with LAPACK
  ``getrf``, using partial pivoting within the block;
* stores each single unknown's diagonal value;
* stores the coupling of each row to earlier blocks.  It is kept either as exact
  nonzeros in compressed sparse row (CSR) form or, for large coupling panels,
  in compressed low-rank form (see :ref:`below <glue-code-solver-lowrank>`).

Solve: block forward substitution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With the residual :math:`\mathbf{r}` permuted into the same order, ``BTS_Solve``
computes, for :math:`b = 1, \dots, K`:

.. math::

   \mathbf{x}_b = J_{bb}^{-1}\Big(\mathbf{r}_b - \sum_{c < b} J_{bc}\,\mathbf{x}_c\Big)

No inverse is formed; :math:`J_{bb}^{-1}` denotes a solve with the diagonal
block.  The solve works at two levels:

* **Across blocks:** a single forward sweep from block 1 to block *K*.  Each
  block uses only the solutions of earlier blocks, because the blocks above the
  diagonal are exactly zero.  No backward sweep across blocks is needed.
* **Within a diagonal block:**

  * for a single unknown, a division by its diagonal entry;
  * for a coupled block, a solve with the block's LU factors (LAPACK
    ``getrs``).  This is a forward substitution through *L* followed by a
    backward substitution through *U*, both limited to the unknowns of that
    block.

In terms of the three regions above, the sweep is:

.. math::

   \mathbf{x}_U &= J_{UU}^{-1}\,\mathbf{r}_U \\
   \mathbf{x}_C &= J_{CC}^{-1}\big(\mathbf{r}_C - J_{CU}\,\mathbf{x}_U\big) \\
   \mathbf{x}_D &= J_{DD}^{-1}\big(\mathbf{r}_D - J_{DU}\,\mathbf{x}_U
                   - J_{DC}\,\mathbf{x}_C\big)

* :math:`J_{UU}` and :math:`J_{DD}` consist only of single-unknown blocks, so
  applying :math:`J_{UU}^{-1}` and :math:`J_{DD}^{-1}` is ordinary forward
  substitution, one unknown at a time.
* Applying :math:`J_{CC}^{-1}` continues the same forward sweep through the
  blocks of the core.  For each coupled block :math:`[C_k]`, the coupling ``*``
  to earlier core blocks is subtracted and the block is then solved with its LU
  factors, forward and backward.  Any single unknowns between coupled blocks
  are divisions.

Backward substitution therefore occurs only inside coupled blocks, never across
blocks.  The solution is finally permuted back to the original order.

Numerical properties:

* The blocks above the diagonal are exactly zero, so no pivoting across blocks
  is needed and no Schur complement is formed.  Block forward substitution is
  backward stable when each diagonal-block solve is.
* The upstream and downstream unknowns never enter an LU factorization.  If the
  coupled blocks :math:`[C_1], \dots, [C_m]` contain :math:`n_1, \dots, n_m`
  unknowns, the factorization cost is proportional to
  :math:`\sum_k n_k^3`, instead of :math:`n^3` for a dense factorization of all
  :math:`n` unknowns.  Likewise, the cost of the LU solves in each iteration is
  proportional to :math:`\sum_k n_k^2` instead of :math:`n^2`.
* Diagonal entries do not affect the block structure.  A zero diagonal entry
  inside a coupled block is handled by pivoting within the block.  The finest
  block structure is obtained when the diagonal of **J** has no zeros, as is the
  case for the TC Jacobian.
* ``BTS_Solve`` does not modify the factorization.

.. _glue-code-solver-lowrank:

Coupling panels and low-rank compression
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once the LU factorization is limited to the coupled blocks, most of the
remaining work in each solve is applying the coupling to earlier unknowns:
:math:`J_{CU}`, the ``*`` entries in the core, and :math:`[J_{DU}\; J_{DC}]`.
This coupling is applied panel by panel:

* **Run:** a sequence of consecutive single-unknown blocks, in solve order,
  whose rows do not reference each other.  Every row in a run depends only on
  unknowns solved before the run starts.
* **Panel:** a group of rows whose coupling to earlier unknowns is applied
  together, as one matrix-vector product, before their diagonal blocks are
  solved.  A panel is either a run or the rows of one coupled block.  The
  coupling inside a coupled block belongs to its diagonal block and is handled
  by its LU factors.
* **Panel coupling matrix:** the horizontal strip of the matrix to the left
  of the panel's diagonal block(s).  It contains the panel's rows and only the
  earlier columns those rows reference.  It can be applied in one step because
  every unknown it references is solved before any row of the panel.

In the matrix above:

* each coupled block :math:`[C_k]` is a panel.  Its coupling matrix is its rows
  of :math:`J_{CU}` together with the ``*`` entries to the left of
  :math:`[C_k]`;
* the downstream region forms one or more runs.  When the downstream unknowns do
  not depend on each other (:math:`J_{DD}` diagonal), the entire region is a
  single run, with panel coupling matrix :math:`[J_{DU}\; J_{DC}]`.  Otherwise a
  new run starts at each row that references an earlier row of the current run.

Coupling matrices between modules are often nearly full but of low numerical
rank.  The coupling between the substructure and the hydrodynamic mesh
variables is one example.  Each panel coupling matrix :math:`\mathbf{G}` that
meets all of the following:

* at least 4096 nonzeros;
* at least 16 rows;
* at least 16 columns

is gathered into a dense matrix.  It is compressed with a truncated QR
factorization with column pivoting, :math:`\mathbf{G} \approx \mathbf{Q}\mathbf{W}`:

* **Q** has orthonormal columns (rows × rank);
* **W** is rank × columns;
* the factorization stops at the smallest rank whose discarded part meets the
  tolerance.

The panel is then applied with two matrix-vector products (BLAS ``gemv``):
first :math:`\mathbf{W}\mathbf{x}`, then the product with **Q**.  The compressed
form is used only if all of the following hold:

* :math:`\|\mathbf{G} - \mathbf{Q}\mathbf{W}\|_F \le 10^{-13}\,\|\mathbf{G}\|_F`,
  verified after compression by rebuilding **G** from **Q** and **W**;
* the rank is at most one eighth of the smaller dimension of **G**;
* the estimated solve cost is at least halved.

Otherwise the exact nonzeros are kept.  The tolerance is a fixed internal
constant.  It is several orders of magnitude below the accuracy of the
finite-difference Jacobian entries, and below the rounding error that a dense LU
factorization introduces for systems of this size.  Compression is therefore not
expected to affect convergence behavior.

Fallback to dense LU
^^^^^^^^^^^^^^^^^^^^

``BTS_Factor`` uses a dense LU factorization of the full **J**
(``LAPACK_getrf``/``LAPACK_getrs``) in the following cases:

* the system is small (fewer than 64 unknowns);
* the block form is not expected to be at least 25% cheaper per solve, e.g.
  when **J** is nearly full or consists of a single coupled block;
* a single-unknown diagonal entry is numerically zero relative to the largest
  entry in its row;
* a coupled diagonal block is exactly singular.

A singular **J** is reported by the dense factorization.  When the fallback is
triggered by a numerically zero pivot or a singular block, a one-time ``Info``
message is printed.  The parameter ``LinSolveMethod`` in ``FAST_Solver.f90`` can
be set to ``BTS_Method_Dense`` at compile time to always use the dense
factorization, for example for comparison.

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
