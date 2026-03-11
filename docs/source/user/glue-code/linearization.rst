.. _glue-code-linearization:

Linearization
=============

OpenFAST can linearise the full multi-physics system about a periodic (or
static) operating point to produce continuous-time, first-order state-space
matrices of the form

.. math::

   \dot{\mathbf{x}} &= A\,\mathbf{x} + B\,\mathbf{u} \\
   \mathbf{y}       &= C\,\mathbf{x} + D\,\mathbf{u}

together with the coupling matrices *dUdu* (input-to-input feed-through) and
*dUdy* (output-to-input coupling).  The linearization engine lives in
``modules/openfast-library/src/FAST_ModGlue.f90``.

.. contents::
   :local:
   :depth: 2

User inputs for linearization
------------------------------

The following parameters appear in the main OpenFAST input file (``*.fst``)
under the **Linearization** section.

.. list-table::
   :header-rows: 1
   :widths: 25 12 63

   * - Parameter
     - Type
     - Description
   * - ``Linearize``
     - logical
     - Master switch.  Set to ``True`` to enable all linearization
       functionality.  When ``False`` all other linearization parameters are
       ignored.
   * - ``CalcSteady``
     - logical
     - When ``True``, OpenFAST first runs the simulation forward until the
       outputs at each target azimuth converge from one rotor revolution to the
       next (steady-state trimming), then performs linearization at each
       azimuth.  When ``False``, linearization is performed at user-specified
       absolute simulation times (``LinTimes``).
   * - ``TrimCase``
     - integer
     - Controller degree of freedom trimmed during ``CalcSteady`` to achieve
       periodic steady state.

       * ``1`` – yaw
       * ``2`` – generator torque
       * ``3`` – collective blade pitch
   * - ``TrimTol``
     - real
     - RMS convergence tolerance on normalised output error across one
       rotor revolution.  Trimming stops when the error falls below this
       value.   Typical value: ``1.0e-5``.
   * - ``TrimGain``
     - real
     - Proportional gain used by the built-in trim controller.
       Units are rad/(rad/s) for yaw/pitch cases and N·m/(rad/s) for
       the torque case.
   * - ``Twr_Kdmp``
     - real
     - Artificial tower damping coefficient (N/(m/s)) added during the
       ``CalcSteady`` run to help damp transients and reach steady state
       faster.  Set to 0 to disable.
   * - ``Bld_Kdmp``
     - real
     - Artificial blade damping coefficient (N/(m/s)) during ``CalcSteady``.
   * - ``NLinTimes``
     - integer
     - Number of linearization time points per rotor revolution (or number of
       equally spaced absolute time instants when ``CalcSteady=False``).
       Must be ≥ 1.  For a periodic model at least 12 azimuths are typically
       needed to resolve the per-revolution variation.
   * - ``LinTimes``
     - real array
     - Absolute simulation times (seconds) at which to linearise when
       ``CalcSteady=False``.  Length must equal ``NLinTimes``.  Ignored when
       ``CalcSteady=True``.
   * - ``LinInputs``
     - integer
     - Controls which input variables appear in the **B** and **D** matrices.

       * ``0`` (``LIN_NONE``) – no inputs; produces state matrix only.
       * ``1`` (``LIN_STANDARD``) – inputs flagged ``VF_Linearize`` by the
         module (default set by each module's ``InitVars``).
       * ``2`` (``LIN_ALL``) – all module inputs including debug ones.
   * - ``LinOutputs``
     - integer
     - Controls which output variables appear in the **C** and **D** matrices.

       * ``0`` (``LIN_NONE``) – no outputs.
       * ``1`` (``LIN_STANDARD``) – ``WriteOutput`` channels only
         (``VF_WriteOut`` flag).
       * ``2`` (``LIN_ALL``) – all module outputs.
   * - ``LinOutJac``
     - logical
     - When ``True`` (requires ``LinInputs=LinOutputs=2``), the full module
       Jacobian matrices are written to the linearization output file for
       debugging.
   * - ``LinOutMod``
     - logical
     - When ``True``, per-module ``.lin`` files are written in addition to the
       full-system file.

Module support for linearization
----------------------------------

Modules that appear in the linearization variable ordering (set in
``ModGlue_Init``) are:

InflowWind → SeaState → ServoDyn → ElastoDyn → BeamDyn → AeroDyn →
HydroDyn → SubDyn → MAP++ → MoorDyn

A module that is not in this ordered list causes a fatal error if
``Linearize=True``.

Variable selection
------------------

During ``ModGlue_Init``, the ``VF_Linearize`` flag is applied to variables
according to the ``LinInputs`` and ``LinOutputs`` settings:

* **States (x)**: the ``VF_Linearize`` flag is always set on all continuous
  state variables of every participating module.
* **Inputs (u)**:

  * ``LIN_NONE`` → flag cleared on all input variables.
  * ``LIN_STANDARD`` → keeps whatever ``VF_Linearize`` flag was set in the
    module's ``InitVars``; module developers choose the *standard* input set.
  * ``LIN_ALL`` → flag set on all input variables.
  * Variables with ``VF_NoLin`` always have ``VF_Linearize`` cleared,
    regardless of the above setting.

* **Outputs (y)**:

  * ``LIN_NONE`` → flag cleared on all output variables.
  * ``LIN_STANDARD`` → flag set only on outputs that also carry ``VF_WriteOut``.
  * ``LIN_ALL`` → flag set on all output variables.
  * Variables with ``VF_NoLin`` are always excluded.

The combined variable set is assembled into a ``ModGlueType`` named ``Lin``
via ``ModGlue_CombineModules``.

Steady-state trimming (``CalcSteady``)
---------------------------------------

When ``CalcSteady=True``, ``ModGlue_CalcSteady`` is called at each time step
to detect periodicity:

1. The module outputs tagged ``VF_Linearize`` (excluding ``VF_WriteOut``) are
   collected into a buffer indexed by azimuth angle.
2. After each complete revolution the outputs at each of the ``NLinTimes``
   azimuth targets are compared against the previous revolution via the
   normalised RMS error:

   .. math::

      \varepsilon = \sqrt{\frac{1}{N} \sum_{i=1}^{N}
        \left(\frac{y_i^{\rm current} - y_i^{\rm previous}}{r_i}\right)^2}

   where :math:`r_i = \max(y_{i,\rm max} - y_{i,\rm min},\, 0.01)` is the
   output range from the current revolution (with a floor to avoid division 
   by near-zero).

3. When :math:`\varepsilon < \texttt{TrimTol}`, ``FoundSteady=True`` and
   linearization at all ``NLinTimes`` azimuths proceeds automatically.

4. If the simulation reaches within approximately two revolutions of ``TMax``
   without converging, a warning is issued and linearization is forced.

The azimuth interpolation between buffer samples uses the extrapolation
routines from ``MV_ExtrapInterp`` (supports constant, linear, and quadratic
schemes depending on the number of available samples).

linearization at an operating point
-------------------------------------

``ModGlue_Linearize_OP`` assembles the full-system matrices at a single
operating point (time / azimuth):

1. **Module Jacobians**: for each module,
   ``FAST_JacobianPInput`` and ``FAST_JacobianPContState`` are called to
   compute the per-module sub-matrices *dYdu*, *dXdu*, *dYdx*, *dXdx* by
   central-difference finite differentiation.  The perturbation magnitudes are
   taken from each variable's ``Perturb`` field (see :ref:`glue-code-modvar`).

2. **Operating point extraction**: ``FAST_GetOP`` packs the current states,
   inputs, and outputs into the linearization arrays
   (``ModGlue%Lin%x``, ``%u``, ``%y``).

3. **Coupling matrices**: the input-output coupling matrices *dUdu* and *dUdy*
   are assembled from the mesh-mapping Jacobians to account for the fact that
   some module inputs are functions of other modules' outputs.

4. **Full-system assembly**: the per-module sub-matrices are placed into the
   combined glue-level matrices using the ``iGlu`` index ranges stored in
   each ``ModVarType``.

5. **Output**: ``ModGlue_CalcWriteLinearMatrices`` writes the ``.lin`` file
   containing:

   * Operating point values (**x_op**, **u_op**, **y_op**)
   * linearization channel names (from ``LinNames``)
   * Derivative order indicators (``VF_DerivOrder1``, ``VF_DerivOrder2``)
   * Rotating-frame flags (``VF_RotFrame``)
   * Full-system matrices **A**, **B**, **C**, **D**, **dUdu**, **dUdy**
   * Per-module matrices (if ``LinOutMod=True``)
   * Full Jacobians (if ``LinOutJac=True``)

Output file format
-------------------

Each linearization call produces a file named
``<RootName>.<N>.lin`` where *N* is the linearization index (1 … ``NLinTimes``).
The file is a plain-text ASCII file that can be read by the
`openfast_io <https://github.com/OpenFAST/openfast_io>`_ Python library or the
`pyFAST <https://github.com/OpenFAST/pyFAST>`_ post-processing tools.

Key fields in the file header:

* ``Rotor_Speed`` – rotor speed at linearization time (RPM)
* ``Azimuth`` – blade-1 azimuth at linearization time (deg)

Variable naming conventions
----------------------------

In linearization output files each channel label follows the pattern:

``<ModAbbr> <MeshName> <Field> [, component [, node [, unit]]]``

Examples:

* ``ED BlPitch1, rad`` – ElastoDyn individual blade-1 pitch state
* ``AD B1N001Fx force, node 1, N`` – AeroDyn blade 1 node 1 X-force input
* ``BD_1 B1TipTDxr translation displacement, node 10, m`` – BeamDyn instance 1

Module developers should ensure that the ``Name`` argument to ``MV_AddVar`` /
``MV_AddMeshVar`` and the entries in ``LinNames`` follow this convention for
consistency with post-processing tools.

Module developer responsibilities
-----------------------------------

To participate in linearization a module must:

1. Call ``MV_AddVar`` / ``MV_AddMeshVar`` with appropriate ``VF_Linearize``
   flags and supply ``LinNames`` for all variables that may appear in the
   standard linearization set.

2. Implement ``<Mod>_JacobianPInput`` and ``<Mod>_JacobianPContState``
   subroutines (or supply analytical Jacobians through the registry).  The
   glue code calls these via the ``FAST_JacobianPInput`` /
   ``FAST_JacobianPContState`` wrappers in ``FAST_Funcs.f90``.

3. Implement ``<Mod>_GetOP`` (via the registry) to extract the operating-point
   values of states, inputs, and outputs.

4. Mark variables that should **not** participate in linearization with
   ``VF_NoLin``.

5. Mark variables in the rotating reference frame with ``VF_RotFrame`` so that
   multi-blade coordinate (MBC) transformations applied by post-processing
   tools are aware of these variables.
