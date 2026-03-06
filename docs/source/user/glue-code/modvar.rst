.. _glue-code-modvar:

Module Variables (``ModVar``)
=============================

The ``ModVar`` module (``modules/nwtc-library/src/ModVar.f90``) provides the data
structures and subroutines that allow each physics module to declare its
continuous states, inputs, and outputs in a way that the glue code can
manipulate generically—without knowing the internals of any particular module.

The complete type hierarchy—from ``DatLoc`` and ``ModVarType`` through
``ModVarsType`` and ``ModDataType`` up to the top-level ``ModGlueType``—is
illustrated in :numref:`fig-modvar-types`.

.. contents::
   :local:
   :depth: 2

Data structures
---------------

``DatLoc``
~~~~~~~~~~

A ``DatLoc`` (data location) is a small structure used to uniquely identify
where a particular variable lives inside a module's derived-type hierarchy:

.. code-block:: fortran

   TYPE :: DatLoc
     INTEGER :: Num   ! Mesh number (or 0 for non-mesh scalars)
     INTEGER :: i1    ! First index
     INTEGER :: i2    ! Second index
     INTEGER :: i3    ! Third index
     INTEGER :: i4
     INTEGER :: i5
   END TYPE DatLoc

A ``DatLoc`` value is created once per variable by the registry-generated
``DatLoc()`` constructor calls inside each module's ``InitVars`` subroutine and
passed to ``MV_AddVar`` / ``MV_AddMeshVar``.  The glue code stores this value
inside each ``ModVarType`` and uses ``MV_EqualDL`` to match variables across
source and destination modules when setting up mesh mappings.

``ModVarType``
~~~~~~~~~~~~~~

Describes a single variable (or group of variables for a mesh field):

.. code-block:: fortran

   TYPE :: ModVarType
     INTEGER      :: Field       ! Field type (FieldForce, FieldTransDisp, ...)
     INTEGER      :: Nodes       ! Number of nodes (mesh variables only)
     INTEGER      :: Num         ! Total number of scalar values
     INTEGER      :: Flags       ! Bit-mask of VF_* flags
     INTEGER      :: DerivOrder  ! 0=disp/orientation, 1=velocity, 2=acceleration
     INTEGER      :: iLoc(2)     ! [start, end] in module-local array
     INTEGER      :: iGlu(2)     ! [start, end] in glue-level array
     INTEGER      :: iq(2)       ! [start, end] in solver state (q) array
     INTEGER      :: iLB, iUB    ! User array bounds (for array-valued scalars)
     INTEGER      :: j, k, m, n  ! Additional user-defined indices
     REAL(R8Ki)   :: Perturb     ! Perturbation size for Jacobian finite differences
     TYPE(DatLoc) :: DL          ! Data location
     CHARACTER    :: Name        ! Human-readable variable name
     CHARACTER(:), ALLOCATABLE :: LinNames(:)  ! Per-value linearization labels
   END TYPE ModVarType

``ModVarsType``
~~~~~~~~~~~~~~~

Holds all variables for one module, partitioned into three groups:

.. code-block:: fortran

   TYPE :: ModVarsType
     INTEGER                             :: Nx, Nu, Ny
     TYPE(ModVarType), ALLOCATABLE :: x(:)   ! Continuous states
     TYPE(ModVarType), ALLOCATABLE :: u(:)   ! Inputs
     TYPE(ModVarType), ALLOCATABLE :: y(:)   ! Outputs
   END TYPE ModVarsType

``ModDataType``
~~~~~~~~~~~~~~~

Top-level container that the glue code holds for each module instance:

.. code-block:: fortran

   TYPE :: ModDataType
     CHARACTER    :: Abbr      ! Module abbreviation ("ED", "BD", ...)
     INTEGER      :: iMod      ! Index in glue module array
     INTEGER      :: ID        ! Module_ED, Module_BD, ...
     INTEGER      :: Ins       ! Instance number
     INTEGER      :: iRotor    ! Rotor index (0 = all rotors)
     INTEGER      :: SubSteps  ! Module sub-steps per solver step
     INTEGER      :: Category  ! Bit-mask of MC_* coupling flags
     REAL(R8Ki)   :: DT        ! Module time step
     TYPE(ModVarsType)  :: Vars
     TYPE(ModLinType)   :: Lin
   END TYPE ModDataType

Variable flags (``VF_*``)
-------------------------

Flags are combined via ``IOR`` and tested with ``MV_HasFlagsAll`` /
``MV_HasFlagsAny``.

.. list-table::
   :header-rows: 1
   :widths: 25 10 65

   * - Flag
     - Value
     - Meaning
   * - ``VF_None``
     - 0
     - No flags set; used as a wildcard that matches any variable.
   * - ``VF_Mesh``
     - 1
     - Variable is a mesh field (set automatically by ``MV_AddMeshVar``).
   * - ``VF_Line``
     - 2
     - Mesh is a *line* mesh (loads per unit length); linearization labels get
       ``/m`` suffix.
   * - ``VF_RotFrame``
     - 4
     - Variable lives in the rotating reference frame.
   * - ``VF_Linearize``
     - 8
     - Variable is included in the full-system linearization.
   * - ``VF_ExtLin``
     - 16
     - Variable is included in extended linearization output.
   * - ``VF_SmallAngle``
     - 32
     - Use small-angle approximation when computing orientation differences for
       linearization.
   * - ``VF_2PI``
     - 64
     - Scalar angle with range [0, 2π] (e.g. generator azimuth).
   * - ``VF_WriteOut``
     - 256
     - Output variable associated with a ``WriteOutput`` channel.
   * - ``VF_Solve``
     - 512
     - Variable participates in the tight-coupling Jacobian or input-output
       convergence solve.  Set automatically by ``FAST_SolverInit``;
       module developers should not set this flag manually.
   * - ``VF_AeroMap``
     - 1024
     - Variable used in aeromap computation.
   * - ``VF_Mapping``
     - 8192
     - Variable participates in a module-to-module transfer mapping.
   * - ``VF_NoLin``
     - 16384
     - Explicitly excludes a variable from both linearization and the solver
       (overrides ``VF_Linearize`` and ``VF_Solve``).

Field types (``Field*``)
------------------------

Used in the ``Field`` member of ``ModVarType`` and in ``MV_AddMeshVar``'s
``Fields`` argument:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Constant
     - Meaning
   * - ``FieldForce``
     - Nodal force (3 components per node, N)
   * - ``FieldMoment``
     - Nodal moment (3 components per node, N·m)
   * - ``FieldTransDisp``
     - Translational displacement (m)
   * - ``FieldOrientation``
     - Orientation, stored internally as unit-quaternion parameters (rad)
   * - ``FieldTransVel``
     - Translational velocity (m/s)
   * - ``FieldAngularVel``
     - Angular velocity (rad/s)
   * - ``FieldTransAcc``
     - Translational acceleration (m/s²)
   * - ``FieldAngularAcc``
     - Angular acceleration (rad/s²)
   * - ``FieldAngularDisp``
     - Angular displacement (rad)
   * - ``FieldScalar``
     - Generic scalar values

Convenience arrays defined in ``ModVar.f90``:

* ``LoadFields``   = ``[FieldForce, FieldMoment]``
* ``TransFields``  = ``[FieldTransDisp, FieldTransVel, FieldTransAcc]``
* ``AngularFields``= ``[FieldOrientation, FieldAngularVel, FieldAngularAcc, FieldAngularDisp]``
* ``MotionFields`` = all translational and angular motion fields

Adding module variables
-----------------------

Each module that participates in the glue code must implement an
``InitVars`` (or equivalent) subroutine that populates a ``ModVarsType``
structure by calling the ``MV_Add*`` subroutines documented below.  This
subroutine is called during the module's ``_Init`` routine and the resulting
``Vars`` is passed immediately to ``MV_AddModule``.

``MV_AddVar``
~~~~~~~~~~~~~

Adds a single (possibly multi-element) scalar variable to a variable array.

.. code-block:: fortran

   subroutine MV_AddVar(VarAry, Name, Field, DL, &
                        Num, iAry, jAry, kAry, &
                        Flags, DerivOrder, Perturb, LinNames, Active)

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Argument
     - Intent
     - Description
   * - ``VarAry``
     - ``INOUT``
     - Allocatable array of ``ModVarType``; the new variable is appended.
   * - ``Name``
     - ``IN``
     - Human-readable name used in debug output and linearization labels.
   * - ``Field``
     - ``IN``
     - Field type constant (``FieldScalar``, ``FieldTransDisp``, etc).
   * - ``DL``
     - ``IN``
     - ``DatLoc`` identifying where the data live in the module's derived type.
   * - ``Num``
     - ``IN`` (optional)
     - Number of scalar values.  Defaults to 1.  If 0, the call is a no-op.
   * - ``iAry``
     - ``IN`` (optional)
     - Starting lower-bound index if the data are stored in an array.
   * - ``jAry``, ``kAry``
     - ``IN`` (optional)
     - Second and third array indices (for 2-D or 3-D arrays).
   * - ``Flags``
     - ``IN`` (optional)
     - Initial ``VF_*`` flag bit-mask.  Defaults to ``VF_None``.
   * - ``DerivOrder``
     - ``IN`` (optional)
     - Override the automatically-inferred derivative order (0, 1, or 2).
   * - ``Perturb``
     - ``IN`` (optional)
     - Finite-difference perturbation magnitude.  A good default is roughly
       1 % of the expected variable magnitude.  For mesh fields the default
       computed inside ``ModVarType_Init`` is used when this is omitted.
   * - ``LinNames``
     - ``IN`` (optional)
     - Array of per-value linearization channel labels (length = ``Num``).
       **Required** for non-mesh, non-scalar variables.
   * - ``Active``
     - ``IN`` (optional)
     - Set to ``.false.`` to conditionally skip adding the variable.

**Example** – registering the generator torque input of ServoDyn:

.. code-block:: fortran

   call MV_AddVar(Vars%u, 'Generator torque command', FieldScalar, &
                  DatLoc(SrvD_u_GenTrq), &
                  Flags=VF_Linearize, &
                  Perturb=1.0e3_R8Ki, &        ! N*m
                  LinNames=['SrvD GenTrq, N*m'])

``MV_AddMeshVar``
~~~~~~~~~~~~~~~~~

Adds all requested mesh fields for a single ``MeshType`` to a variable array.
It is a convenience wrapper around ``MV_AddVar`` that iterates over the
``Fields`` argument and skips fields that are absent from the committed mesh.

.. code-block:: fortran

   subroutine MV_AddMeshVar(VarAry, Name, Fields, DL, Mesh, &
                            Flags, Perturbs, Active, iVar)

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Argument
     - Intent
     - Description
   * - ``VarAry``
     - ``INOUT``
     - Allocatable array to append the new variable(s) to.
   * - ``Name``
     - ``IN``
     - Base name for all fields on this mesh.
   * - ``Fields``
     - ``IN``
     - Integer array of field-type constants; use the ``LoadFields``,
       ``MotionFields``, etc. convenience parameters where appropriate.
   * - ``DL``
     - ``IN``
     - ``DatLoc`` for this mesh within the module's data type.
   * - ``Mesh``
     - ``INOUT``
     - The committed ``MeshType``.  Its ``ID`` field is set to identify it in
       mesh-mapping operations.  The subroutine returns without adding anything
       if the mesh has not been committed.
   * - ``Flags``
     - ``IN`` (optional)
     - Extra ``VF_*`` flags added on top of ``VF_Mesh``.
   * - ``Perturbs``
     - ``IN`` (optional)
     - Array of perturbation values, one per entry in ``Fields``.
   * - ``Active``
     - ``IN`` (optional)
     - Conditionally disable the entire mesh variable registration.
   * - ``iVar``
     - ``OUT`` (optional)
     - Returns the assigned mesh ``ID`` so the caller can store it for later
       field look-ups.

**Example** – registering ElastoDyn's blade-root output motion mesh:

.. code-block:: fortran

   call MV_AddMeshVar(Vars%y, 'BladeRootMotion', MotionFields, &
                      DatLoc(ED_y_BladeRootMotion, i), &   ! blade i
                      p%BladeRootMotion(i), &
                      Flags=VF_Linearize + VF_RotFrame)

``MV_AddModule``
~~~~~~~~~~~~~~~~

After a module's ``InitVars`` subroutine is complete, the caller registers the
module with the glue code using ``MV_AddModule``.

.. code-block:: fortran

   subroutine MV_AddModule(ModDataAry, ModID, ModAbbr, Instance, &
                           ModDT, SolverDT, Vars, Linearize, &
                           ErrStat, ErrMsg, iRotor)

.. list-table::
   :header-rows: 1
   :widths: 20 10 70

   * - Argument
     - Intent
     - Description
   * - ``ModDataAry``
     - ``INOUT``
     - Allocatable array of ``ModDataType``; the new entry is appended.
   * - ``ModID``
     - ``IN``
     - Module identifier constant (``Module_ED``, ``Module_BD``, etc.).
   * - ``ModAbbr``
     - ``IN``
     - Short abbreviation string used in output labels (``"ED"``, ``"BD"``).
   * - ``Instance``
     - ``IN``
     - Instance number (1-based).  Most modules have a single instance.
   * - ``ModDT``
     - ``IN``
     - Module time step (seconds).  Must be an exact integer divisor of
       ``SolverDT``.
   * - ``SolverDT``
     - ``IN``
     - Solver (global) time step.
   * - ``Vars``
     - ``IN``
     - Populated ``ModVarsType`` from the module's ``InitVars`` call.
   * - ``Linearize``
     - ``IN``
     - Whether linearization is enabled.  When ``.false.``, ``LinNames``
       arrays are deallocated to save memory.
   * - ``ErrStat`` / ``ErrMsg``
     - ``OUT``
     - Error status and message.
   * - ``iRotor``
     - ``IN`` (optional)
     - Rotor number for multi-rotor turbines (0 = all rotors).

**Sub-stepping logic**: if ``ModDT < SolverDT``, ``MV_AddModule`` calculates
``ModData%SubSteps = NINT(SolverDT/ModDT)`` and validates that the module DT
divides the solver DT exactly.  An error is returned if ``ModDT > SolverDT``.

**Typical call sequence inside FAST_Subs.f90**:

.. code-block:: fortran

   ! Module computes its own Vars in Init
   call ED_Init(InitInp, u, p, ..., InitOut, ErrStat, ErrMsg)
   ! Register with glue code
   call MV_AddModule(m%ModData, Module_ED, 'ED', 1, p%DT, p_FAST%DT, &
                     InitOut%Vars, p_FAST%Linearize, ErrStat, ErrMsg, iRotor=1)

``MV_InitVarsJac``
~~~~~~~~~~~~~~~~~~

Called inside each module's ``InitVars`` after all ``MV_AddVar`` /
``MV_AddMeshVar`` calls are complete.  It assigns the module-local ``iLoc``
index ranges to each variable and allocates the ``ModJacType`` working arrays
used during Jacobian calculations.

.. code-block:: fortran

   subroutine MV_InitVarsJac(Vars, Jac, Linearize, ErrStat, ErrMsg)

Perturbation values
-------------------

Every variable carries a ``Perturb`` value used for central-difference
finite-differencing when building module-level Jacobians.  The ``Perturb``
argument to ``MV_AddVar`` / ``MV_AddMeshVar`` should be chosen so that the
resulting output change is large enough to distinguish from numerical noise
but small enough to stay in the linear regime.  Typical values:

* Translational displacement: ``1.0e-4`` m
* Rotational (orientation): ``2.0e-5`` rad
* Translational velocity: ``1.0e-3`` m/s
* Angular velocity: ``2.0e-4`` rad/s
* Translational acceleration: ``1.0e-2`` m/s²
* Force: ``1.0e1`` N
* Moment: ``1.0e1`` N·m
* Generic scalar: context-dependent

The ``UJacSclFact`` input parameter (see :ref:`glue-code-solver-inputs`) is a
global conditioning factor that the solver applies to load variables in the
Jacobian to improve matrix conditioning when force/moment magnitudes are very
different from state magnitudes.

Orientation representation
--------------------------

Orientations are **not** stored or manipulated as direction cosine matrices
(DCMs) inside the glue-variable arrays.  Instead, a compact three-component
unit-quaternion parameterization is used:

.. math::

   \mathbf{q}_p = [q_1, q_2, q_3] \quad \text{where } q_0 = \sqrt{1 - q_1^2 - q_2^2 - q_3^2}

This parameterization avoids the redundancy in a full DCM and enables
straightforward finite-differencing via quaternion composition
(``quat_compose``).  Conversion utilities exported from ``ModVar`` include
``dcm_to_quat``, ``quat_to_dcm``, ``quat_compose``, ``quat_inv``,
``quat_to_rvec``, ``rvec_to_quat``, ``wm_to_quat``, and ``quat_to_wm``.

When computing orientation differences for Jacobian rows, ``MV_ComputeDiff``
computes the relative rotation between the negative and positive-perturbation
quaternions and converts it to a rotation vector (small-angle approximation
or full Rodrigues formula depending on ``VF_SmallAngle``).
