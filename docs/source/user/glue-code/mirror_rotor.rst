.. _glue-code-mirror-rotor:

Mirrored (counter-clockwise) rotors
===================================

By convention an OpenFAST turbine rotor turns **clockwise when viewed from
upwind**.  Setting ``MirrorRotor = T`` in the main OpenFAST input file runs the
same turbine as its **mirror image**, so the rotor turns **counter-clockwise
viewed from upwind**, without changing any of the input files that describe it.

.. contents::
   :local:
   :depth: 2

.. _glue-code-mirror-rotor-input:

User input
----------

``MirrorRotor`` is read from the **Feature Switches and Flags** section of the
main OpenFAST input file (``*.fst``) and holds one flag per rotor:

.. code-block:: text

          1   NRotors         - Number of rotors in turbine (-)
          F   MirrorRotor     - Flag to reverse rotor rotation direction [1 to NRotors] {F=Normal, T=Mirror}

Nothing else changes.  The ElastoDyn, AeroDyn and blade input files, the airfoil
polars and the controller are all supplied exactly as they would be for the
clockwise machine, and the turbine described by those files is mirrored
internally.

Airfoil tables in particular are used **verbatim** — there is no polar
transformation, and the airfoil coordinate files are not modified.

.. _glue-code-mirror-rotor-concept:

What the flag does
------------------

Mirroring is the reflection :math:`S = \mathrm{diag}(1, -1, 1)` about the rotor
:math:`xz` plane.  Under that reflection

.. math::

   \mathbf{p}' = S\,\mathbf{p}, \qquad
   \mathbf{v}' = S\,\mathbf{v}, \qquad
   \boldsymbol{\omega}' = -S\,\boldsymbol{\omega}, \qquad
   R' = S\,R\,S

for positions, true vectors such as force and velocity, pseudovectors such as
moment and angular velocity, and direction cosine matrices respectively.

The reflection is applied **at module boundaries**, with one exception noted
below for the BeamDyn blade description.  The physics kernels —
the blade-element momentum solver, the unsteady aerodynamics and dynamic wake
models, the airfoil interpolation, and the structural finite elements — are
never told the rotor is mirrored.  They continue to solve the equivalent
clockwise problem, which is why the polars are used unchanged and why the
angle of attack, inflow angle and lift and drag coefficients come out
numerically identical to the clockwise machine.

Quantities crossing a boundary are converted on the way in and back on the way
out.  Inside ElastoDyn the azimuth and rotor speed states are the **physical**
ones, so a mirrored rotor really does have a negative shaft speed about the
:math:`+x` axis.

.. _glue-code-mirror-rotor-where:

Where the mirror is applied
---------------------------

.. list-table::
   :header-rows: 1
   :widths: 22 12 66

   * - Module
     - Changed
     - Where the transformation happens
   * - Glue code
     - yes
     - Reads ``MirrorRotor`` and distributes it to ElastoDyn, AeroDyn and
       BeamDyn; holds the restrictions listed below; and carries the whole
       ElastoDyn-to-ServoDyn presentation layer
   * - ServoDyn
     - **no**
     - Nothing internal.  Signals are converted in the glue code, and the
       generator and brake torque signs are applied inside ElastoDyn.  A
       Bladed-style controller is likewise untouched
   * - AeroDyn
     - yes
     - Four places, all at the edges — see the table below
   * - BeamDyn
     - yes
     - The blade **input data**, not the boundary.  See
       :ref:`glue-code-mirror-rotor-beamdyn`
   * - ElastoDyn
     - yes
     - Initial states, blade-pitch geometry, the gearbox, and the outputs
   * - SimplifiedElastoDyn
     - yes
     - The same pattern as ElastoDyn
   * - InflowWind
     - **no**
     - The flag mirrors the turbine, not the environment; a reflected
       turbulence box is supplied by the user
   * - BEMT, UnsteadyAero, DBEMT, AirfoilInfo
     - **no**
     - The physics kernels solve the equivalent clockwise problem
   * - AeroDyn driver
     - yes
     - Per-turbine flag; mirrors the prescribed hub kinematics and pitch
   * - SimplifiedElastoDyn driver
     - yes
     - Flag pass-through only

Within AeroDyn:

.. list-table::
   :header-rows: 1
   :widths: 34 46 20

   * - What
     - Where
     - Direction
   * - Blade twist and sweep
     - On read, after the cant angle is derived from the twist
     - In
   * - Rotor speed, blade pitch, toe angle, in-plane inflow, blade angular rate
     - ``SetInputsForBEMT``
     - In
   * - Skew-aligned disk frame
     - ``DiskAvgValues`` — the disk normal is a pseudovector and needs an
       explicit flip
     - In
   * - Hub and airfoil loads
     - Load conversion out of the blade-element solver
     - Out
   * - Output channels
     - Aerodynamic power, tangential force, and the lateral force, moment and
       induction coefficients
     - Out

.. _glue-code-mirror-rotor-beamdyn:

BeamDyn blades
--------------

BeamDyn is the one place where the mirror is **not** confined to a module
boundary.  A BeamDyn blade is described by a reference line and by full
:math:`6 \times 6` stiffness and mass matrices that couple bending, extension,
shear and torsion, and those cross-couplings carry a handedness of their own.
Presenting mirrored motion to an unmirrored blade would not give the mirrored
result.

The blade data is therefore transformed once, as it is read.  The key-point
:math:`y` coordinates and the structural twist are negated, and each matrix is
transformed by :math:`T\,M\,T` with

.. math::

   T = \mathrm{diag}(1, -1, 1, -1, 1, -1)

which is the reflection written in BeamDyn's ordering of three translational
followed by three rotational degrees of freedom.  In practice an entry changes
sign if exactly one of its two indices is a :math:`y` translation, an :math:`x`
rotation or a :math:`z` rotation.  The transform is its own inverse, and it
preserves the polar-inertia constraint that BeamDyn validates on input.

The blade input file itself still describes the **clockwise** blade and is
supplied unchanged, exactly as for ElastoDyn.  The finite-element solver is
still never told the rotor is mirrored.

.. _glue-code-mirror-rotor-conventions:

Output conventions
------------------

Two families of output channel behave differently, and the distinction matters
only for a mirrored rotor.

**Channels named after a rotor or drivetrain quantity** report it in the
**rotor's own convention** — positive when the rotor turns the way it was
designed to turn, whichever way that is.  ``RotSpeed`` is positive for a
normally operating rotor whether or not it is mirrored, and negative only when
the rotor is genuinely running backwards.

**Channels carrying an explicit axis suffix** (``*Mxa``, ``*Vxa``, ``*Axa``,
``*Pxa`` and their ``*xs`` counterparts) report the **physical** component about
that axis, so they change sign under the mirror.

Several names used to be aliases of a single value and are now separate
channels.  For a clockwise rotor the two are numerically identical, so no
existing model or output file changes.

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Rotor convention
     - Physical
     - Quantity
   * - ``RotSpeed``
     - ``LSSTipVxa``, ``LSSTipVxs``, ``LSSTipV``
     - Rotor angular speed
   * - ``RotAccel``
     - ``LSSTipAxa``, ``LSSTipAxs``, ``LSSTipA``
     - Rotor angular acceleration
   * - ``Azimuth``
     - ``LSSTipPxa``, ``LSSTipPxs``, ``LSSTipP``
     - Rotor azimuth
   * - ``RotTorq``, ``LSShftTq``
     - ``LSShftMxa``, ``LSShftMxs``, ``LSSGagMxa``, ``LSSGagMxs``
     - Low-speed shaft torque

``GenSpeed``, ``GenAccel``, ``HSShftV``, ``HSShftA``, ``HSShftTq`` and
``HSShftPwr`` follow the rotor convention, so the generator side reads positive
during normal operation regardless of rotation direction.  ``RotPwr`` and
``RotThrust`` are unchanged by the mirror: power is a product of two quantities
that both flip, and thrust is along the mirror axis.

.. _glue-code-mirror-rotor-signs:

Measured sign table
-------------------

The table below is **measured**, not asserted.  It is produced by running the
same model with ``MirrorRotor = F`` and ``MirrorRotor = T`` and comparing every
output channel.  Note that mirroring reverses the order in which the blades
sweep past a given point, so under non-axisymmetric inflow — shear, yaw, or
shaft tilt — blade 2 of the mirrored rotor corresponds to blade 3 of the
clockwise one.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Behaviour
     - Channels
   * - Identical
     - ``RotSpeed``, ``RotAccel``, ``Azimuth``, ``RotTorq``, ``LSShftTq``,
       ``GenSpeed``, ``GenAccel``, ``HSShftTq``, ``HSShftPwr``, ``RotPwr``,
       ``RotThrust``, ``LSShftFxa``, ``LSShftFza``, ``LSSTipMya``,
       ``YawBrFxp``, ``YawBrFzp``, ``YawBrMyp``, ``TwrBsFxt``, ``TwrBsMyt``,
       ``OoPDefl*``, ``TipDxc*``, ``TipDzc*``, ``RootFxc*``, ``RootFzc*``,
       ``RootMyc*``, ``BldPitch*``, ``GenTq``, ``GenPwr``, ``HSSBrTq``,
       the BeamDyn ``B*RootFxr``, ``B*RootFzr``, ``B*RootMyr``, ``B*TipTDxr``,
       ``B*TipTDzr``, ``B*TipRDyr`` families, and the AeroDyn ``*Alpha``,
       ``*Theta``, ``*Phi``, ``*Cl``, ``*Cd``, ``*Fn`` families
   * - Sign-flipped
     - ``LSShftMxa``, ``LSSTipVxa``, ``LSSTipAxa``, ``LSSGagMxa``,
       ``LSShftFya``, ``LSSTipMza``, ``YawBrFyp``, ``YawBrMxp``, ``YawBrMzp``,
       ``TwrBsFyt``, ``TwrBsMxt``, ``TwrBsMzt``, ``IPDefl*``, ``TipDyc*``,
       ``RootFyc*``, ``RootMxc*``, ``RootMzc*``,
       the BeamDyn ``B*RootFyr``, ``B*RootMxr``, ``B*RootMzr``, ``B*TipTDyr``,
       ``B*TipRDxr``, ``B*TipRDzr`` families, and the AeroDyn ``*Ft``,
       ``*Cy``, ``*Vindy`` families
   * - Mirrored angle
     - ``LSSTipPxa``, ``LSSGagPxa``

.. _glue-code-mirror-rotor-control:

Controllers
-----------

ServoDyn is not told that the rotor has been reversed.  Everything crossing its
boundary is presented in the **clockwise convention**, so an unmodified
controller — including a Bladed-style DLL such as ROSCO — sees exactly what it
would see on a clockwise machine and behaves identically.  Rotor speed, blade
pitch, generator and brake torque, shaft azimuth and the blade root moments are
all converted; the controller needs no mirrored copy and no new input.

Yaw is the exception, and it is deliberate.  Yaw acts about the vertical axis in
the inertial frame, so it is **not** a rotor-convention quantity and is left
alone.  The yaw angle, the wind direction and hence the yaw error all stay
physically correct, which means a yaw controller still points the nacelle into
the real wind rather than into its mirror image.

.. _glue-code-mirror-rotor-asymmetric:

Inputs that are not mirrored
----------------------------

The flag mirrors the **turbine**.  It does not mirror the environment, the
control setpoints, or any other input that happens to be one-sided.  Anything in
the list below is left exactly as written, and must be mirrored by hand if the
intent is to reproduce the mirror image of a clockwise simulation:

- initial or fixed nacelle yaw (``NacYaw``) and the neutral yaw position
  (``YawNeut``);
- wind direction and horizontal shear in a uniform wind file — vertical shear is
  symmetric about the mirror plane and needs no change;
- a full-field turbulence box, which has to be reflected in :math:`y`;
- prescribed force and moment time series for a structural control, where the
  lateral force and the roll and yaw moments change sign;
- lateral geometry such as ``NacCMyn``;
- the furl geometry of a furling turbine — the tail boom, tail fin and the
  rotor- and tail-furl axes are all offset to one side.  Setting
  ``MirrorRotor`` with ``Furling = True`` raises a warning and continues: the
  tail is modelled as a drag force applying a moment about the yaw axis, and
  that calculation does not depend on which way the rotor turns, since tail
  interaction with the wake is not modelled.  The combination has not been
  verified against a mirror pair, however, because doing so requires mirroring
  the furl input file as well.

This matters most when verifying the mirror: leaving one of these unmirrored
looks exactly like a sign error in the code.

.. _glue-code-mirror-rotor-limits:

Limitations
-----------

``MirrorRotor = T`` currently produces a fatal error when combined with any of
the following.  Each restriction is removed as that part of the code is worked
through.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Not yet supported with
     - Notes
   * - Linearization
     - ``Linearize = T``
   * - Steady-state solver
     - ``CompAeroMaps = T``
   * - SimplifiedElastoDyn
     - ``CompElast = 3``
   * - AeroDisk, ExtLoads
     - ``CompAero = 1`` or ``3``
   * - OLAF free vortex wake
     - ``Wake_Mod = 3`` in the AeroDyn input file
   * - AeroAcoustics
     - ``CompAA = True`` in the AeroDyn input file

.. _glue-code-mirror-rotor-verification:

Verification
------------

The mirror is verified by running a model twice, once clockwise and once
mirrored, and requiring **every** output channel to resolve to one of: identical,
exactly sign-flipped, a mirrored angle, or below the numerical noise floor.
Anything else indicates that two quantities have been combined while expressed in
different frames.

The check is repeated across a matrix of conditions, since any single condition
leaves most of the sign map untested — rigid and flexible blades, ElastoDyn and
BeamDyn blades, fixed and free drivetrain, vertical shear, positive and negative
nacelle yaw, fixed and free yaw, and combinations of those.  At the AeroDyn
module level the same comparison is run over blade pitch, wind speed, tip-speed
ratio, the propeller-brake state, shaft tilt, precone, both BEM models, dynamic
wake, and four unsteady-aerodynamic models.

Three of those comparisons are kept as regression cases, each paired with the
clockwise model it mirrors:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Case
     - What it covers
   * - ``5MW_Land_noDLL_Steady_MirrorRotor``
     - Steady wind, no controller, ElastoDyn blades
   * - ``5MW_Land_BD_noDLL_Steady_MirrorRotor``
     - The same with BeamDyn blades
   * - ``AWT_WSt_StartUp_HighSpShutDown_MirrorRotor``
     - The high-speed-shaft brake taking the rotor down through zero speed,
       which is the one torque signed by the direction of rotation
   * - ``5MW_Land_DLL_WTurb_MirrorRotor``
     - Turbulence and a Bladed-style controller, ElastoDyn blades
   * - ``5MW_Land_BD_DLL_WTurb_MirrorRotor``
     - The same with BeamDyn blades

The two turbulent cases are the ones that demonstrate the controller claim.  The
DISCON library is used completely unchanged, and blade pitch, generator torque,
generator power, generator speed and rotor speed all come out identical between
the clockwise and mirrored runs.  They read a ``y``-reflected copy of the
turbulence box, for the reason given above.

For a clockwise rotor every mirror-related expression reduces to a multiplication
by ``+1``, so existing regression baselines reproduce bit-for-bit.
