.. _glue-code-mirror-rotor:

Mirrored (counter-clockwise) rotors
===================================

By convention an OpenFAST turbine rotor turns **clockwise when viewed from
upwind**.  Setting ``MirrorRotor = T`` in the main OpenFAST input file makes the
rotor turn **counter-clockwise viewed from upwind** instead, without changing any
of the input files that describe it.

The purpose is to model a counter-clockwise rotor, not to mirror a whole turbine.
The reflection is applied internally to the rotor and to the drivetrain
quantities attached to it.  The tower, nacelle, support structure, mooring and
the environment are left alone, and are described exactly as they really are.

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
   * - AeroDisk
     - yes
     - The rotor speed and tip-speed ratio going into the coefficient table,
       the skew-aligned disk triad, and the moment components coming out.  The
       triad is built so that all three basis vectors mirror as true vectors,
       which makes the disk-frame force components invariant and every
       disk-frame moment component change sign
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

The table below is **measured**, not asserted.  It is generated by running every
registered mirrored case and comparing each channel against its clockwise
counterpart, so it records behaviour that is actually observed rather than
behaviour that is expected from reading the code.  The cases contributing to it
are the ElastoDyn + AeroDyn pair, the BeamDyn pair, the marine-turbine pair, the
yawed AeroDisk pair, the two-rotor aerodynamic driver case and the twin-rotor
semisubmersible.

Note that mirroring reverses the order in which the blades sweep past a given
point, so under non-axisymmetric inflow — shear, yaw, or shaft tilt — blade 2 of
the mirrored rotor corresponds to blade 3 of the clockwise one.  Mooring lines
exchange in the same way, and which line pairs with which depends on the layout.

It records how the two runs of a **symmetric** comparison relate to each other,
which is how the implementation is verified.  It is not a claim that these
channels change sign whenever the flag is set: in an ordinary simulation the
tower, support structure and inflow are not mirrored, so a quantity such as
``TwrBsMxt`` is simply the response of an unchanged structure to a
counter-clockwise rotor.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Behaviour
     - Channels
   * - Identical
     - ``RotSpeed``, ``RotAccel``, ``Azimuth``, ``RotTorq``, ``LSShftTq``,
       ``RotPwr``, ``RotThrust``, ``GenSpeed``, ``GenAccel``, ``GenTq``,
       ``GenPwr``, ``HSShftTq``, ``HSShftPwr``, ``LSShftFxa``, ``LSShftFza``,
       ``LSSTipMya``, ``LSSGagMya``, ``YawBrFxp``, ``YawBrFzp``, ``YawBrMyp``,
       ``TwrBsFxt``, ``TwrBsFzt``, ``TwrBsMyt``, ``TTDspFA``, ``TwrTpTDxi``,
       ``TwrTpTDzi``, ``NcIMUTAxs``, ``NcIMUTAzs``, ``PtfmSurge``,
       ``PtfmHeave``, ``PtfmPitch``, ``HydroFxi``, ``HydroFzi``,
       ``OoPDefl*``, ``TipDxc*``, ``RootFxc*``, ``RootFzc*``, ``RootMyc*``,
       ``RootFxb*``, ``RootFzb*``, ``RootMyb*``;
       the BeamDyn ``B*RootFxr``, ``B*RootFzr``, ``B*RootMyr``, ``B*TipTDxr``,
       ``B*TipTDzr``, ``B*TipRDyr``, ``B*FldFz`` families;
       the AeroDyn ``*Alpha``, ``*Theta``, ``*Phi``, ``*Fn``, ``*Fl``,
       ``*Fd``, ``*Fx``, ``*Cl``, ``*Cd``, ``*Cx``, ``*Cn``, ``*Vrel``,
       ``*Vindx``, ``*AxInd``, ``*TnInd``, ``*Gam``, ``*VUndx``,
       ``*VDisx`` families, in both the module (``B1N001Fn``) and nodal
       (``AB1N001Fn``) forms;
       the per-blade ``B*AeroPwr``;
       the AeroDisk ``ADFx``, ``ADFy``, ``ADFz``, ``ADFxi``, ``ADFzi``,
       ``ADMyi``, ``ADCp``, ``ADCq``, ``ADCt``, ``ADPower``, ``ADSkew``,
       ``ADTSR``, ``ADVRel`` channels;
       ``RtAeroFxh``, ``RtAeroFzh``, ``RtAeroMyh``, ``RtAeroPwr``,
       ``RtAeroCp``, ``RtAeroCt``, ``RtArea``, ``RtSkew``, ``RtTSR``,
       ``RtSpeed``, ``RtVAvgxh``, ``RtVAvgzh``, ``RtFldFxh``, ``RtFldFzg``,
       ``RtFldMyg``; and for a marine turbine the buoyant ``*Fbn``, ``*Fbs``,
       ``*Mbt``, ``HbFbx``, ``HbFbz``, ``HbMby`` families together with
       ``*SgCav``, ``*SigCr`` and ``*Clrnc``
   * - Sign-flipped
     - ``LSShftMxa``, ``LSShftFya``, ``LSSTipMza``, ``LSSGagMxa``,
       ``LSSGagMza``, ``LSSTipVxa``, ``LSSTipAxa``, ``YawBrFyp``,
       ``YawBrMxp``, ``YawBrMzp``, ``TwrBsFyt``, ``TwrBsMxt``, ``TwrBsMzt``,
       ``TTDspSS``, ``TwrTpTDyi``, ``NcIMUTAys``, ``PtfmSway``, ``PtfmRoll``,
       ``PtfmYaw``, ``HydroFyi``,
       ``IPDefl*``, ``TipDyc*``, ``RootFyc*``, ``RootMxc*``, ``RootMzc*``,
       ``RootFyb*``, ``RootMxb*``, ``RootMzb*``;
       the BeamDyn ``B*RootFyr``, ``B*RootMxr``, ``B*RootMzr``, ``B*TipTDyr``,
       ``B*TipRDxr``, ``B*TipRDzr``, ``B*FldMx`` families;
       the AeroDyn ``*Ft``, ``*Fy``, ``*Cy``, ``*Cm``, ``*Ct``, ``*Vindy``,
       ``*STVy``, ``*Mm`` families, again in both the module and nodal forms;
       the AeroDisk ``ADMx``, ``ADMy``, ``ADMz``, ``ADFyi``, ``ADMxi``,
       ``ADMzi``, ``ADSpeed`` channels;
       ``RtAeroFyh``, ``RtAeroMxh``, ``RtAeroMzh``, ``RtAeroCq``,
       ``RtVAvgyh``, ``RtFldFyh``, ``RtFldMxh``, ``RtFldMzh``; and for a
       marine turbine the buoyant ``*Fbt``, ``*Mbn``, ``*Mbs``, ``HbFby``,
       ``HbMbx``, ``HbMbz`` families
   * - Mirrored angle
     - ``LSSTipPxa``, ``LSSGagPxa``

Every channel that carries a mirror sign in the code is now requested by at
least one registered case, so no entry in this table is inferred from reading
the source.  The coefficient families in particular are measured together with
their mirror-invariant partners — ``*Cy``, ``*Cm`` and ``*Ct`` alongside
``*Cl``, ``*Cd``, ``*Cx`` and ``*Cn`` — so that a sign applied to a whole group
by mistake cannot pass unnoticed.

Across the six pairs, 1002 channels resolve: 677 identical, 324 sign-flipped and
one mirrored angle, with a further 58 below the noise floor in every case and
133 mooring channels set aside because their pairing is layout-specific.  No
channel is unresolved.


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

Nothing outside the rotor is mirrored
-------------------------------------

The flag reverses the **rotor**.  It does not touch the environment, the control
setpoints, the support structure, or anything else that happens to be one-sided.
For ordinary use that is exactly what is wanted: a counter-clockwise rotor, on
the turbine and in the conditions you actually have.  A real mooring spread, a
real wind field and a real yaw setpoint should all be left as they are.

The list below matters only when **verifying** the implementation.  That check
compares a clockwise run against a mirrored one and expects the two to be
reflections of each other, which requires the whole problem — not just the rotor
— to be symmetric about the rotor ``xz`` plane.  For that comparison, and only
for it, these have to be mirrored by hand:

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

Leaving one of these unmirrored during a verification run looks exactly like a
sign error in the code, which is the only reason the list is written down.

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
   * - ExtLoads
     - ``CompAero = 3``
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

Some of those comparisons are kept as regression cases, each paired with the
clockwise model it mirrors — either a case that already existed, or a ``_CW``
case registered alongside it.  They all carry the ctest label ``mirrorrotor``,
so ``ctest -R mirrorrotor`` runs the set:

.. list-table::
   :header-rows: 1
   :widths: 45 55

   * - Case
     - What it covers
   * - ``5MW_Land_noDLL_Steady_MirrorRotor``
     - Steady wind, no controller, ElastoDyn blades.  The baseline pair, pinned
       against ``5MW_Land_noDLL_Steady_CW``
   * - ``5MW_Land_BD_noDLL_Steady_MirrorRotor``
     - The same with BeamDyn blades, against ``5MW_Land_BD_noDLL_Steady_CW``
   * - ``AWT_WSt_StartUp_HighSpShutDown_MirrorRotor``
     - The high-speed-shaft brake taking the rotor down through zero speed,
       which is the one torque signed by the direction of rotation.  Also a
       two-bladed teetering rotor, which needs no blade swap
   * - ``5MW_Land_DLL_WTurb_MirrorRotor``
     - Turbulence and a Bladed-style controller, ElastoDyn blades
   * - ``5MW_Land_BD_DLL_WTurb_MirrorRotor``
     - The same with BeamDyn blades
   * - ``5MW_OC4Semi_WSt_WavesWN_MirrorRotor``
     - A floating platform, exercising HydroDyn, SeaState and MoorDyn beneath a
       mirrored rotor
   * - ``5MW_Land_DLL_WTurb_ADsk_SED_MirrorRotor``
     - SimplifiedElastoDyn and AeroDisk in place of ElastoDyn and AeroDyn
   * - ``5MW_MRSemi_DLL_WSt_WavesIrr_MirrorRotor``
     - A twin-rotor floating machine, with the whole stack solved together, so
       one rotor is mirrored and the other is not
   * - ``MHK_RM1_Floating_Steady_MirrorRotor``
     - A marine turbine, where the blade buoyancy and centre-of-buoyancy offset
       are mirrored too.  Paired with ``MHK_RM1_Floating_Steady_CW``
   * - ``ad_MultipleHAWT_MirrorRotor``
     - The AeroDyn driver rather than the glue code, covering the **nodal**
       output path with two rotors in a single run
   * - ``5MW_Land_ADsk_SED_Yaw_MirrorRotor``
     - A yawed AeroDisk rotor, paired with ``5MW_Land_ADsk_SED_Yaw_CW``.  The
       nacelle yaw is reversed between the two, since yaw acts about the
       inertial vertical and is not mirrored.  It uses a coefficient table with
       the lateral coefficients filled in, because the shipped 5MW table has
       them identically zero and so cannot reach the lateral sign factors at all

The two turbulent cases are the ones that demonstrate the controller claim.  The
DISCON library is used completely unchanged, and blade pitch, generator torque,
generator power, generator speed and rotor speed all come out identical between
the clockwise and mirrored runs.  They read a ``y``-reflected copy of the
turbulence box, for the reason given above.

For a clockwise rotor every mirror-related expression reduces to a multiplication
by ``+1``, so existing regression baselines reproduce bit-for-bit.
