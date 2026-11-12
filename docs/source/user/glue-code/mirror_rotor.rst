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

This page covers everything needed to **use** the flag.  How the mirror is
implemented inside the code, and how it is verified, is documented separately in
:numref:`glue-code-mirror-rotor-implementation`.

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

Nothing else changes.  The ElastoDyn, BeamDyn, AeroDyn and blade input files,
the airfoil polars and the controller are all supplied exactly as they would be
for the clockwise machine, and the turbine described by those files is mirrored
internally.

Airfoil tables in particular are used **verbatim** — there is no polar
transformation, and the airfoil coordinate files are not modified.

.. _glue-code-mirror-rotor-summary:

What to expect
--------------

A mirrored rotor is the exact mirror image of the clockwise machine described
by the input files:

- **The rotor turns the other way**, and the blades are numbered in the
  direction of rotation, so blade 2 of a mirrored three-bladed rotor is found
  120° *behind* blade 1 in the clockwise sense — the mirror image of its
  clockwise position.
- **Aerodynamic performance is identical.**  Power, thrust, angle of attack and
  the lift and drag coefficients match the clockwise machine; lateral forces
  and in-plane moments act in the mirrored directions.
- **The controller needs no change** — see :ref:`glue-code-mirror-rotor-control`.
- **The wake rotates the other way**, in OLAF and in FAST.Farm, and the VTK
  visualisation draws the blade sections as mirror images.
- **Everything outside the rotor behaves as described in its own input files.**
  The inflow, waves, support structure and mooring are not mirrored, so the
  turbine responds as the real machine would: a counter-clockwise rotor in the
  conditions you actually have.

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

A complete, measured classification of every output channel — which are
identical between a clockwise rotor and its mirror image, and which change
sign — is maintained with the implementation documentation
(:ref:`glue-code-mirror-rotor-signs`).

.. _glue-code-mirror-rotor-control:

Controllers
-----------

ServoDyn is not told that the rotor has been reversed.  Everything crossing its
boundary is presented in the **clockwise convention**, so an unmodified
controller — including a Bladed-style DLL such as ROSCO — sees exactly what it
would see on a clockwise machine and behaves identically.  Rotor speed, blade
pitch, generator and brake torque, shaft azimuth and the blade root moments are
all converted; the controller needs no mirrored copy and no new input.

Assumptions of even blade spacing *in the direction of rotation* also hold,
because the blades are numbered that way for either convention.

Yaw is the exception, and it is deliberate.  Yaw acts about the vertical axis in
the inertial frame, so it is **not** a rotor-convention quantity and is left
alone.  The yaw angle, the wind direction and hence the yaw error all stay
physically correct, which means a yaw controller still points the nacelle into
the real wind rather than into its mirror image.

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
   * - AeroAcoustics
     - ``CompAA = True`` in the AeroDyn input file

FAST.Farm is supported: a farm may mix clockwise and mirrored turbines, and a
mirrored turbine's wake swirl reverses.  The restriction that only one rotor
per OpenFAST instance is supported with FAST.Farm is unrelated to the mirror
and remains.  Setting ``MirrorRotor`` on a furling turbine raises a warning:
the tail-furl drag model does not depend on the rotation direction, but the
combination has not been verified.

.. _glue-code-mirror-rotor-further:

Implementation and verification
-------------------------------

Where each sign lives in the code, the BeamDyn blade-matrix transformation, the
free-wake and visualisation details, the measured per-channel sign map, and the
regression cases and guard scripts that verify all of it are documented in
:numref:`glue-code-mirror-rotor-implementation`.
