.. _glue-code-mirror-rotor-sign-map:

Mirror-rotor sign map
=====================

.. warning::
   This page is **generated**.  Do not edit it by hand; regenerate it with

   .. code-block:: bash

      python3 reg_tests/otherTests/emit_sign_table.py --tol 0.005 --emit-map

Every entry below is **measured**, not asserted.  Each registered
clockwise/mirrored pair is run and every output channel is compared against its
counterpart, so this records behaviour that is observed rather than behaviour
expected from reading the source.  A channel seen in more than one pair must
agree across them; a disagreement is reported rather than silently resolved.

The grouped, prose form of the same information is in
:ref:`glue-code-mirror-rotor-verification`, which is the better place to start.
This page exists for looking a single channel up, and for anything that wants to
consume the map as data -- see ``mirror_rotor_sign_map.yaml`` beside this file.

**Read it correctly.**  This records how the two runs of a *symmetric* comparison
relate to each other, which is how the implementation is verified.  It is **not**
a claim that these channels change sign whenever the flag is set: in an ordinary
simulation the tower, support structure and inflow are not mirrored, so a
quantity such as ``TwrBsMxt`` is simply the response of an unchanged structure to
a counter-clockwise rotor.

Blade 1 lies on the mirror plane and blades 2 and 3 exchange, so a mirrored
blade 2 is compared against the clockwise blade 3.  Mooring channels are set
aside because which line pairs with which depends on the layout.

Measured at a relative tolerance of ``0.005`` with a noise floor of ``1e-08``, across 1203 channels.

* **identical** (686) -- v' = v
* **sign-flipped** (328) -- v' = -v
* **mirrored angle** (1) -- v' = -v, wrapped
* **below the noise floor** (55) -- indistinguishable from zero everywhere
* **set aside** (133) -- mooring, pairing is layout-specific

.. list-table::
   :header-rows: 1
   :widths: 30 25 45

   * - Channel
     - Behaviour
     - Measured in
   * - ``AB1N001Alpha``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N001Cd``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Cl``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Cm``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Cn``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Ct``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Cx``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Cy``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Fd``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Fl``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Fn``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Ft``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Fx``
     - below the noise floor
     - AeroDyn nodal outputs
   * - ``AB1N001Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N001Gam``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Mm``
     - below the noise floor
     - AeroDyn nodal outputs
   * - ``AB1N001Phi``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001STVx``
     - below the noise floor
     - OLAF free wake
   * - ``AB1N001STVy``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001STVz``
     - below the noise floor
     - OLAF free wake
   * - ``AB1N001Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N001TnInd``
     - below the noise floor
     - AeroDyn nodal outputs
   * - ``AB1N001VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N001Vindx``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001Vindy``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001VRel``
     - identical
     - OLAF free wake
   * - ``AB1N001Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N001VUndx``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N001VUndy``
     - sign-flipped
     - OLAF free wake
   * - ``AB1N001VUndz``
     - identical
     - OLAF free wake
   * - ``AB1N002Alpha``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N002Cd``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Cl``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Cm``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Cn``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Ct``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Cx``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Cy``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Fd``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Fl``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Fn``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Ft``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N002Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N002Gam``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N002Phi``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002STVx``
     - below the noise floor
     - OLAF free wake
   * - ``AB1N002STVy``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002STVz``
     - below the noise floor
     - OLAF free wake
   * - ``AB1N002Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N002TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N002VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N002Vindx``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002Vindy``
     - sign-flipped
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002VRel``
     - identical
     - OLAF free wake
   * - ``AB1N002Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N002VUndx``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``AB1N002VUndy``
     - sign-flipped
     - OLAF free wake
   * - ``AB1N002VUndz``
     - identical
     - OLAF free wake
   * - ``AB1N003Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N003Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N003Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N003Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N003Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N003Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N003Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N003Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N003Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N003VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N004Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N004Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N004Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N004Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N004Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N004Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N004Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N004Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N004VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N005Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N005Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N005Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N005Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N005Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N005Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N005Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N005Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N005VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N006Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N006Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N006Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N006Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N006Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N006Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N006Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N006Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N006VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N007Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N007Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N007Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N007Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N007Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N007Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N007Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N007Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N007VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N008Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N008Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N008Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N008Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N008Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N008Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N008Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N008Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N008VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N009Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N009Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N009Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N009Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N009Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N009Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N009Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N009Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N009VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N010Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N010Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N010Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N010Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N010Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N010Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N010Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N010Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N010VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N011Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N011Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N011Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N011Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N011Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N011Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N011Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N011Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N011VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N012Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N012Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N012Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N012Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N012Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N012Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N012Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N012Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N012VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N013Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N013Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N013Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N013Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N013Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N013Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N013Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N013Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N013VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N014Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N014Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N014Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N014Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N014Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N014Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N014Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N014Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N014VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N015Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N015Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N015Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N015Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N015Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N015Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N015Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N015Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N015VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N016Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N016Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N016Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N016Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N016Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N016Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N016Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N016Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N016VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N017Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N017Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N017Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N017Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N017Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N017Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N017Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N017Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N017VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N018Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N018Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N018Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N018Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N018Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N018Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N018Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N018Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N018VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N019Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N019Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N019Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N019Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N019Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N019Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N019Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N019Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N019VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N020Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N020Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N020Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N020Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N020Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N020Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N020Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N020Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N020VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N021Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N021Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N021Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N021Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N021Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N021Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N021Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N021Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N021VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N022Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N022Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N022Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N022Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N022Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N022Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N022Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N022Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N022VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N023Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N023Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N023Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N023Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N023Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N023Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N023Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N023Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N023VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N024Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N024Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N024Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N024Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N024Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N024Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N024Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N024Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N024VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N025Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N025Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N025Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N025Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N025Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N025Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N025Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N025Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N025VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N026Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N026Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N026Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N026Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N026Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N026Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N026Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N026Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N026VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N027Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N027Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N027Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N027Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N027Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N027Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N027Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N027Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N027VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N028Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N028Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N028Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N028Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N028Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N028Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N028Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N028Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N028VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N029Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N029Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N029Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N029Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N029Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N029Phi``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N029Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029TnInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029Vindy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N029Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N029VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Alpha``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030AxInd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Cd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Cl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Cm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N030Cn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Ct``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N030Cx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Cy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N030Fd``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Fl``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Fn``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Ft``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N030Fx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Fy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N030Gam``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Mm``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N030Phi``
     - below the noise floor
     - AeroDyn nodal outputs
   * - ``AB1N030STVy``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``AB1N030Theta``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030TnInd``
     - below the noise floor
     - AeroDyn nodal outputs
   * - ``AB1N030VDisx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Vindx``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030Vindy``
     - below the noise floor
     - AeroDyn nodal outputs
   * - ``AB1N030Vrel``
     - identical
     - AeroDyn nodal outputs
   * - ``AB1N030VUndx``
     - identical
     - AeroDyn nodal outputs
   * - ``ADCp``
     - identical
     - AeroDisk, yawed
   * - ``ADCq``
     - identical
     - AeroDisk, yawed
   * - ``ADCt``
     - identical
     - AeroDisk, yawed
   * - ``ADFx``
     - identical
     - AeroDisk, yawed
   * - ``ADFxi``
     - identical
     - AeroDisk, yawed
   * - ``ADFy``
     - identical
     - AeroDisk, yawed
   * - ``ADFyi``
     - sign-flipped
     - AeroDisk, yawed
   * - ``ADFz``
     - identical
     - AeroDisk, yawed
   * - ``ADFzi``
     - identical
     - AeroDisk, yawed
   * - ``ADMx``
     - sign-flipped
     - AeroDisk, yawed
   * - ``ADMxi``
     - sign-flipped
     - AeroDisk, yawed
   * - ``ADMy``
     - sign-flipped
     - AeroDisk, yawed
   * - ``ADMyi``
     - identical
     - AeroDisk, yawed
   * - ``ADMz``
     - sign-flipped
     - AeroDisk, yawed
   * - ``ADMzi``
     - sign-flipped
     - AeroDisk, yawed
   * - ``ADPitch``
     - below the noise floor
     - AeroDisk, yawed
   * - ``ADPower``
     - identical
     - AeroDisk, yawed
   * - ``ADSkew``
     - identical
     - AeroDisk, yawed
   * - ``ADSpeed``
     - sign-flipped
     - AeroDisk, yawed
   * - ``ADSTVx``
     - below the noise floor
     - AeroDisk, yawed
   * - ``ADSTVxi``
     - below the noise floor
     - AeroDisk, yawed
   * - ``ADSTVy``
     - below the noise floor
     - AeroDisk, yawed
   * - ``ADSTVyi``
     - below the noise floor
     - AeroDisk, yawed
   * - ``ADSTVz``
     - below the noise floor
     - AeroDisk, yawed
   * - ``ADSTVzi``
     - below the noise floor
     - AeroDisk, yawed
   * - ``ADTSR``
     - identical
     - AeroDisk, yawed
   * - ``ADVRel``
     - identical
     - AeroDisk, yawed
   * - ``ADVWindx``
     - identical
     - AeroDisk, yawed
   * - ``ADVWindxi``
     - identical
     - AeroDisk, yawed
   * - ``ADVWindy``
     - identical
     - AeroDisk, yawed
   * - ``ADVWindyi``
     - below the noise floor
     - AeroDisk, yawed
   * - ``ADVWindz``
     - below the noise floor
     - AeroDisk, yawed
   * - ``ADVWindzi``
     - below the noise floor
     - AeroDisk, yawed
   * - ``ADYawErr``
     - below the noise floor
     - AeroDisk, yawed
   * - ``ANCHTEN1``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``ANCHTEN2``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``ANCHTEN3``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``ANCHTEN4``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``Azimuth``
     - identical
     - AeroDisk, yawed, AeroDyn nodal outputs, BeamDyn blades, ElastoDyn + AeroDyn, OLAF free wake, Twin-rotor semisubmersible
   * - ``B1AeroPwr``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1FldFz``
     - identical
     - MHK buoyancy
   * - ``B1N1Alpha``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N1Cd``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N1Cl``
     - below the noise floor
     - ElastoDyn + AeroDyn
   * - ``B1N1Cm``
     - below the noise floor
     - ElastoDyn + AeroDyn
   * - ``B1N1Cn``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N1Ct``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N1Cx``
     - below the noise floor
     - ElastoDyn + AeroDyn
   * - ``B1N1Cy``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N1Fd``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N1Fl``
     - below the noise floor
     - ElastoDyn + AeroDyn
   * - ``B1N1Fn``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N1Ft``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N1Fx``
     - below the noise floor
     - ElastoDyn + AeroDyn
   * - ``B1N1Fy``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N1Theta``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N1VIndx``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N1VIndy``
     - below the noise floor
     - ElastoDyn + AeroDyn
   * - ``B1N2Alpha``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N2Cd``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N2Cl``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N2Cm``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N2Cn``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N2Ct``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N2Cx``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N2Cy``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N2Fd``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N2Fl``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N2Fn``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N2Ft``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N2Fx``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N2Fy``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N2Mbn``
     - sign-flipped
     - MHK buoyancy
   * - ``B1N2SgCav``
     - identical
     - MHK buoyancy
   * - ``B1N2Theta``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N2VIndx``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N2VIndy``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N3Alpha``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N3Cd``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N3Cl``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N3Clrnc``
     - identical
     - Twin-rotor semisubmersible
   * - ``B1N3Cm``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N3Cn``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N3Ct``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N3Cx``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N3Cy``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N3Fd``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N3Fl``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N3Fn``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N3Ft``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N3Fx``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N3Fy``
     - sign-flipped
     - ElastoDyn + AeroDyn
   * - ``B1N3SigCr``
     - identical
     - MHK buoyancy
   * - ``B1N3Theta``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``B1N3VIndx``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B1N3VIndy``
     - below the noise floor
     - ElastoDyn + AeroDyn
   * - ``B1N6Mbs``
     - sign-flipped
     - MHK buoyancy
   * - ``B1N7Fbt``
     - sign-flipped
     - MHK buoyancy
   * - ``B1RootFxr``
     - identical
     - BeamDyn blades
   * - ``B1RootFyr``
     - sign-flipped
     - BeamDyn blades
   * - ``B1RootFzr``
     - identical
     - BeamDyn blades
   * - ``B1RootMxr``
     - sign-flipped
     - BeamDyn blades
   * - ``B1RootMyr``
     - identical
     - BeamDyn blades
   * - ``B1RootMzr``
     - sign-flipped
     - BeamDyn blades
   * - ``B1TipRDxr``
     - sign-flipped
     - BeamDyn blades
   * - ``B1TipRDyr``
     - identical
     - BeamDyn blades
   * - ``B1TipRDzr``
     - sign-flipped
     - BeamDyn blades
   * - ``B1TipTDxr``
     - identical
     - BeamDyn blades
   * - ``B1TipTDyr``
     - sign-flipped
     - BeamDyn blades
   * - ``B1TipTDzr``
     - identical
     - BeamDyn blades
   * - ``B2AeroPwr``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B2FldMx``
     - sign-flipped
     - MHK buoyancy
   * - ``B2N3Clrnc``
     - identical
     - Twin-rotor semisubmersible
   * - ``B2N3Mbt``
     - identical
     - MHK buoyancy
   * - ``B2N4Fbn``
     - identical
     - MHK buoyancy
   * - ``B2N5SigCr``
     - identical
     - MHK buoyancy
   * - ``B2N6SgCav``
     - identical
     - MHK buoyancy
   * - ``B2N8Fbs``
     - identical
     - MHK buoyancy
   * - ``B2RootFxr``
     - identical
     - BeamDyn blades
   * - ``B2RootFyr``
     - sign-flipped
     - BeamDyn blades
   * - ``B2RootFzr``
     - identical
     - BeamDyn blades
   * - ``B2RootMxr``
     - sign-flipped
     - BeamDyn blades
   * - ``B2RootMyr``
     - identical
     - BeamDyn blades
   * - ``B2RootMzr``
     - sign-flipped
     - BeamDyn blades
   * - ``B2TipRDxr``
     - sign-flipped
     - BeamDyn blades
   * - ``B2TipRDyr``
     - identical
     - BeamDyn blades
   * - ``B2TipRDzr``
     - sign-flipped
     - BeamDyn blades
   * - ``B2TipTDxr``
     - identical
     - BeamDyn blades
   * - ``B2TipTDyr``
     - sign-flipped
     - BeamDyn blades
   * - ``B2TipTDzr``
     - identical
     - BeamDyn blades
   * - ``B3AeroPwr``
     - identical
     - ElastoDyn + AeroDyn
   * - ``B3N3Clrnc``
     - identical
     - Twin-rotor semisubmersible
   * - ``B3RootFxr``
     - identical
     - BeamDyn blades
   * - ``B3RootFyr``
     - sign-flipped
     - BeamDyn blades
   * - ``B3RootFzr``
     - identical
     - BeamDyn blades
   * - ``B3RootMxr``
     - sign-flipped
     - BeamDyn blades
   * - ``B3RootMyr``
     - identical
     - BeamDyn blades
   * - ``B3RootMzr``
     - sign-flipped
     - BeamDyn blades
   * - ``B3TipRDxr``
     - sign-flipped
     - BeamDyn blades
   * - ``B3TipRDyr``
     - identical
     - BeamDyn blades
   * - ``B3TipRDzr``
     - sign-flipped
     - BeamDyn blades
   * - ``B3TipTDxr``
     - identical
     - BeamDyn blades
   * - ``B3TipTDyr``
     - sign-flipped
     - BeamDyn blades
   * - ``B3TipTDzr``
     - identical
     - BeamDyn blades
   * - ``BldPitch1``
     - below the noise floor
     - AeroDyn nodal outputs, OLAF free wake, Twin-rotor semisubmersible
   * - ``BldPitch2``
     - below the noise floor
     - AeroDyn nodal outputs
   * - ``BldPitch3``
     - below the noise floor
     - AeroDyn nodal outputs
   * - ``BlPitch1``
     - below the noise floor
     - AeroDisk, yawed
   * - ``BlPitch2``
     - below the noise floor
     - AeroDisk, yawed
   * - ``BlPitch3``
     - below the noise floor
     - AeroDisk, yawed
   * - ``Case``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``CON10PX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON10PY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON10PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON11PX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON11PY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON11PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON1FX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON1FY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON1FZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON2FX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON2FY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON2FZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON3FX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON3FY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON3FZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON4FX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON4FY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON4FZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON5FX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON5FY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON5FZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON5PX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON5PY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON5PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON6FX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON6FY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON6FZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON6PX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON6PY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON6PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON7FX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON7FY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON7FZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON7PX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON7PY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON7PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON8FX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON8FY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON8FZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON8PX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON8PY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON8PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON9PX``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON9PY``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``CON9PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``FAIRTEN1``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``FAIRTEN2``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``FAIRTEN3``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``FAIRTEN4``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``GenAcc``
     - below the noise floor
     - AeroDisk, yawed
   * - ``GenAccel``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``GenPwr``
     - identical
     - Twin-rotor semisubmersible
   * - ``GenSpeed``
     - identical
     - AeroDisk, yawed, BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``GenTq``
     - identical
     - Twin-rotor semisubmersible
   * - ``HbFbx``
     - identical
     - MHK buoyancy
   * - ``HbFby``
     - sign-flipped
     - MHK buoyancy
   * - ``HbFbz``
     - identical
     - MHK buoyancy
   * - ``HbMbx``
     - sign-flipped
     - MHK buoyancy
   * - ``HbMby``
     - identical
     - MHK buoyancy
   * - ``HbMbz``
     - sign-flipped
     - MHK buoyancy
   * - ``HSShftPwr``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``HSShftTq``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``HWindSpeedX``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``HWindSpeedY``
     - below the noise floor
     - AeroDyn nodal outputs, OLAF free wake
   * - ``HWindSpeedZ``
     - below the noise floor
     - AeroDyn nodal outputs, OLAF free wake
   * - ``HydroFxi``
     - identical
     - MHK buoyancy
   * - ``HydroFyi``
     - sign-flipped
     - MHK buoyancy
   * - ``HydroFzi``
     - identical
     - MHK buoyancy
   * - ``IPDefl1``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``IPDefl2``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``IPDefl3``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``L1N10PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N11PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N12PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N13PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N14PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N15PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N16PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N17PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N18PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N19PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N1PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N20PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N21PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N22PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N23PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N24PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N25PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N26PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N27PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N28PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N29PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N2PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N30PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N31PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N32PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N33PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N34PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N35PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N36PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N37PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N38PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N39PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N3PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N40PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N4PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N5PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N6PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N7PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N8PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L1N9PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N10PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N11PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N12PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N13PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N14PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N15PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N16PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N17PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N18PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N19PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N1PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N20PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N21PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N22PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N23PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N24PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N25PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N26PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N27PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N28PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N29PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N2PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N30PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N31PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N32PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N33PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N34PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N35PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N36PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N37PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N38PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N39PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N3PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N40PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N4PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N5PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N6PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N7PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N8PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``L3N9PZ``
     - set aside, layout-specific
     - MHK buoyancy
   * - ``LSSGagMxa``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``LSSGagMya``
     - identical
     - Twin-rotor semisubmersible
   * - ``LSSGagMza``
     - sign-flipped
     - Twin-rotor semisubmersible
   * - ``LSShftFxa``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``LSShftFya``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``LSShftFza``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``LSShftMxa``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``LSShftTq``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``LSSTipAxa``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``LSSTipMya``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``LSSTipMza``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``LSSTipPxa``
     - mirrored angle
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``LSSTipVxa``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``NcFbx``
     - below the noise floor
     - MHK buoyancy
   * - ``NcFby``
     - below the noise floor
     - MHK buoyancy
   * - ``NcFbz``
     - below the noise floor
     - MHK buoyancy
   * - ``NcIMUTAxs``
     - identical
     - Twin-rotor semisubmersible
   * - ``NcIMUTAys``
     - sign-flipped
     - Twin-rotor semisubmersible
   * - ``NcIMUTAzs``
     - identical
     - Twin-rotor semisubmersible
   * - ``NcMbx``
     - below the noise floor
     - MHK buoyancy
   * - ``NcMby``
     - below the noise floor
     - MHK buoyancy
   * - ``NcMbz``
     - below the noise floor
     - MHK buoyancy
   * - ``OoPDefl1``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn, MHK buoyancy, Twin-rotor semisubmersible
   * - ``OoPDefl2``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``OoPDefl3``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``PtfmHeave``
     - identical
     - AeroDyn nodal outputs, MHK buoyancy, OLAF free wake, Twin-rotor semisubmersible
   * - ``PtfmPitch``
     - identical
     - AeroDyn nodal outputs, MHK buoyancy, OLAF free wake, Twin-rotor semisubmersible
   * - ``PtfmRoll``
     - sign-flipped
     - AeroDyn nodal outputs, MHK buoyancy, OLAF free wake, Twin-rotor semisubmersible
   * - ``PtfmSurge``
     - identical
     - AeroDyn nodal outputs, MHK buoyancy, OLAF free wake, Twin-rotor semisubmersible
   * - ``PtfmSway``
     - sign-flipped
     - AeroDyn nodal outputs, MHK buoyancy, OLAF free wake, Twin-rotor semisubmersible
   * - ``PtfmYaw``
     - sign-flipped
     - AeroDyn nodal outputs, MHK buoyancy, OLAF free wake, Twin-rotor semisubmersible
   * - ``RootFxb1``
     - identical
     - Twin-rotor semisubmersible
   * - ``RootFxc1``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RootFyb1``
     - sign-flipped
     - Twin-rotor semisubmersible
   * - ``RootFyc1``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RootFzb1``
     - identical
     - Twin-rotor semisubmersible
   * - ``RootFzc1``
     - identical
     - Twin-rotor semisubmersible
   * - ``RootMxb1``
     - sign-flipped
     - Twin-rotor semisubmersible
   * - ``RootMxc1``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RootMxc2``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``RootMxc3``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``RootMyb1``
     - identical
     - Twin-rotor semisubmersible
   * - ``RootMyc1``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RootMyc2``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``RootMyc3``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``RootMzb1``
     - sign-flipped
     - Twin-rotor semisubmersible
   * - ``RootMzc1``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RotAcc``
     - below the noise floor
     - AeroDisk, yawed
   * - ``RotAccel``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``RotPwr``
     - identical
     - AeroDisk, yawed, BeamDyn blades, ElastoDyn + AeroDyn
   * - ``RotSpeed``
     - identical
     - AeroDisk, yawed, AeroDyn nodal outputs, BeamDyn blades, ElastoDyn + AeroDyn, OLAF free wake, Twin-rotor semisubmersible
   * - ``RotThrust``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``RotTorq``
     - identical
     - AeroDisk, yawed, BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RtAeroCp``
     - identical
     - AeroDyn nodal outputs, ElastoDyn + AeroDyn
   * - ``RtAeroCq``
     - sign-flipped
     - AeroDyn nodal outputs, ElastoDyn + AeroDyn
   * - ``RtAeroCt``
     - identical
     - AeroDyn nodal outputs, ElastoDyn + AeroDyn
   * - ``RtAeroFxh``
     - identical
     - AeroDyn nodal outputs, BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RtAeroFyh``
     - sign-flipped
     - AeroDyn nodal outputs, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RtAeroFzh``
     - identical
     - AeroDyn nodal outputs, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RtAeroMxh``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RtAeroMyh``
     - identical
     - AeroDyn nodal outputs, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RtAeroMzh``
     - sign-flipped
     - ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``RtAeroPwr``
     - identical
     - AeroDyn nodal outputs, BeamDyn blades, ElastoDyn + AeroDyn
   * - ``RtArea``
     - identical
     - AeroDyn nodal outputs, Twin-rotor semisubmersible
   * - ``RtFldFxh``
     - identical
     - MHK buoyancy
   * - ``RtFldFyh``
     - sign-flipped
     - MHK buoyancy
   * - ``RtFldFzg``
     - identical
     - MHK buoyancy
   * - ``RtFldMxh``
     - sign-flipped
     - MHK buoyancy
   * - ``RtFldMyg``
     - identical
     - MHK buoyancy
   * - ``RtFldMzh``
     - sign-flipped
     - MHK buoyancy
   * - ``RtSkew``
     - identical
     - AeroDyn nodal outputs, BeamDyn blades, ElastoDyn + AeroDyn
   * - ``RtSpeed``
     - identical
     - AeroDyn nodal outputs, ElastoDyn + AeroDyn
   * - ``RtTSR``
     - identical
     - AeroDyn nodal outputs, BeamDyn blades, ElastoDyn + AeroDyn
   * - ``RtVAvgxh``
     - identical
     - AeroDyn nodal outputs, Twin-rotor semisubmersible
   * - ``RtVAvgyh``
     - sign-flipped
     - AeroDyn nodal outputs
   * - ``RtVAvgzh``
     - identical
     - AeroDyn nodal outputs
   * - ``ShearExp``
     - identical
     - AeroDyn nodal outputs, OLAF free wake
   * - ``Time``
     - identical
     - AeroDisk, yawed, AeroDyn nodal outputs, BeamDyn blades, ElastoDyn + AeroDyn, MHK buoyancy, OLAF free wake
   * - ``TipDxc1``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``TipDyc1``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn
   * - ``TTDspFA``
     - identical
     - Twin-rotor semisubmersible
   * - ``TTDspSS``
     - sign-flipped
     - Twin-rotor semisubmersible
   * - ``TTDspTwst``
     - below the noise floor
     - Twin-rotor semisubmersible
   * - ``TwN1Fbx``
     - below the noise floor
     - MHK buoyancy
   * - ``TwN1Mbx``
     - below the noise floor
     - MHK buoyancy
   * - ``TwN2Mby``
     - below the noise floor
     - MHK buoyancy
   * - ``TwN3Fby``
     - below the noise floor
     - MHK buoyancy
   * - ``TwN3Mbz``
     - below the noise floor
     - MHK buoyancy
   * - ``TwN4Fbz``
     - below the noise floor
     - MHK buoyancy
   * - ``TwrBsFxt``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``TwrBsFyt``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``TwrBsFzt``
     - identical
     - Twin-rotor semisubmersible
   * - ``TwrBsMxt``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``TwrBsMyt``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``TwrBsMzt``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``TwrTpTDxi``
     - identical
     - MHK buoyancy
   * - ``TwrTpTDyi``
     - sign-flipped
     - MHK buoyancy
   * - ``TwrTpTDzi``
     - identical
     - MHK buoyancy
   * - ``TwstDefl1``
     - below the noise floor
     - Twin-rotor semisubmersible
   * - ``Wave1Elev``
     - identical
     - MHK buoyancy
   * - ``Wind1VelX``
     - identical
     - AeroDisk, yawed, AeroDyn nodal outputs, BeamDyn blades, ElastoDyn + AeroDyn
   * - ``Wind1VelY``
     - below the noise floor
     - AeroDisk, yawed, AeroDyn nodal outputs, BeamDyn blades, ElastoDyn + AeroDyn
   * - ``Wind1VelZ``
     - below the noise floor
     - AeroDisk, yawed, AeroDyn nodal outputs, BeamDyn blades, ElastoDyn + AeroDyn
   * - ``Yaw``
     - sign-flipped
     - AeroDisk, yawed, AeroDyn nodal outputs, OLAF free wake
   * - ``YawBrFxp``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``YawBrFyp``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``YawBrFzp``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``YawBrMxp``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``YawBrMyp``
     - identical
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``YawBrMzp``
     - sign-flipped
     - BeamDyn blades, ElastoDyn + AeroDyn, Twin-rotor semisubmersible
   * - ``YawRate``
     - below the noise floor
     - AeroDisk, yawed
