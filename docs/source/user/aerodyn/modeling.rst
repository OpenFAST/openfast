.. _ad_modeling:

Modeling Considerations
=======================


AeroDyn was designed as an extremely flexible tool for modeling a
wide-range of aerodynamic conditions and turbine configurations. This
section provides some general guidance to help you construct models that
are compatible with AeroDyn.

Please refer to the theory of Section 7 for detailed information about
the implementation approach we have followed in AeroDyn.

Standalone AeroDyn Driver
-------------------------

The standalone AeroDyn driver code is very useful for computing turbine
aerodynamics independent of aero-elastic coupling. The standalone
AeroDyn driver code essentially replaces the functionality previously
available in the separate wind turbine rotor-performance tool WT\_Perf.
For example, the standalone AeroDyn driver code can be used to compute
the surfaces of power coefficient (C\ :sub:`P`), thrust coefficient
(C\ :sub:`T`), and/or torque coefficient (C\ :sub:`Q`) as a function of
tip-speed ratio (TSR) and blade-pitch angle for a given rotor. Moreover,
the standalone AeroDyn driver code is more powerful than WT\_Perf in
that the standalone AeroDyn driver can capture time-varying dynamics as
a result of nacelle-yaw error, shaft tilt, and/or wind shear.

Environmental Conditions
------------------------

For air, typical values for ``AirDens``, ``KinVisc``,
``SpdSound``, and ``Patm`` are around 1.225 kg/m\ :sup:`3`, 1.460E-5
m\ :sup:`2`/s, 340.3 m/s, and 101,325 Pa, respectively. For seawater,
typical values for ``AirDens``, ``KinVisc``, and ``Pvap`` are
around 1025 kg/m\ :sup:`3`, 1.004E-6 m\ :sup:`2`/s, and 2000 Pa,
respectively.

Temporal and Spatial Discretization
-----------------------------------

For accuracy and numerical stability, we recommend that ``DTAero`` be
set such that there are at least 200 azimuth steps per rotor revolution.
However, when AeroDyn is coupled to FAST, FAST may require time steps
much smaller than this rule of thumb. If UA is enabled while using very
small time steps, you may need to recompile AeroDyn in double precision
to avoid numerical problems in the UA routines.

For the blade and tower spatial discretization, using higher number of
analysis nodes will result in a more accurate solution at the expense of
longer computational time. When AeroDyn is coupled to FAST, the blade
and tower analysis node discretization may be independent from the
discretization of the nodes in the structural modules.

We recommend that ``NumBlNds`` be between 10 and 20 to balance
accuracy with computational expense for the rotor aerodynamic load
calculation. It may be beneficial to use a finer resolution of nodes
where large gradients are expected in the aerodynamic loads e.g. near
the blade tip. Aerodynamic imbalances are possible through the use of
geometrical differences between each blade.

When the tower potential-flow (``TwrPotent`` > 0), tower shadow
(``TwrShadow`` = TRUE), and/or the tower aerodynamic load
(``TwrAero`` = TRUE) models are enabled, we also recommend that
``NumTwrNds`` be between 10 and 20 to balance accuracy with
computational expense. Normally the local elevation of the tower node
above ground (or above MSL for offshore wind turbines or above the
seabed for MHK turbines) (``TwrElev``), must be entered in
monotonically increasing order from the lowest (tower-base) to the
highest (tower-top) elevation. However, when AeroDyn is coupled to FAST,
the tower-base node in AeroDyn cannot be set lower than the lowest point
where wind is specified in the InflowWind module. To avoid truncating
the lower section of the tower in AeroDyn, we recommend that the wind be
specified in InflowWind as low to the ground (or MSL for offshore wind
turbines or above the seabed for MHK turbines) as possible (this is a
particular issue for full-field wind file formats).

Model Options Under Operational and Parked/Idling Conditions
------------------------------------------------------------

To model an operational rotor, we recommend to include induction
(``WakeMod = 1``) and UA (``AFAeroMod = 2``). Normally, the Pitt and
Peters skewed-wake (``SkewMod = 2``), Prandtl tip-loss (``TipLoss
= TRUE``), Prandtl hub-loss (``HubLoss = TRUE``), and tangential
induction (``TanInd = TRUE``) models should all be enabled, but
``SkewMod = 2`` is invalid for very large yaw errors (much greater
than 45 degrees). The nonlinear solve in the BEM solution is in terms of the
inflow angle, but ``IndToler`` represents the tolerance of the
nondimensional residual equation, with no physical association possible;
we recommend setting ``IndToler`` to ``DEFAULT``.

*While all of the UA models are documented in this manual, the original
B-L model is not yet functional. Testing has shown that the González and
Minnema/Pierce models produce reasonable hysteresis of the normal force,
tangential force, and pitching-moment coefficients if the UA model
parameters are set appropriately for a given airfoil, Reynold’s number,
and/or Mach number. However, the results will differ a bit from earlier
versions of AeroDyn, (which was based on the Minnema/Pierce extensions
to B-L) even if the default UA model parameters are used, due to
differences in the UA model logic between the versions. We recommend
that users run test cases with uniform inflow and fixed yaw error (e.g.,
through the standalone AeroDyn driver) to examine the accuracy of the
normal force, tangential force, and pitching-moment coefficient
hysteresis and to adjust the UA model parameters appropriately.*

To model a parked or idling rotor, we recommend to disable induction
(``WakeMod = 0``) and UA (``AFAeroMod = 1``), in which case the
inflow velocity and angle are determined purely geometrically and the
airfoil data is determined statically.

The direct aerodynamic load on the tower often dominates the aerodynamic
load on the rotor for parked or idling conditions above the cut-out wind
speed, in which case we recommend that ``TwrAero = TRUE``. Otherwise,
``TwrAero = FALSE`` may be satisfactory.

We recommend to include the influence of the tower on the fluid local to
the blade for both operational and parked/idling rotors. We recommend
that ``TwrPotent`` > 0 for upwind rotors and that ``TwrPotent = 2``
or ``TwrShadow = TRUE`` for downwind rotors.

Linearization
-------------


When coupled to FAST, AeroDyn can be linearized as part of the
linearization of the full coupled solution. When induction is enabled
(``WakeMod = 1``), we recommend to base the linearized solution on the
frozen-wake assumption, by setting ``FrozenWake = TRUE``. The UA
models are not set up to support linearization, so, UA must be disabled
during linearization by setting ``AFAeroMod = 1``.
