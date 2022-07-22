.. _ad_intro:

Introduction
============

AeroDyn is a time-domain wind turbine aerodynamics module that is coupled in the
OpenFAST multi-physics engineering tool to enable aero-elastic simulation of
horizontal-axis turbines.  AeroDyn can also be driven as a standalone code to
compute wind turbine aerodynamic response uncoupled from OpenFAST.  When coupled
to OpenFAST, AeroDyn can also be linearized as part of the linearization of the
full coupled solution (linearization is not available in standalone mode).
AeroDyn was originally developed for modeling wind turbine aerodynamics.
However, the module equally applies to the hydrodynamics of marine hydrokinetic
(MHK) turbines (the terms “wind turbine”, “tower”, “aerodynamics” etc.  in this
document imply “MHK turbine”, “MHK support structure”, “hydrodynamics” etc. for
MHK turbines).  Additional physics important for MHK turbines, not applicable to
wind turbines, computed by AeroDyn include a cavitation check.  This
documentation pertains version of AeroDyn in the OpenFAST github repository.
The AeroDyn version released of OpenFAST 1.0.0 is most closely related to
AeroDyn version 15 in the legacy version numbering.  AeroDyn version 15 was a
complete overhaul from earlier version of AeroDyn.  AeroDyn version 15 and newer
follows the requirements of the FAST modularization framework. 

AeroDyn calculates aerodynamic loads on both the blades and tower.
Aerodynamic calculations within AeroDyn are based on the principles of
actuator lines, where the three-dimensional (3D) flow around a body is
approximated by local two-dimensional (2D) flow at cross sections, and
the distributed pressure and shear stresses are approximated by lift
forces, drag forces, and pitching moments lumped at a node in a 2D cross
section. Analysis nodes are distributed along the length of each blade
and tower, the 2D forces and moment at each node are computed as
distributed loads per unit length, and the total 3D aerodynamic loads
are found by integrating the 2D distributed loads along the length. When
AeroDyn is coupled to OpenFAST, the blade and tower analysis node
discretization may be independent from the discretization of the nodes
in the structural modules. The actuator line approximations restrict the
validity of the model to slender structures and 3D behavior is either
neglected, captured through corrections inherent in the model (e.g.,
tip-loss, hub-loss, or skewed-wake corrections), or captured in the
input data (e.g., rotational augmentation corrections applied to airfoil
data).

AeroDyn assumes the turbine geometry consists of a one-, two-, or
three-bladed rotor atop a single tower. While the undeflected tower is
assumed to be straight and vertical, an undeflected blade may consider
out-of-plane curvature and in-plane sweep. For blades, the 2D cross
sections where the aerodynamic analysis take place may follow the
out-of-plane curvature, but in-plane sweep is assumed to be accomplished
by shearing, rather than rotation of the 2D cross section. Aerodynamic
imbalances are possible through the use of geometrical differences
between each blade.

When AeroDyn is coupled to OpenFAST, AeroDyn receives the instantaneous
(possibly displaced/deflected) structural position, orientation, and
velocities of analysis nodes in the tower, hub, and blades. As with
curvature and sweep, the 2D cross sections where the blade aerodynamic
analysis takes place will follow the out-of-plane deflection, but
in-plane deflection is assumed to be accomplished by shearing, rather
than rotation of the 2D cross section. AeroDyn also receives the local
freestream (undisturbed) fluid velocities at the tower and blade nodes.
(Fluid and structural calculations take place outside of the AeroDyn
module and are passed as inputs to AeroDyn by the driver code.) The
fluid and structural motions are provided at each coupling time step and
then AeroDyn computes the aerodynamic loads on the blade and tower nodes
and returns them back to OpenFAST as part of the aero-elastic calculation.
In standalone mode, the inputs to AeroDyn are prescribed by a simple
driver code, without aero-elastic coupling.

AeroDyn consists of four submodels: (1) rotor wake/induction, (2) blade
airfoil aerodynamics, (3) tower influence on the fluid local to the
blade nodes, and (4) tower drag. Nacelle, hub, and tail-vane fluid
influence and loading, aeroacoustics, and wake and array effects between
multiple turbines in a wind plant, are not yet available in AeroDyn v15
and newer.

For operating wind and MHK turbine rotors, AeroDyn calculates the
influence of the wake via induction factors based on the quasi-steady
Blade-Element/Momentum (BEM) theory, which requires an iterative
nonlinear solve (implemented via Brent’s method). By quasi-steady, it is
meant that the induction reacts instantaneously to loading changes. The
induction calculation, and resulting inflow velocities and angles, are
based on flow local to each analysis node of each blade, based on the
relative velocity between the fluid and structure (including the effects
of local inflow skew, shear, turbulence, tower flow disturbances, and
structural motion, depending on features enabled). The Glauert’s
empirical correction (with Buhl’s modification) replaces the linear
momentum balance at high axial induction factors. In the BEM solution,
Prandtl tip-loss, Prandtl hub-loss, and Pitt and Peters skewed-wake are
all 3D corrections that can optionally be applied. When the skewed-wake
correction is enabled, it is applied after the BEM iteration.
Additionally, the calculation of tangential induction (from the angular
momentum balance), the use of drag in the axial-induction calculation,
and the use of drag in the tangential-induction calculation are all
terms that can optionally be included in the BEM iteration (even when
drag is not used in the BEM iteration, drag is still used to calculate
the nodal loads once the induction has been found). The wake/induction
calculation can be bypassed altogether for the purposes of modeling
rotors that are parked or idling, in which case the inflow velocity and
angle are determined purely geometrically. During linearization analyses
with AeroDyn coupled to OpenFAST and BEM enabled, the wake can be assumed to
be frozen (i.e., the axial and tangential induces velocities, :math:`-V_x a` and :math:`V_y a'`, are
fixed at their operating-point values during linearization) or the
induction can be recalculated during linearization using BEM theory.
Dynamic wake that accounts for induction dynamics as a result of
transient conditions are not yet available in AeroDyn v15 and newer.

The blade airfoil aerodynamics can be steady or unsteady, except in the
case that a cavitation check is requested for MHK, in which case only
steady aerodynamics are supported. In the steady model, the supplied
static airfoil data — including the lift force, drag force, and optional
pitching moment and minimum pressure coefficients versus angle of attack
(AoA) — are used directly to calculate nodal loads. The
`AirfoilPrep <https://nwtc.nrel.gov/AirFoilPrep>`__ preprocessor can be
used to generate the needed static airfoil data based on uncorrected 2D
data (based, e.g., on airfoil tests in a wind tunnel or
`XFoil <http://web.mit.edu/drela/Public/web/xfoil/>`__), including
features to blend data between different airfoils, apply 3D rotational
augmentation, and extrapolate to high AoA. The unsteady airfoil
aerodynamic (UA) models account for flow hysteresis, including unsteady
attached flow, trailing-edge flow separation, dynamic stall, and flow
reattachment. The UA models can be considered as 2D dynamic corrections
to the static airfoil response as a result of time-varying inflow
velocities and angles. Three semi-empirical UA models are available: the
original theoretical developments of Beddoes-Leishman (B-L), extensions
to the B-L developed by González, and extensions to the B-L model
developed by Minnema/Pierce. **While all of the UA models are documented
in this manual, the original B-L model is not yet functional. Testing
has shown that the González and Minnema/Pierce models produce reasonable
hysteresis of the normal force, tangential force, and pitching-moment
coefficients if the UA model parameters are set appropriately for a
given airfoil, Reynolds number, and/or Mach number. However, the
results will differ a bit from earlier versions of AeroDyn, (which was
based on the Minnema/Pierce extensions to B-L) even if the default UA
model parameters are used, due to differences in the UA model logic
between the versions. We recommend that users run test cases with
uniform wind inflow and fixed yaw error (e.g., through the standalone
AeroDyn driver) to examine the accuracy of the normal force, tangential
force, and pitching-moment coefficient hysteresis and to adjust the UA
model parameters appropriately.** The airfoil-, Reynolds-, and
Mach-dependent parameters of the UA models may be derived from the
static airfoil data. These UA models are valid for small to moderate AoA
under normal rotor operation; the steady model is more appropriate under
parked or idling conditions. The static airfoil data is always used in
the BEM iteration; when UA is enabled, it is applied after the BEM
iteration and after the skewed-wake correction. The UA models are not
set up to support linearization, so, UA must be disabled during
linearization analyses with AeroDyn coupled to OpenFAST. The interpolation
of airfoil data based on Reynolds number or aerodynamic-control setting
(e.g., flaps) is not yet available in AeroDyn v15 and newer.

The influence of the tower on the fluid flow local to the blade is based
on a potential-flow and/or a tower-shadow model. The potential-flow
model uses the analytical potential-flow solution for flow around a
cylinder to model the tower dam effect on upwind rotors. In this model,
the freestream (undisturbed) flow at each blade node is disturbed based
on the location of the blade node relative to the tower and the tower
diameter, including lower velocities upstream and downstream of the
tower, higher velocities to the left and right of the tower, and
cross-stream flow. The Bak correction can optionally be included in the
potential-flow model, which augments the tower upstream disturbance and
improves the tower wake for downwind rotors based on the tower drag
coefficient. The tower shadow model can also be enabled to account for
the tower wake deficit on downwind rotors. This model includes an axial
flow deficit on the freestream fluid at each blade node dependent on the
location of the blade node relative to the tower and the tower diameter
and drag coefficient, based on the work of Powles. Both tower-influence
models are quasi-steady models, in that the disturbance is applied
directly to the freestream fluid at the blade nodes without dynamics,
and are applied within the BEM iteration.

The aerodynamic load on the tower is based directly on the tower
diameter and drag coefficient and the local relative fluid velocity
between the freestream (undisturbed) flow and structure at each tower
analysis node (including the effects of local shear, turbulence, and
structural motion, depending on features enabled). The tower drag load
calculation is quasi-steady and independent from the tower influence on
flow models.

The primary AeroDyn input file defines modeling options, environmental
conditions (except freestream flow), airfoils, tower nodal
discretization and properties, as well as output file specifications.
Airfoil data properties are read from dedicated inputs files (one for
each airfoil) and include coefficients of lift force, drag force, and
optional pitching moment and minimum pressure versus AoA, as well as UA
model parameters. (Minimum pressure coefficients versus AoA are also
included in the airfoil input files in case that a cavitation check is
requested.) Blade nodal discretization, geometry, twist, chord, and
airfoil identifier are likewise read from separate input files (one for
each blade).

:numref:`ad_input` describes the AeroDyn input files. 
:numref:`ad_output` discusses the
output files generated by AeroDyn; these include an echo file, summary
file, and the results file. 
:numref:`ad_modeling` provides modeling guidance when
using AeroDyn. 
Example input files are included in :numref:`ad_input_files`. A summary of
available output channels are found :numref:`ad_output_channels`. 
