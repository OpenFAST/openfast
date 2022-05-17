HydroDyn User Guide and Theory Manual
=====================================

.. toctree::
   :maxdepth: 1

   getting_started.rst
   input_files.rst
   output_files.rst
   modeling_considerations.rst
   future_work.rst
   zrefs.rst
   appendix.rst

HydroDyn is a time-domain hydrodynamics module that has been coupled
into the OpenFAST wind turbine
multi-physics engineering tool to enable aero-hydro-servo-elastic
simulation of offshore wind turbines. HydroDyn is applicable to both
fixed-bottom and floating offshore substructures. The current release of
HydroDyn integrates with OpenFAST through the FAST modularization framework.
HydroDyn can also be driven as a standalone code to compute hydrodynamic
loading uncoupled from OpenFAST.

In addition to this documentation, the following materials including presentation slides,
development plans, and publications are made available for reference.
Note that some of these may be outdated and pertain to older versions
of HydroDyn.

- :download:`Computation of Wave Loads Under Multidirectional Sea States for Floating Offshore Wind Turbines <https://www.nrel.gov/docs/fy14osti/61161.pdf>`
- :download:`Effects of Second-Order Hydrodynamic Forces on Floating Offshore Wind Turbines <https://www.nrel.gov/docs/fy14osti/60966.pdf>`
- :download:`State-Space Realization of the Wave-Radiation Force within FAST <https://www.nrel.gov/docs/fy13osti/58099.pdf>`
- :download:`Dynamics of Offshore Floating Wind Turbines—Model Development and Verification <https://dx.doi.org/10.1002/we.347>`
- :download:`Dynamics Modeling and Loads Analysis of an Offshore Floating Wind Turbine <https://www.nrel.gov/docs/fy08osti/41958.pdf>`
- :download:`Draft Implementation Plan - Changes in HydroDyn to Support Time-Varying Buoyancy Loads on Morison Members <../../../OtherSupporting/HydroDyn/HydroDyn_Plan_TCF_Morison.docx>`
- :download:`Implementation Plan - Modifications to State-Space Modules in HydroDyn to Support Multiple WAMIT Bodies <../../../OtherSupporting/HydroDyn/HydroDyn_Plan_TCF_NBodyStateSpace.docx>`
- :download:`Implementation Plan (Revised) - Changes in HydroDyn to Support Multiple WAMIT Bodies <../../../OtherSupporting/HydroDyn/HydroDyn_Plan_TCF_NBody.docx>`
- :download:`Implementation Plan - 2nd-order Forces Within HydroDyn <../../../OtherSupporting/HydroDyn/HydroDyn_2ndOrderForces_Plan.pdf>`
- :download:`Implementation Plan - 2nd-order Wave Kinematics Within HydroDyn <../../../OtherSupporting/HydroDyn/WAVE2_document.pdf>`
- :download:`Plan for Adding Wave Stretching to HydroDyn <../../../OtherSupporting/HydroDyn/HydroDyn_WaveStretching_Plan.docx>`
- :download:`Breaking Wave Modeling Approach for FAST <../../../OtherSupporting/HydroDyn/Breaking_Wave_Modeling_Approach_for_FAST.docx>`


HydroDyn allows for multiple approaches for calculating the hydrodynamic
loads on a structure:

* Potential-flow theory solution
* Strip-theory solution
* Hybrid combination of the tower

Waves generated internally
within HydroDyn can be regular (periodic) or irregular (stochastic) and
long-crested (unidirectional) or short-crested (wave energy spread
across a range of directions). Wave elevations or full wave kinematics
can also be generated externally and used within HydroDyn. Internally,
HydroDyn generates waves analytically for finite depth using first-order
(linear Airy) or first plus second-order wave theory :cite:`SharmaDean:1981`
with the option to include directional spreading, but wave
kinematics are only computed in the domain between the flat seabed and
still-water level (SWL) and no wave stretching or higher order wave
theories are included. The second-order hydrodynamic implementations
include time-domain calculations of difference- (mean- and slow-drift-)
and sum-frequency terms. To minimize computational expense, Fast Fourier
Transforms (FFTs) are applied in the summation of all wave frequency
components.

The potential-flow solution is applicable to substructures or members of
substructures that are large relative to a typical wavelength. The
potential-flow solution involves either frequency-to-time-domain
transforms or fluid-impulse theory (FIT). In the former, potential-flow
hydrodynamic loads include linear hydrostatic restoring, the added mass
and damping contributions from linear wave radiation (including
free-surface memory effects), and the incident-wave excitation from
first- and second-order diffraction (Froude-Kriloff and scattering). The
hydrodynamic coefficients (first and second order) required for the
potential-flow solution are frequency dependent and must be supplied by
a separate frequency-domain panel code (e.g., WAMIT) from a
pre-computation step. The radiation memory effect can be calculated
either through direct time-domain convolution or through a linear
state-space approach, with a state-space model derived through the
`SS_Fitting <https://www.nrel.gov/wind/nwtc/ss-fitting.html>`_
preprocessor. The second-order terms can be derived from the full
difference- and sum-frequency quadratic transfer functions (QTFs) or the
difference-frequency terms can be estimated via Standing et al.’s :cite:`Standing:1987`
extension to Newman’s approximation, based only on first-order
coefficients. The use of FIT is not yet documented in this manual.

The strip-theory solution may be preferable for substructures or members
of substructures that are small in diameter relative to a typical
wavelength. Strip-theory hydrodynamic loads can be applied across
multiple interconnected members, each with possible incline and taper,
and are derived directly from the undisturbed wave and current
kinematics at the undisplaced position of the substructure. The
strip-theory loads include the relative form of Morison’s equation for
the distributed fluid-inertia, added-mass, and viscous-drag components.
Additional distributed load components include axial loads from tapered
members and static buoyancy loads. Hydrodynamic loads are also applied
as lumped loads on member endpoints (called joints). It is also possible
to include flooding or ballasting of members, and the effects of marine
growth. The hydrodynamic coefficients required for this solution come
through user-specified dynamic-pressure, added-mass, and viscous-drag
coefficients.

For some substructures and sea conditions, the hydrodynamic loads from a
potential-flow theory should be augmented with the loads brought about
by flow separation.  For this, the viscous-drag component of the
strip-theory solution may be included with the potential-flow theory
solution.  Another option available is to supply a global damping matrix
(linear or quadratic) to the system to represent this effect.

When HydroDyn is coupled to OpenFAST, HydroDyn receives the position,
orientation, velocities, and accelerations of the (rigid or flexible)
substructure at each coupling time step and then computes the
hydrodynamic loads and returns them back to OpenFAST. At this time,
OpenFAST’s ElastoDyn structural-dynamics module assumes for a floating platform
that the substructure (floating platform) is a six degree-of-freedom
(DOF) rigid body. For fixed-bottom offshore wind turbines, OpenFAST’s SubDyn
module allows for structural flexibility of multi-member substructures
and the coupling to HydroDyn includes hydro-elastic effects.

The primary HydroDyn input file defines the substructure geometry,
hydrodynamic coefficients, incident wave kinematics and current,
potential-flow solution options, flooding/ballasting and marine growth,
and auxiliary parameters. The geometry of strip-theory members is
defined by joint coordinates of the undisplaced substructure in the
global reference system, with the origin at the intersection of the
undeflected tower centerline with mean sea level (MSL). A member
connects two joints; multiple members can use a common joint. The
hydrodynamic loads are computed at nodes, which are the resultant of
member refinement into multiple (**MDivSize** input) elements (nodes are
located at the ends of each element), and they are calculated by the
module. Member properties include outer diameter, thickness, and
dynamic-pressure, added-mass and viscous-drag coefficients. Member
properties are specified at the joints; if properties change from one
joint to the other, they will be linearly interpolated for the inner
nodes.

See :ref:`hd-getting-started` for details on how to download or
compile the HydroDyn and OpenFAST software executables, as well as
instructions for running HydroDyn standalone and coupled to OpenFAST.
