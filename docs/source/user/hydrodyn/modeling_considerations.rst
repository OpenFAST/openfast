.. _hd-modeling-considerations:

Modeling Considerations
=======================
HydroDyn was designed as an extremely flexible tool for modeling a
wide-range of hydrodynamic conditions and substructures. This section
provides some general guidance to help you construct models that are
compatible with HydroDyn.

.. TODO The Theory section does not yet exist.
.. Please refer to the theory of Section 7 for detailed information about
.. HydroDyn’s coordinate systems, and the implementation approach we have
.. followed in HydroDyn.

Waves
~~~~~
Waves generated internally within HydroDyn can be regular (periodic) or
irregular (stochastic) and long-crested (unidirectional) or
short-crested (with wave energy spread across a range of directions).
Internally, HydroDyn generates waves analytically for finite depth using
first-order (linear Airy) or first- plus second-order wave theory
[Sharma and Dean, 1981] with the option to include directional
spreading, but wave kinematics are only computed in the domain between
the flat seabed and SWL and no wave stretching or higher order wave
theories are included. Modeling unidirectional sea states is often
overly conservative in engineering design. Enabling the second-order
terms allows one to capture some of the nonlinearities of real surface
waves, permitting more accurate modeling of sea states and the
associated wave loads at the expense of greater computational effort
(mostly at HydroDyn initialization). The magnitude and frequency content
of second-order hydrodynamic loads can excite structural natural
frequencies, leading to greater ultimate and fatigue loads than can be
predicted solely using first-order theory. Sum-frequency effects are
important to the loading of stiff fixed-bottom structures and for the
springing and ringing analysis of TLPs. Difference-frequency (mean-drift
and slow-drift) effects are important to the analysis of compliant
structures, including the motion analysis and mooring loads of
catenary-moored floating platforms (spar buoys and semi-submersibles).

When modeling irregular sea states, we recommend that **WaveTMax** be
set to at least 1 hour (3600 s) and that **WaveDT** be a value in the
range between 0.1 and 1.0 s to ensure sufficient resolution of the wave
spectrum and wave kinematics. When HydroDyn is coupled to FAST,
**WaveDT** may be specified arbitrarily independently from the glue code
time step of FAST. (The wave kinematics and hydrodynamic loads will be
interpolated in time as necessary.)

Wave directional spreading is implemented in HydroDyn via the
equal-energy method, which assumes that the directional spreading
spectrum is the product of a frequency spectrum and a spreading function
i.e. *S*\ (*ω*,\ *β*) = *S*\ (*ω*)\ *D*\ (*β*). Directional spreading is
not permitted when using Newman’s approximation of the second-order
difference-frequency potential-flow loads.

When second-order terms are optionally enabled, the second-order terms
are calculated using the first-order wave-component amplitudes and extra
energy is added to the wave spectrum (at the difference and sum
frequencies). The second-order terms cannot be computed without also
including the first-order terms.

It is important to set proper wave cut-off frequencies to minimize
computational expense and to ensure that the wave kinematics and
hydrodynamic loads are realistic. HydroDyn gives the user six
user-defined cut-off frequencies—\ **WvLowCOff** and **WvHiCOff** for
the low- and high-frequency cut-offs of first-order wave components,
**WvLowCOffD** and **WvHiCOffD** for the low- and high-frequency
cut-offs of second-order difference-frequency wave components, and
**WvLowCOffS** and **WvHiCOffS** for low- and high-frequency cut-offs of
second-order sum-frequency wave components—none of which have default
settings. The second-order cut-offs apply directly to the physical
difference and sum frequencies, not the two individual first-order
frequency components of the difference and sum frequencies. Because the
second-order terms are calculated using the first-order wave-component
amplitudes, the second-order cut-off frequencies are used in conjunction
with the first-order cut-off frequencies. However, the second-order
cut-off frequencies are not used by Newman’s approximation of the
second-order difference-frequency potential-flow loads, which are
derived solely from first-order effects.

For the first-order wave-component cut-off frequencies, **WvLowCOff**
may be set lower than the low-energy limit of the first-order wave
spectrum to minimize computational expense. Setting a proper upper
cut-off frequency (**WvHiCOff**) also minimizes computational expense
and is important to prevent nonphysical effects when approaching of the
breaking-wave limit and to avoid nonphysical wave forces at high
frequencies (i.e., at short wavelengths) when using a strip-theory
solution.

When enabling second-order potential-flow theory, a setting of
**WvLowCOffD** = 0 is advised to avoid eliminating the mean-drift term
(second-order wave kinematics do not have a nonzero mean). **WvHiCOffD**
need not be set higher than the peak-spectral frequency of the
first-order wave spectrum (*ω\ p* = 2\ *π*/**WaveTp**) to minimize
computational expense. **WvLowCOffS** need not be set lower than the
peak-spectral frequency of the first-order wave spectrum (*ω\ p* =
2\ *π*/**WaveTp**) to minimize computational expense. Setting a proper
upper cut-off frequency (**WvHiCOffS**) also minimizes computational
expense and is important to (1) ensure convergence of the second-order
summations, (2) avoid unphysical "bumps" in the wave troughs, (3)
prevent nonphysical effects when approaching of the breaking-wave limit,
and (4) avoid nonphysical wave forces at high frequencies (i.e., at
short wavelengths) when using a strip-theory solution.

For all models with internally generated wave data, if you want to run
different time-domain incident wave realizations for given boundary
conditions (of significant wave height, and peak-spectral period, etc.),
you should change one or both wave seeds (**WaveSeed(1)** and
**WavedSeed(2)**) between simulations.

Wave elevations or full wave kinematics can also be generated externally
and used within HydroDyn.

**WaveMod** = 5 allows the use of externally generated wave-elevation
time series, which is useful if you want HydroDyn to simulate specific
wave transient events where the wave-elevation time series is known a
priori e.g. to match wave-elevation measurements taken from a wave tank
or open-ocean test. Internally, HydroDyn will compute an FFT of the
provided wave-elevation time series to store the amplitudes and phases
of each frequency component, and use those in place of a wave energy
spectrum and random seeds to internally derive the hydrodynamic loads in
the potential-flow solution or the wave kinematics used in the
strip-theory solution. The wave-elevation time series specified is
assumed to be of first order and long-crested, but is not checked for
physical correctness. The time series must be at least **WaveTMax** in
length and not less than the total simulation time and the time step
must match **WaveDT**. When second-order terms are optionally enabled,
the second-order terms are calculated using the wave-component
amplitudes derived from the provided wave-elevation time series and
extra energy is added to the wave energy spectrum (at the difference and
sum frequencies). Using higher order wave data may produce erroneous
results; alternatively, **WvLowCOff** and **WvHiCOff** can be used to
filter out energy outside of the first-order wave energy range. The
wave-elevation time series output by HydroDyn will only match the
specified time series identically if the second-order terms are disabled
and the cut-off frequencies are outside the range of wave energy.

**WaveMod** =6 allows the use of full externally generated wave
kinematics for use with the strip-theory solution (but not the
potential-flow solution), completely bypassing HydroDyn’s internal wave
models. This feature is useful if you want HydroDyn to make use of wave
kinematics data derived outside of HydroDyn a priori e.g. from a
separate numerical tool, perhaps bypassing some of HydroDyn’s internal
wave modeling limitations. To use this feature, it is the burden of the
user to generate wave kinematics data at each of HydroDyn’s time steps
and analysis nodes. HydroDyn will not interpolate the data; as such,
when HydroDyn is coupled to FAST, **WaveDT** must equal the glue code
time step of FAST. Before generating the wave kinematics data
externally, users should identify all of the internal analysis nodes by
running HydroDyn and generating the summary file—see :numref:`hd-summary-file`. The
fluid domain at each time step are specified by the use of numeric
values and nonnumeric strings in the wave data input files. The wave
kinematics data specified are not limited to the domain between a flat
seabed and SWL and may consider wave stretching, higher-order wave
theories, or an uneven seabed. The specified wave kinematics data are
not processed (filtered, etc.) or checked for physical correctness. The
wave kinematics output by HydroDyn should match the specified data
identically.

You can generate up to 9 wave elevation outputs (at different points on
the SWL plane) when HydroDyn is coupled to FAST or a large grid of wave
elevations when running HydroDyn standalone. While the second-order
effects are included when enabled, the wave elevations output from
HydroDyn will only include the second-order terms when the second-order
wave kinematics are enabled.

Strip-Theory Model Discretization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
A user will define the geometry of a structure modeled with strip theory
in HydroDyn using joints and members. Members in HydroDyn are assumed to
be straight circular (and possibly tapered) cylinders. Members can be
further subdivided using **MDivSize**, which HydroDyn will internally
use to subdivide members into multiple elements (and nodes). HydroDyn
may further refine the geometry at the free surface, flat seabed,
marine-growth region, and filled-fluid free surface.

.. TODO 7.5.2 is the theory section which does not yet exist. The rules HydroDyn uses for refinement may be found in Section 7.5.2.

Due to the exponential decay of hydrodynamic loads with depth, a higher
resolution near the water free surface is required to capture
hydrodynamic loading as waves oscillate about SWL. It is recommended,
for instance, that the HydroDyn discretization not exceed element
lengths of 0.5 m in the region of the free surface (5 to 10 m above and
below SWL), 1.0 m between 25 and 50 m depth, and 2.0 m in deeper waters.
When HydroDyn is coupled to SubDyn through FAST for the analysis of
fixed-bottom systems, it is recommended that the length ratio between
elements of HydroDyn and SubDyn not exceed 10 to 1.

.. _hd-domain-for-strip-theory:

Domain for Strip-Theory Hydrodynamic Load Calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Part of the automated geometry refinement mentioned in the above section
deals with splitting of input members into sub-elements such that both
of the resulting nodes at the element ends lie within the discrete
domains described in the following sections.

Distributed Loads
-----------------

Inertia, Added Mass, Buoyancy, Marine-Growth Weight, Marine-Growth Mass Inertia
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

These loads are generated at a node as long as **PropPot** = FALSE, the
*Z*-coordinate is in the range [–**WtrDpth**,\ **MSL2SWL**], and the
element the node is connected to is in the water. When **WaveMod** = 6,
the domain is determined by the use of numeric values and nonnumeric
strings in the wave data input files.

Viscous Drag
++++++++++++
These loads are generated at a node as long as the *Z*-coordinate is in
the range [–**WtrDpth**, **MSL2SWL**] and the element the node is
connected to is in the water. When **WaveMod** = 6, the domain is
determined by the use of numeric values and nonnumeric strings in the
wave data input files.

Filled Buoyancy, Filled Mass Inertia
++++++++++++++++++++++++++++++++++++
These loads are generated at a node as long as the *Z*-coordinate is in
the range [–**WtrDpth**, **FillFSLoc**] and the element the node is
connected to is in the filled fluid.

Lumped Loads
------------
Lumped loads at member ends (axial effects) are only calculated at
user-specified joints, and not at joints HydroDyn may automatically
create as part its solution process.
For example, if you want axial effects at a
marine-growth boundary, you must explicitly set a joint at that
location.

.. TODO 7.5.2 is the theory section which does not yet exist.
.. (see Section 7.5.2 for differences
.. between the input-file discretization and the simulation
.. discretization)

Added Mass, Inertia, Buoyancy
+++++++++++++++++++++++++++++
These loads are generated at a node as long as **PropPot** = FALSE and
the *Z*-coordinate is in the range [–**WtrDpth**,\ **MSL2SWL**]. When
**WaveMod** = 6, the domain is determined by the use of numeric values
and nonnumeric strings in the wave data input files.

Axial Drag
++++++++++
These loads are generated at a node as long as the *Z*-coordinate is in
the range [–**WtrDpth**,\ **MSL2SWL**]. When **WaveMod** = 6, the domain
is determined by the use of numeric values and nonnumeric strings in the
wave data input files.

Filled Buoyancy
+++++++++++++++
These loads are generated at a node as long as the *Z*-coordinate is in
the range [–**WtrDpth**,\ **FillFSLoc**]

Strip-Theory Hydrodynamic Coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The strip-theory solution of HydroDyn is dependent, among other factors,
on user-specified hydrodynamic coefficients, including viscous-drag
coefficients, **Cd**, added-mass coefficients, **Ca**, and
dynamic-pressure coefficients, **Cp**, for transverse and axial (**Ax**)
loads distributed along members and for axial lumped loads at member
ends (joints). There are no default settings for these coefficients in
HydroDyn. In general, these coefficients are dependent on many factors,
including Reynold’s number (Re), Keulegan-Carpenter number (KC), surface
roughness, substructure geometry, and location relative to the free
surface, among others. In practice, the coefficients are (1) selected
from tables derived from measurements of flow past cylinders, (2)
calculated through high-fidelity computational fluid dynamics (CFD)
solutions, or (3) tuned to match experimental results. A value of 1.0 is
a plausible guess for all coefficients in the absence of any other
information.

While the strip-theory solution assumes circular cross sections, the
hydrodynamic coefficients can include shape corrections; however, there
is no distinction made in HydroDyn between different transverse
directions.

Please note that added-mass coefficients in HydroDyn influence both the
added-mass loads and the scattering component of the fluid-inertia
loads. For the coefficients associated with transverse loads distributed
along members, note that :math:`C_{P} + C_{A} = C_{M}`,
the inertia coefficient. For the
distributed loads along members, there are separate set of hydrodynamic
coefficients both with and without marine growth (**MG**).

Impact of Substructure Motions on Loads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In general, HydroDyn assumes that structural motions of the substructure
are small, such that (1) small-angle assumptions apply to structural
rotations, (2) the frequency-to-time-domain-based potential-flow
solution can be split into uncoupled hydrostatic, radiation, and
diffraction solutions, and (3) the hydrodynamic loads dependent on wave
kinematics (both from diffraction loads in the potential-flow solution
and from the fluid-inertia and viscous-drag loads in the strip-theory
solution) can be computed using wave kinematics solved at the
undisplaced position of the substructure (the wave kinematics are not
recomputed at the displaced position). Nevertheless, HydroDyn uses the
substructure motions in the following calculations:

*  The structural displacements of the WRP are used in the calculation
   of the hydrostatic loads (i.e., the change in buoyancy with
   substructure displacement) in the potential-flow solution.

*  The structural velocities and accelerations of the WRP are used in
   the calculation of the wave-radiation loads (i.e., the radiation
   memory effect and added mass) in the potential-flow solution.

*  The structural displacements and velocities of the WRP are used in
   the calculation of the additional platform loads (via the Platform
   Additional Stiffness and Damping).

*  The structural velocities of the substructure nodes are used in the
   calculation of the viscous-drag loads in the strip-theory solution
   (e.g., the relative form of Morison’s equation is applied).

*  The structural accelerations of the substructure nodes are used in
   the calculation of the added-mass, marine-growth mass inertia, and
   filled-fluid mass inertia loads in the strip-theory solution.

*  When coupled to FAST, the hydrodynamic loads computed by HydroDyn are
   applied to the displaced position of the substructure (i.e., the
   displaced platform in ElastoDyn and/or the displaced substructure
   in SubDyn), but are based on wave kinematics at the undisplaced
   position.

Platform Additional Stiffness and Damping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
HydroDyn allows the user to apply additional loads to the platform (in
addition to other hydrodynamic terms calculated by HydroDyn), by
including a 6x1 static load vector (preload) (**AddF0**), a 6x6 linear
restoring matrix (**AddCLin**), a 6x6 linear damping matrix
(**AddBLin**), and a 6x6 quadratic drag matrix (**AddBQuad**). These
terms can be used, e.g., to model a linearized mooring system, to
augment strip-theory members with a linear hydrostatic restoring matrix
(see :numref:`hd-modeling-hydrostatic-restoring-strip-theory`), or to
"tune" HydroDyn to match damping to
experimental results, such as free-decay tests. While likely most useful
for floating systems, these matrices can also be used for fixed-bottom
systems; in both cases, the resulting load is applied at the WRP, which
when HydroDyn is coupled to FAST, get applied to the platform in
ElastoDyn (bypassing SubDyn for fixed-bottom systems).

Fixed-Bottom Substructures
~~~~~~~~~~~~~~~~~~~~~~~~~~
When modeling a fixed-bottom system, the use of a strip-theory (Morison)
only model is recommended. When HydroDyn is coupled to FAST, SubDyn is
used for the substructure structural dynamics.

All members that are embedded into the seabed (e.g., through piles or
suction buckets) must have a joint that is located below the water
depth. For example, if the water depth is set to 20 m, and you are
modeling a fixed-bottom monopile, then the bottom-most joint needs to
have a *Z*-coordinate such that m. This configuration avoids having
HydroDyn apply static pressure loads on the bottom of the structure.

Gravity-based foundations should be modeled such that the lowest
joint(s) are located exactly at the prescribed water depth. In other
words, the lowest *Z*-coordinate should be set to m if the water depth
is set to 20 m. This configuration allows for static pressure loads to
be applied at the bottom of the gravity-base structure.

Floating Platforms
~~~~~~~~~~~~~~~~~~
When modeling a floating system, you may use potential-flow theory only,
strip-theory (Morison) only, or a hybrid model containing both.

Potential-flow theory based on frequency-to-time-domain transforms is
enabled when **PotMod** is set to 1. In this case, you must run WAMIT
(or equivalent) in a pre-processing step and HydroDyn will use the WAMIT
output files—see :numref:`hd-modeling-floating-systems-potential-flow`
for guidance. For a potential-flow-only
model, do not create any strip-theory joints or members in the input
file. The WAMIT model should account for all of the members in the
floating substructure, and Morison’s equation is neglected in this case.

For a strip-theory-only model, set **PotMod** to FALSE and create one or
more strip-theory members in the input file. Marine growth and nonzero
**MSL2SWL** (the offset between still-water and mean-sea level) may only
be included in strip-theory-only models.

A hybrid model is formed when both **PotMod** is TRUE and you have
defined one or more strip-theory members. The potential-flow model
created can consider all of the Morison members in the floating
substructure, or just some. Specify whether certain members of the
structure are considered in the potential-flow model by setting the
**PropPot** flag for each member.
The state of the **PropPot** flag for a given member determines which
components of the strip-theory equations are applied.

.. TODO 7.5.1 is the theory section which does not yet exist.
.. As detailed in Section 7.5.1,

When using either the strip-theory-only or hybrid approaches, filled
fluid (flooding or ballasting) may be added to the strip-theory members.
Also, the hydrostatic restoring matrix must be entered manually for the
strip-theory members—see :numref:`hd-modeling-hydrostatic-restoring-strip-theory` for guidance.

Please note that current-induced water velocity only induces
hydrodynamic loads in HydroDyn through the viscous-drag terms (both
distributed and lumped) of strip-theory members. Current is not used in
the potential-flow solution. Thus, modeling the effects of current
requires the use of a strip-theory-only or hybrid approach.

Undisplaced Position for Floating Systems
-----------------------------------------
The HydroDyn model (geometry, etc.) is defined about the undisplaced
position of the substructure. For floating systems, it is important for
solution accuracy for the undisplaced position to coincide with the
static-equilibrium position in the platform-heave (vertical) direction
in the absence of loading from wind, waves, and current. As such, the
undisplaced position of the substructure should be defined such that the
external buoyancy from displaced water balances with the weight of the
system (including the weight of the rotor-nacelle assembly, tower and
substructure) and mooring system pretension following the equation
below. In this equation, is the water density, is gravity, is the
undisplaced volume of the floating platform (found in the HydroDyn
summary file), is the total mass of the system (found in the ElastoDyn
summary), and is the mooring system pretension (found in e.g. the MAP
summary file). The effects of marine growth, filled fluid (flooding
and/or ballasting), and the additional static force (**AddFX0**) should
also be taken into consideration in this force balance, where
appropriate.

.. math::
   :label: FloatingForceBalance

   \rho g V_{0} - m_{Total} g - T_{Mooring} = 0

Initial Conditions for Floating Systems
---------------------------------------
Because the initial conditions used for dynamic simulations typically
have an effect on the response statistics during the beginning of the
simulation period, an appropriate amount of initial data should be
eliminated from consideration in any post-processing analysis. This
initial condition solution is more important for floating offshore wind
turbines because floating systems typically have long natural periods of
the floating substructure and low damping. The appropriate time to
eliminate should be chosen such that initial numeric transient effects
have sufficiently decayed and the floating substructure has reached a
quasi-stationary position. To decrease this initial time in each
simulation, it is suggested that the initial conditions of the model
(especially blade-pitch angle, rotor speed, substructure surge, and
substructure pitch in ElastoDyn) be initialized according to the
specific prevalent wind, wave, current, and operational conditions.

.. _hd-modeling-hydrostatic-restoring-strip-theory:

Hydrostatic Restoring for Strip-Theory Members of Floating Systems
------------------------------------------------------------------
One notable absence from the list calculations in HydroDyn that make use
of substructure motions—see :numref:`hd-domain-for-strip-theory`—is that the substructure
buoyancy in the strip-theory solution is not recomputed based on the
displaced position of the substructure. While the change in buoyancy is
likely negligible for fixed-bottom systems, for floating systems modeled
using a strip-theory solution, the change in buoyancy with displacement
is likely important and should not be neglected. In this latter case,
the user should manually calculate the 6x6 linear hydrostatic restoring
matrix associated with the strip-theory members and enter this as the
additional linear restoring (stiffness) matrix, **AddCLin**. (The static
buoyancy of the strip-theory members is automatically calculated and
applied within HydroDyn.)

In its most general form, the 6x6 linear hydrostatic restoring matrix of
a floating platform is given by the equation below.

.. math::
   :label: HydrostraticRestoringMatrix

   \text{AddCLin} =
   \left[
      \begin{array}{cccccc}
         0 & 0 & 0 & 0 & 0 & 0 \\
         0 & 0 & 0 & 0 & 0 & 0 \\
         0 & 0 & \rho g A_{0} & \rho g \iint_{A_{0}} ydA & -\rho g \iint_{A_{0}} xdA & 0 \\
         0 & 0 & \rho g \iint_{A_{0}} ydA & \rho g \iint_{A_{0}} y^2dA + \rho g V_{0} z_{b} - m_{mg}gz_{mg} - m_{f}gz_{f} & -\rho g \iint_{A_{0}} xydA & -\rho g V_{0} x_{b} + m_{mg}gx_{mg} + m_{f}gx_{f} \\
         0 & 0 & -\rho g \iint_{A_{0}} xdA & -\rho g \iint_{A_{0}} xydA & \rho g \iint_{A_{0}} x^2dA + \rho g V_{0} z_{b} - m_{mg}gz_{mg} - m_{f}gz_{f} & -\rho g V_{0} y_{b} + m_{mg}gy_{mg} + m_{f}gy_{f} \\
         0 & 0 & 0 & 0 & 0 & 0
      \end{array}
   \right]

where:

* :math:`\rho` is water density (kg/m\ :sup:`3`)

* :math:`g` is gravity (m/s\ :sup:`2`)

* :math:`A_{0}` is undisplaced waterplane area of platform (m\ :sup:`2`)

* :math:`V_{0}` is undisplaced volume of platform (m\ :sup:`3`)

* :math:`(x_{b}, y_{b}, z_{b})` is coordinates of the center of buoyancy of the undisplaced platform (m)

* :math:`m_{mg}` is total mass of marine growth (kg)

* :math:`(x_{mg}, y_{mg}, z_{mg})` is coordinates of the center of mass of the undisplaced marine growth mass (m)

* :math:`m_{f}` is total mass of ballasting/flooding (kg)

* :math:`(x_{f}, y_{f}, z_{f})` is coordinates of the center of mass of the undisplaced filled fluid (flooding or ballasting) mass (m)

The equation above can be simplified when the floating platform has one
or more planes of symmetry. That is,
:math:`\iint_{A_{0}} ydA = 0`,
:math:`\iint_{A_{0}} xydA = 0`,
:math:`y_{b} = 0`,
:math:`y_{mg} = 0`,
:math:`y_{f} = 0`,
and if the :math:`x-z` plane of the platform is a symmetry plane.
Likewise,
:math:`\iint_{A_{0}} xdA = 0`,
:math:`\iint_{A_{0}} xydA = 0`,
:math:`x_{b} = 0`,
:math:`x_{mg} = 0`,
:math:`x_{f} = 0`,
and if the :math:`y-z` plane of the platform is a symmetry plane.

The undisplaced coordinates of the center of buoyancy,
:math:`(x_{b}, y_{b}, z_{b})`, center of
marine-growth mass, :math:`(x_{mg}, y_{mg}, z_{mg})`, and center of
filled-fluid mass, :math:`(x_{f}, y_{f}, z_{f})`, are in the
global inertial-frame coordinate system. Most of these parameters can be
derived from data found in the HydroDyn summary file. While the equation
above makes use of several area integrals, the integrals can often be
easily estimated by hand for platforms composed of one or more circular
members piercing the waterplane (still-water free surface).

The waterplane area of the undisplaced platform, :math:`A_{0}`, affects the
hydrostatic load because the displaced volume of the fluid changes with
changes in the platform displacement. Similarly, the location of the
center of buoyancy of the platform affects the hydrostatic load because
its vector position changes with platform displacement and because the
cross product of the buoyancy force with the vector position produces
hydrostatic moments about the WRP. :math:`A_{0}`, :math:`V_{0}`, and
:math:`(x_{b}, y_{b}, z_{b})` should be based on the
external volume of the platform, including marine-growth thickness. The
marine-growth mass and filled-fluid mass also have a direct effect of
the hydrostatic restoring because of the moments produced about the WRP.

In classical marine hydrostatics, the effects of body weight are often
lumped with the effects of hydrostatics when defining the
hydrostatic-restoring matrix; for example, when it is defined in terms
of metacentric heights. However, when HydroDyn is coupled to FAST, the
body-weight terms (other than the marine-growth and filled-fluid mass
within HydroDyn) are automatically accounted for by ElastoDyn, and so,
are not included here.

.. _hd-modeling-floating-systems-potential-flow:

Floating Systems Modeled with Potential Flow
--------------------------------------------
Frequency-dependent hydrodynamic coefficients are needed before running
the potential-flow solution in HydroDyn using **PotMod** = 1. An
external pre-processing tool should be used to generate the appropriate
frequency-dependent hydrodynamic coefficients. The naming in this manual
has focused on WAMIT :cite:`LeeNewman:2006`, but other frequency-domain wave-body
interaction panel codes can be used that produce similar data. However,
in the end, the WAMIT format is what is expected by HydroDyn.

For the first-order potential-flow solution, HydroDyn requires data from
the WAMIT files with *.1, .3*, and *.hst* extensions. When creating
these files, one should keep in mind:

*  The *.1* file must contain the 6×6 added-mass matrix at infinite
   frequency (period = zero). Additionally, the *.1* file must contain
   the 6×6 damping matrix over a large range from low frequency to high
   frequency (the damping should approach zero at both ends of the
   range). A range of 0.0 to 5.0 rad/s with a discretization of 0.05
   rad/s is suggested.

*  The .\ *3* file must contain the first-order wave-excitation
   (diffraction) loads (3 forces and 3 moments) per unit wave amplitude
   across frequencies and directions where there is wave energy. A range
   of 0.0 to 5.0 rad/s with a discretization of 0.05 rad/s is suggested
   and the direction should be specified across the desired range—the
   full direction range of (-180 to 180] degrees with a discretization
   of 10 degrees is suggested. While the .\ *3* file contains both the
   magnitude/phase and real/imaginary components of the first-order
   wave-excitation loads, only the latter are used by HydroDyn.

*  The .\ *hst* file should account for the restoring provided by
   buoyancy, but not the restoring provided by body mass or moorings.
   (The hydrostatic file is not frequency dependent.) An important thing
   to keep in mind is that the pitch and roll restoring of a floating
   body depends on the vertical distance between the center of buoyancy
   and center of mass of the body. In WAMIT, the vertical center of
   gravity (VCG) is used to determine the pitch and roll restoring
   associated with platform weight, and WAMIT will include these effects
   in the restoring matrix that it outputs (the *.hst* file). However,
   the ElastoDyn module of FAST intrinsically accounts for the platform
   weight’s influence on the pitch and roll restoring if the platform
   weight and center-of-mass location are defined appropriately. To
   avoid double booking these terms, it is important to neglect these
   terms in WAMIT. This can be achieved by setting VCG to zero when
   solving the first-order problem in WAMIT.

The second-order WAMIT files only need to pre-calculated if a
second-order potential-flow option is enabled in HydroDyn. For the
second-order mean-drift solution, or for Standing et al.’s extension to
Newman’s approximation to the mean- and slow-drift solution, HydroDyn
requires WAMIT files with .\ *7*, .\ *8*. .\ *9*, .\ *10d*, .\ *11d*, or
.\ *12d* extensions. For the second-order full difference-frequency
solution of the mean- and slow-drift terms, HydroDyn requires WAMIT
files with .\ *10d*, .\ *11d*, or .\ *12d* extension. For the
second-order full sum-frequency solution, HydroDyn requires WAMIT files
with .\ *10s*, .\ *11s*, or .\ *12s* extensions. When creating any of
these files, one should keep in mind:

*  The second-order frequency-domain solution is dependent on
   first-order body motions, whose accuracy is impacted by properly
   setting the 6×6 rigid-body mass matrix and center of gravity of the
   complete floating wind system and the 6×6 mooring system restoring
   matrix. So, while the body center of gravity and mooring stiffness
   should be zeroed when creating the first-order WAMIT files, they
   should not be zeroed when creating the second-order WAMIT files.
   (Thus, obtaining the first-order and second-order WAMIT files
   requires distinct WAMIT runs.)

*  The .\ *7*, .\ *8*, and .\ *9* files contain the diagonal of the
   difference-frequency QTF, based on the first-order potential-flow
   solution. The files contain the second-order mean-drift loads (3
   forces and 3 moments) per unit wave amplitude squared at each
   first-order wave frequency and pair of wave directions, across a
   range of frequencies and a range of direction pairs. While the
   .\ *7*, .\ *8*, and .\ *9* files contains both the magnitude/phase
   and real/imaginary components of the second-order wave-excitation
   loads, only the latter are used by HydroDyn.

*  The *10d*, .\ *11d*, and .\ *12d*, or .\ *10s*, .\ *11s*, and
   .\ *12s* files contain the full difference- and sum-frequency QTFs,
   respectively, based on the first-order or first- plus second-order
   potential-flow solutions. The files contain the second-order
   wave-excitation (diffraction) loads (3 forces and 3 moments) per unit
   wave amplitude squared at each pair of first-order wave frequencies
   and directions, across a range of frequency and direction pairs.
   While the *10d*, .\ *11d*,.\ *12d*, .\ *10s*, .\ *11s*, and .\ *12s*
   files contains both the magnitude/phase and real/imaginary components
   of the second-order wave-excitation loads, only the latter are used
   by HydroDyn.

*  The frequencies and directions in the WAMIT files do not need to be
   evenly spaced.

*  The discretization of the first set of directions does not need to be
   the same as the discretization of the second set of directions;
   however, the matrix of direction pairs must be fully populated (not
   sparse). Both sets of directions should span across the desired
   range—the full direction range of (-180 to 180] degrees with a
   discretization of 10 degrees is suggested.

*  The frequencies should span the range where there is first-order wave
   energy and the frequency discretization should be such that the
   differences and sums between pairs of frequencies span the range
   where there is second-order wave energy. A range of 0.25 to 2.75
   rad/s with a discretization of 0.05 rad/s is suggested.

*  Second-order hydrodynamic theory dictates that difference-frequency
   QTFs are conjugate symmetric between frequency pairs and
   sum-frequency QTFs are symmetric between frequency pairs. Due to this
   symmetry, the QTFs (the *10d*, .\ *11d*, or .\ *12d*, .\ *10s*,
   .\ *11s*, and .\ *12s* files) may be upper triangular, lower
   triangular, a mix of upper and lower triangular terms, or full;
   however, after applying the symmetry, the matrix of frequency pairs
   must be fully populated (not sparse). When an element of the QTF is
   supplied together with its symmetric pairing, HydroDyn will warn the
   user if the QTF is not properly symmetric.
