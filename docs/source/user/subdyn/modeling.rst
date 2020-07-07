.. _sd_modeling-considerations:

Modeling Considerations
=======================

SubDyn was designed as a flexible tool for modeling a wide range of
substructures for both land-based and offshore applications. This
section provides some general guidance to help construct models that are
compatible with SubDyn.

Please refer to the theory in Section 6 for detailed information about
SubDyn’s coordinate systems, and the theoretical approach we have
followed in SubDyn.

Model Discretization
--------------------

SubDyn allows for the specification of arbitrary multimember structure
geometries. The user defines the geometry of a structure in SubDyn using
joints and members. Specifically, the user specifies a list of joints
that represent the endpoints of beams, and the connectivity between one
or more members at each joint. Members and their cross-sectional
properties are then defined between two joints. Members can be further
subdivided into multiple (**NDiv**) elements to increase the model
resolution. Nodes, where the numerical calculations take place, are
located at the endpoints of each element. To keep the mesh as uniform as
possible when using **NDiv**, the initial member definition should
also have a roughly uniform mesh. For tapered members, we recommend
setting **NDiv** > 1. Improper discretization of the members may
decrease the accuracy of the model.

When SubDyn is coupled to FAST, the joints and members need not match
between HydroDyn and SubDyn—FAST’s mesh-mapping utility handles the
transfer of motion and loads across meshes in a physically relevant
manner :cite:`sprague2014`, but consistency between the joints and
members in HydroDyn and SubDyn is advised.

For offshore applications, because of the exponential decay of
hydrodynamic loads with depth, HydroDyn requires higher resolution near
the water free surface to properly capture loads as waves oscillate
about the still water level (SWL). We recommend that the HydroDyn
discretization not exceed element lengths of 0.5 m in the region of the
free surface (5 to 10 m above and below SWL), 1.0 m between 25- and 50-m
depth, and 2.0 m in deeper waters.

When SubDyn is hydro-elastically coupled to HydroDyn through FAST for
the analysis of fixed-bottom offshore systems, we recommend that the
length ratio between elements of SubDyn and HydroDyn not exceed 10 to 1.
As such, we recommend that the SubDyn discretization not exceed element
lengths of 5 m in the region of the free surface, 10 m down to 25- to
50-m depth, and 20 m in deeper waters. These are not absolute rules, but
rather a good starting point that will likely require refinement for a
given substructure. Additional considerations for SubDyn discretization
include aspects that will impact structural accuracy, such as member
weight, substructure modes and/or natural frequencies, load transfer,
tapered members, and so on.

Members in SubDyn are assumed to be straight circular (and possibly
tapered) cylinders. The use of more generic cross-sectional shapes will
be considered in a future release.

Foundations
-----------

There are two methods that can be used to model foundation flexibility
or soil-structure interaction in SubDyn. The first method makes us of
the SSI stiffness and mass matrices at the partially restrained bottom
joints as described in Sections 3.3.4, 3.4, and 6. The second method
mimics the flexibility of the foundation through the apparent (or
effective) fixity (AF) length approach, which idealizes a pile as a
cantilever beam that has properties that are different above and below
the mudline. The beam above the mudline should have the real properties
(i.e., diameter, thickness, and material) of the pile. The beam below
the mudline is specified with effective properties and a fictive length
(i.e., the distance from the mudline to the cantilevered base) that are
tuned to ensure that the overall response of the pile above the mudline
is the same as the reality. The response can only be identical under a
particular set of conditions; however, it is common for the properties
of the fictive beam to be tuned so that the mudline displacement and
rotation would be realistic when loaded by a mudline shear force and
bending moment that are representative of the loading that exists when
the offshore wind turbine is operating under normal conditions.

Note that in HydroDyn, all members that are embedded into the seabed
(e.g., through piles or suction buckets) must have a joint that is
located below the water depth. In SubDyn, the bottom joint(s) will be
considered clamped or partially restrained and therefore need not be
located below the seabed when not applying the AF approach. For example,
if the water depth is set to 20 m, and the user is modeling a
fixed-bottom monopile with a rigid foundation, then the bottom-most
joint in SubDyn can be set at *Z* = -20 m; HydroDyn, however, needs to
have a Z-coordinate such that *Z* < -20 m. This configuration avoids
HydroDyn applying static and dynamic pressure loads from the water on
the bottom of the structure. When the AF approach is applied, the
bottom-most joint in SubDyn should be set at *Z* < -20 m.


Member Overlap
--------------

As mentioned earlier, the current version of SubDyn is incapable of
treating the overlap of members at the joints, resulting in an
overestimate of the mass and potentially of the structure stiffness. One
strategy to overcome this shortcoming employs virtual members to
simulate the portion of each member within the overlap at a joint. The
virtual members should be characterized by low self-mass and high
stiffness. This can be achieved by introducing virtual joints at the
approximate intersection of the finite-sized members, and then
specifying additional members from these new joints to the original
(centerline) joints. The new virtual members then use reduced material
density and increased Young’s and shear moduli. Care is advised in the
choice of these parameters as they may render the system matrix
singular. Inspection of the eigenvalue results in the summary file
should confirm whether acceptable approximations have been achieved.

.. _TowerTurbineCpling:

Substructure Tower/Turbine Coupling 
-----------------------------------

When SubDyn is coupled to FAST, the 6 DOFs of the platform in ElastoDyn
must be enabled to couple loads and displacements between the turbine
and the substructure. The platform reference-point coordinates in
ElastoDyn should also be set equal to the TP reference-point’s
coordinates (commonly indicating either the tower-base flange location,
or TP centroid, or TP center of mass) that the user may have set in the
stand-alone mode for checking the SubDyn model. A rigid connection
between the SubDyn interface joints and TP reference point (:math:`{\equiv}` platform
reference point) is assumed.

For full lattice support structures or other structures with no
transition piece, the entire support structure up to the yaw bearing may
be modeled within SubDyn. Modeling the tower in SubDyn as opposed to
ElastoDyn, for example, allows the ability to include more than the
first two fore-aft and side-to-side bending modes, thus accounting for
more general flexibility of the tower and its segments; however, for
tubular towers, the structural model in ElastoDyn tends to be more
accurate because ElastoDyn considers geometric nonlinearities not
treated in SubDyn. When modeling full-lattice towers using SubDyn, the
platform reference point in ElastoDyn can be located at the yaw bearing;
in this case, the tower-bending DOFs in ElastoDyn should be disabled.

If FAST is run with SubDyn but not HydroDyn, the water depth will be
automatically set to 0 m. This will influence the calculation of the
reaction loads. Reactions are always provided at the assumed mudline,
therefore, they would not be correctly located for an offshore turbine
as a result. Thus, it is recommended that HydroDyn always be enabled
when modeling bottom-fixed offshore wind turbines.

ElastoDyn also needs tower mode shapes specified (coefficients of
best-fit sixth-order polynomials), derived using appropriate tower-base
boundary conditions. They can be derived with an appropriate software
(finite-element analysis, energy methods, or analytically) and by making
use of the SubDyn-derived equivalent substructure stiffness and mass
matrices (the **KBBt** and **MBBt** matrices found in the SubDyn summary
file) to prescribe the boundary conditions at the base of the tower.

For instance, using NREL’s BModes software, the SubDyn-obtained matrices
can be used in place of the hydrodynamic stiffness (**hydro\_K**) and mass
matrices (**hydro\_M**) (**mooring\_K** can be set to zero). By setting
the **hub\_conn** boundary condition to two (free-free), BModes will
calculate the mode shapes of the tower when tower cross-sectional
properties are supplied. To obtain eigenmodes that are compatible with
the FAST modal treatment of the tower (i.e., no axial or torsional modes
and no distributed rotational-inertia contribution to the eigenmodes),
the tower-distributed properties should be modified accordingly in
BModes (e.g., by reducing mass moments of inertia towards zero and by
increasing torsional and axial stiffness while assuring convergence of
the results; see also
`https://wind.nrel.gov/forum/wind/viewtopic.php?f=4&t=742 <https://wind.nrel.gov/forum/wind/viewtopic.php?f=4&t=742>`__).

The rotational inertia of the undeflected tower about its centerline is
not currently accounted for in ElastoDyn. Thus, when the nacelle-yaw DOF
is enabled in ElastoDyn there will not be any rotational inertia of the
platform-yaw DOF (which rotates the tower about its centerline) when
both the platform-yaw inertia in ElastoDyn is zero and the tower is
undeflected. To avoid a potential division-by-zero error in ElastoDyn
when coupled to SubDyn, we recommend setting the platform-yaw inertia
(**PtfmYIner**) in ElastoDyn equal to the total rotational inertia of
the undeflected tower about its centerline. Note that the platform mass
and inertia in ElastoDyn can be used to model heavy and rigid transition
pieces that one would not want to model as a flexible body in either the
ElastoDyn tower or SubDyn substructure models.

***Damping of the Guyan modes:***

There are three ways to specify the damping associated with the motion
of the interface node.

1. SubDyn Guyan damping matrix using Rayleigh damping
2. SubDyn Guyan damping matrix using user defined 6x6 matrix
3. HydroDyn additional linear damping matrix (**AddBLin**)

The specificaiton of the Guyan damping matrix in SubDyn is discussed in :numref:`SD_DampingSpecifications`.


**Old:**

The C-B method assumes no damping for the interface modes. This is
equivalent to having six undamped rigid-body DOFs at the TP reference
point in the absence of aerodynamic or hydrodynamic damping. Experience
has shown that negligible platform-heave damping can cause numerical
problems when SubDyn is coupled to FAST. One way to overcome this
problem is to augment overall system damping with an additional linear
damping for the platform-heave DOF. This augmentation can be achieved
quite easily by calculating the damping from Eq. :eq:`damping` and specifying this
as the (3,3) element of HydroDyn’s additional linear damping matrix,
**AddBLin**. Experience has shown that a damping ratio of 1% of
critical (:math:`{\zeta=0.01}`) is sufficient. In Eq. :eq:`damping`, :math:`{K_{33}^{(SD)}}` is the equivalent heave stiffness
of the substructure (the (3,3) element of the **KBBt** (i.e., :math:`{\tilde{K}_{BB}}`) matrix
found in the SubDyn summary file, see also Section 6), :math:`{M_{33}^{(SD)}}` is the equivalent
heave mass of the substructure (the (3,3) element of the **MBBt**
(i.e., :math:`{\tilde{M}_{BB}}`) matrix found in the SubDyn summary file, see also Section 6),
and :math:`{M^{(ED)}}` is the total mass of the rotor, nacelle, tower, and TP (found in the
ElastoDyn summary file).

.. math:: :label: damping
   	
   	C_{33}^{(HD)} = 2 \zeta \sqrt{ K_{33}^{(SD)} \left( M_{33}^{(SD)}+M^{(ED)} \right)}  


To minimize extraneous excitation of the platform-heave DOF, it is
useful to set the initial platform-heave displacement to its natural
static-equilibrium position, which can be approximated by Eq. :eq:`ptfmheave`, where
is the magnitude of gravity. *PtfmHeave* from Eq. :eq:`ptfmheave` should be
specified in the initial conditions section of the ElastoDyn input file.

.. math:: :label: ptfmheave
   	
   	PtfmHeave = -\dfrac{ \left( M_{33}^{(SD)}+M^{(ED)} \right) g}{K_{33}^{(SD)}}   
   




Self-Weight Calculations
------------------------

SubDyn will calculate the self-weight of the members and apply
appropriate forces and moments at the element nodes. Lumped masses will
also be considered as concentrated gravity loads at prescribed joints.
The array of self-weight forces can be seen in the summary file if the
code is compiled with DEBUG compiler directives. In general, SubDyn
assumes that structural motions of the substructure are small, such that
(1) small-angle assumptions apply to structural rotations and (2) the
so-called P- :math:`{\Delta}` effect is negligible, and therefore undeflected node
locations are used for self-weight calculations.

Note On Other Load Calculations
-------------------------------

When SubDyn is coupled to HydroDyn through FAST, the hydrodynamic loads,
which include buoyancy, marine-growth weight, and wave and current
loads, will be applied to the effective, deflected location of the nodes
by the mesh-mapping routines in the glue code. Those loads, however, are
based on wave kinematics at the undeflected position (see Jonkman et al.
2014 for more information).

.. _CBguide:

Craig-Bampton Guidelines 
------------------------

When SubDyn is coupled with FAST, it is important to choose a sufficient
number of C-B modes, ensuring that the vibrational modes of the coupled
system are properly captured by the coupled model. We recommend that all
modes up to at least 2-3 Hz be captured; wind, wave, and turbine
excitations are important for frequencies up to 2-3 Hz. Eigenanalysis of
the linearized, coupled system will make checking this condition
possible and aid in the selection of the number of retained modes;
however, the linearization process has yet to be implemented in FAST v8.
Until full-system linearization is made available, experience has shown
that it is sufficient to enable all C-B modes up to 10 Hz (the natural
frequencies of the C-B modes are written to the SubDyn summary file). If
SIM (see Section :numref:`sim`) is not enabled, in addition to capturing physical
modes up to a given frequency, the highest C-B mode must include the
substructure axial modes so that gravity loading from self-weight is
properly accounted for within SubDyn. This inclusion likely requires
enabling a high number of C-B modes, reducing the benefit of the C-B
reduction. Thus, we recommend employing the C-B reduction with SIM
enabled. Because of the fixed-fixed treatment of the substructure
boundary conditions in the C-B reduction, the C-B modes will always have
higher natural frequencies than the physical modes.

Integration Time Step Guidelines
--------------------------------

Another consideration when creating SubDyn input files is the time step
size. SubDyn offers three explicit time-integrators --- the fourth-order
Runge-Kutta (RK4), fourth-order Adams-Bashforth (AB4), fourth-order
Adams-Bashforth-Moulton (ABM4) methods --- and the implicit second-order
Adams-Moulton (AM2) method. Users have the option of using the global
time step from the glue code or an alternative SubDyn-unique time step
that is an integer multiple smaller than the glue-code time step.
It is essential that a small enough time step is used to ensure solution
accuracy (by providing a sufficient sampling rate to characterize all
key frequencies of the system), numerical stability of the selected
explicit time-integrator, and that the coupling with FAST is numerically
stable.

For the RK4 and ABM4 methods, we recommend that the SubDyn time step
follow the relationship shown in Eq. :eq:`dtmax`, where :math:`{f_{max}}` is the higher of (1) the
highest natural frequency of the retained C-B modes and (2) the highest
natural frequency of the physical modes when coupled to FAST. Although
the former can be obtained from the SubDyn summary file, the latter is
hard to estimate before the full-system linearization of the coupled
FAST model is realized. Until then, experience has shown that the
highest physical mode when SubDyn is coupled to FAST is often the
platform-heave mode of ElastoDyn, with a frequency given by Eq. :eq:`freq`,
where the variables are defined in Section 5.3.

.. math:: :label: dtmax
   	
   	dt_{max} = \dfrac{1}{10 f_{max}} 
   	
.. math:: :label: freq 
   	
   	f= \dfrac{1}{2\pi} \sqrt{\dfrac{K_{33}^{(SD)}}{ M_{33}^{(SD)}+M^{(ED)}}}  
   	
For the AB4 method, the recommended time step is half the value given by
Eq. :eq:`dtmax`.

For AM2, being implicit, the required time step is not driven by natural
frequencies within SubDyn, but should still be chosen to ensure solution
accuracy and that the coupling to FAST is numerically stable.

