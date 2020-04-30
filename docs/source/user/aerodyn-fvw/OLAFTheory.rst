.. _OLAF_Theory:

OLAF Theory
===========

This section details the FVW method and provides an overview of the
computational method, followed by a brief explanation of its integration
with OpenFAST.

.. _sec:FVW:

Free Vortex Wake Model
----------------------

.. _sec:vorticityformulation:

Introduction - Vorticity formulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The vorticity equation for incompressible homogeneous flows in the
absence of non-conservative force is given by
Eq. :eq:`eq:vorticityconservationincompr`

.. math::
   \begin{aligned}
   \frac{d\vec{\omega}}{dt} = \frac{\partial\vec{\omega}}{\partial{t}} + \underbrace{(\vec{u} \cdot \nabla)}_{\text{convection}}\vec{\omega} = \underbrace{(\vec{\omega}\cdot\nabla)\vec{u}}_{\text{strain}} +\underbrace{\nu\Delta\vec{\omega}}_{\text{diffusion}}
   \end{aligned}
   :label: eq:vorticityconservationincompr


Here, :math:`\vec{\omega}` is the vorticity, :math:`\vec{u}` is the
velocity, and :math:`\nu` is the viscosity. In free vortex wake methods,
the vorticity equation is used to describe the evolution of the wake
vorticity. Different approximations are introduced to ease its
resolution, such as projecting the vorticity onto a discrete number of
vortex elements (here vortex filaments), and separately treating the
convection and diffusion steps, known as viscous-splitting. Several
complications arise from the method, in particular, the discretization
requires a regularization of the vorticity field (or velocity field) to
ensure a smooth approximation.

The forces exerted by the blades onto the flow are expressed in
vorticity formulation as well. This vorticity is bound to the blade and
has a circulation associated with the lift force. A lifting-line
formulation is used here to model the bound vorticity.

The different models of the free vortex code implemented are described
in the following sections.

.. _sec:discretization:

Discretization - Projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The numerical method uses a finite number of states to model the
continuous vorticity distribution. To achieve this, the vorticity
distribution is projected onto basis function which will be referred to
as vortex elements. Vortex filaments are here used as elements that
represents the vorticity field. A vortex filaments is delimited by two
points and hence assumes a direction formed by these two points. A
vorticity tube is oriented along the unit vector :math:`\vec{e}_x` of
cross section :math:`dS` and length :math:`l`. It can then be
approximated by a vortex filament of length :math:`l` oriented along the
same direction. The total vorticity of the tube and the vortex filaments
are the same and related by:

.. math::
   \begin{aligned}
       \vec{\omega} \,  dS  = \vec{\Gamma}
   \end{aligned}
   :label: OmegaGamma

where :math:`\vec{\Gamma}` is the circulation intensity of the vortex
filament. If the vorticity tubes are complex and occupy a large volumes,
the projection onto vortex filaments is difficult and the projection
onto vortex particle is more adapted. Yet, assuming the wake is confined
to a thin vorticity layer which defines a velocity jump of know
direction, it is possible to approximate the wake vorticity sheet as a
mesh of vortex filaments. This is the basis of vortex filament wake
methods. Vortex filaments are a singular representation of the vorticity
field, since they occupy a line instead of a volume. To better represent
the vorticity field, the filaments are “inflated”, a process referred to
as regularization (see :numref:`sec:Regularization`). The
regularization of the vorticity field also regularizes the velocity
field and avoids the singularities that would otherwise occur.

.. _sec:vortconv:

Vortex Convection
~~~~~~~~~~~~~~~~~

The governing equation of motion for a vortex filament is given by the
convection equation of a Lagrangian marker:

.. math::
   \frac{d\vec{r}}{dt}=\vec{V}(\vec{r},t)
   :label: VortFilCart

where :math:`\vec{r}` is the position of a Lagrangian marker. The
Lagrangian convection of the filaments stretches the filaments and thus
automatically accounts for the strain part of the vorticity equation.

.. _sec:vortconvPolar:

Vortex Convection in Polar Coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The governing equation of motion for a vortex filament is given by:

.. math:: 
   \frac{d\vec{r}(\psi,\zeta)}{dt}=\vec{V}[\vec{r}(\psi,\zeta),t]
   :label: VortFil

Using the chain rule, Eq. :eq:`VortFil` is rewritten as:

.. math::
   \frac{\partial\vec{r}(\psi,\zeta)}{\partial\psi}+\frac{\partial\vec{r}(\psi,\zeta)}{\partial\zeta}=\frac{\vec{V}[\vec{r}(\psi,\zeta),t]}{\Omega}
   :label: VortFil_expanded

where :math:`d\psi/dt=\Omega` and
:math:`d\psi=d\zeta` (:cite:`Leishman02_1`). Here,
:math:`\vec{r}(\psi,\zeta)` is the position vector of a Lagrangian
marker, and :math:`\vec{V}[\vec{r}(\psi,\zeta)]` is the velocity.

At present, first-order forward Euler method is used to numerically
solve the left-hand side of
Eq. :eq:`VortFil_expanded` for the vortex-filament
location [**IntMethod=5**]. This is an explicit method solved using
Eq. :eq:`Euler`.

.. math::
   \vec{r}(\psi+\Delta\psi_i,\zeta+\Delta\zeta)  = \vec{r}(\psi,\zeta) + \vec{V}(\psi,\zeta) \Delta t
   :label: Euler

Induced Velocity and Velocity Field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The velocity term on the right-hand side of
Eq. :eq:`VortFilCart` is a nonlinear function of the
vortex position, representing a combination of the freestream and
induced velocities (:cite:`Hansen08_1`). The induced
velocities at point :math:`\vec{x}`, caused by each straight-line
filament, are computed using the Biot-Savart law, which considers the
locations of the Lagrangian markers and the intensity of the vortex
elements (:cite:`Leishman02_1`):

.. math::
   d\vec{v}(\vec{x})=\frac{\Gamma}{4\pi}\frac{d\vec{l}\times\vec{r}}{r^3}
   :label: BiotSavart

Here, :math:`\Gamma` is the circulation strength of the filament,
:math:`\vec{dl}` is an elementary length along the filament, and
:math:`\vec{r}` is the vector between a point on the filament and the
control point :math:`\vec{x}`, and :math:`r=|\vec{r}|` is the norm of
the vector. The integration of the Biot-Savart law along the filament
length, delimited by the points :math:`\vec{x}_1` and :math:`\vec{x}_2`
leads to:

.. math::
   \begin{aligned}
     \vec{v}(\vec{x}) 
     =  F_\nu \frac{\Gamma}{4\pi} \frac{(r_1+r_2)}{r_1r_2(r_1r_2+\vec{r}_1\cdot\vec{r}_2)  }\vec{r}_1\times\vec{r}_2
   \end{aligned}
   :label: eq:BiotSavartSegment

with :math:`\vec{r}_1= \vec{x}-\vec{x}_1` and
:math:`\vec{r}_2= \vec{x}-\vec{x}_2`. The factor :math:`F_\nu` is a
regularization parameter that will be discussed in
:numref:`sec:RegularizationFunction`. The filament length
is noted :math:`r_0`, where :math:`\vec{r}_0= \vec{x}_2-\vec{x}_1`. The
distance orthogonal to the filament is:

.. math::
   \begin{aligned}
      \rho = \frac{|\vec{r}_1\times\vec{r}_2|}{r_0}
   \end{aligned}

The velocity at any point of the domain is obtained by superposition of
the velocity induced by all vortex filaments, and by superposition of
the main flow, :math:`\vec{V}_0`, (here assumed divergence free):

.. math::
   \begin{aligned}
    \vec{V}(\vec{x}) = \vec{V}_0 +  \sum_{k} \vec{v}_k(\vec{x})
   \end{aligned}

where the sum is over all the vortex filaments, each of intensity
:math:`\Gamma_k`. The intensity of each filament is determined by
spanwise and time changes of the bound circulation, as discussed in
:numref:`sec:circ`.

.. _sec:Regularization:

Regularization
~~~~~~~~~~~~~~

Regularization and viscous diffusion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The singularity that occurs in Eq. :eq:`BiotSavart`
greatly affects the numerical accuracy of vortex methods. By
regularizing the “1-over-r” kernel of the Biot-Savart law, it is
possible to obtain a numerical method that converges to the
Navier-Stokes equations. The regularization is used to improve the
regularity of the discrete vorticity field, as compared to the “true”
continuous vorticity field. This regularization is usually obtained by
convolution with a smooth function. In this case, the regularization of
the vorticity field and the velocity field are the same. Some
engineering models also perform regularization by directly introducing
additional terms in the denominator of the Biot-Savart velocity kernel.
The factor, :math:`F_\nu`, was introduced in
Eq. :eq:`eq:BiotSavartSegment` to account for
this regularization.

In the convergence proofs of vortex methods, regularization and viscous
diffusion are two distinct aspects. It is yet common practice in vortex
filament methods to blur the notion of regularization with the notion of
viscous diffusion. Indeed, for a physical vortex filament, viscous
effects prevent the singularity from occurring and diffuse the vortex
strength with time. The circular zone where the velocity drops to zero
around the vortex is referred to as the vortex core. An increase of
length of the vortex segment will result in a decrease of the vortex
core radius, and conversely for a decrease of length. Diffusion, on the
other hand, continually spreads the vortex radially.

Because of the previously mentioned analogy, practitioners of vortex
filament methods often refer to regularization as “viscous-core” models
and regularization parameters as “core-radii.” Additionally, viscous
diffusion is often introduced by modifying the regularization parameter
in space and time instead of solving the diffusion from the vorticity
equation. The distinction is made explicit in this document when
clarification is required, but a loose terminology is used when the
context is clear enough.

Determination of the regularization parameter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The regularization parameter is both a function of the physics being
modelled (blade boundary layer and wake) and the choice of
discretization. Contributing factors are the chord length, the boundary
layer height, and the volume that each vortex filament is approximating.
Currently the choice is left to the user (**RegDetMethod\ =0**).
Empirical results for a rotating blade are found in the work of
Gupta (:cite:`Gupta06_1`). As a guideline, the
regularization parameter may be chosen as twice the average spanwise
discretization of the blade. The current implementation will implement
this guideline when the user chooses **RegDetMethod\ =1**. Further
refinement of this option will be considered in the future.

.. _sec:RegularizationFunction:

Regularization functions implemented
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Several regularization functions have been
developed (:cite:`Rankine58_1,Scully75_1,Vatistas91_1`).
At present, five options are available: (1) No correction, (2) the
Rankine method, (3) the Lamb-Oseen method, (4) the Vatistas method, or
(5) the denominator offset method. If no correction method is used,
[**RegFunction=0**], :math:`F_\nu=1`. The remaining methods are detailed
in the following sections. The regularization parameter
(**WakeRegParam**) is noted :math:`r_c` and the distance to the filament
is written :math:`\rho`. The different functions are compared in
:numref:`FilamentRegularization`.

Rankine
'''''''

The Rankine method (:cite:`Rankine58_1`) is the simplest
regularization model. With this method, the Rankine vortex has a finite
core with a solid body rotation near the vortex center and a potential
vortex away from the center. If this method is used,
[**RegFunction=1**], the viscous core correction is given by
Eq. :eq:`rankine`.

.. math::
       F_\nu= \begin{cases} \rho^2/r_c^2 & 0 < \rho < 1 \\
       1 & \rho > 1 \end{cases}
   :label: rankine

Here, :math:`r_c` is the viscous core radius of a vortex filament,
detailed in :numref:`sec:corerad`.

Lamb-Oseen
''''''''''

If this method is used, [**RegFunction=2**], the viscous core correction
is given by Eq. :eq:`lamboseen`.

.. math::
   F_\nu= \bigg[1-\text{exp}(-\frac{\rho^2}{r_c^2})\bigg]
   :label: lamboseen

Vatistas
''''''''

If this method is used, [**RegFunction=3**], the viscous core correction
is given by Eq. :eq:`vatistas`.

.. math::
   F_\nu
   = \frac{\rho^2}{(\rho^{2n}+r_c^{2n})^{1/n}}
   = \frac{(\rho/r_c)^2}{(1 + (\rho/r_c)^{2n})^{1/n}}
   :label: vatistas

Here, :math:`\rho` is the distance from a vortex segment to an arbitrary
point (:cite:`Abedi16_1`). Research from rotorcraft
applications suggests a value of :math:`n=2`, which is used in this
work (:cite:`Bagai93_1`).

Denominator offset/cut-off
''''''''''''''''''''''''''

If this method is used, [**RegFunction=4**], the viscous core correction
is given by Eq. :eq:`denom`

.. math::
   \begin{aligned}
     \vec{v}(\vec{x}) 
     =   \frac{\Gamma}{4\pi} \frac{(r_1+r_2)}{r_1r_2(r_1r_2+\vec{r}_1\cdot\vec{r}_2) + r_c^2  r_0^2} \vec{r}_1\times\vec{r}_2
   \end{aligned}
   :label: denom

Here, the singularity is removed by introducing an additive factor in
the denominator of
Eq. :eq:`eq:BiotSavartSegment`, proportional to
the filament length :math:`r_0`. In this case, :math:`F_\nu=1`. The
method is found in the work of van Garrel
(:cite:`Garrel03_1`).

.. _sec:corerad:

Time Evolution of the Regularization Parameter–Core Spreading Method
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are four available methods by which the regularization parameter
may evolve with time: (1) constant value, (2) stretching, (3) wake age,
or (4) stretching and wake age. The three latter methods blend the
notion of viscous diffusion with the notion of regularization. The
notation :math:`r_{c0}` used in this section corresponds to input file
parameter value .

Constant
''''''''

If a constant value is selected, [**WakeRegMethod=0**], the value of
:math:`r_c` remains unchanged for all Lagrangian markers throughout the
simulation and taken as the value given with the parameter in meters.

.. math::
   r_c(\zeta) = r_{c0}
   :label: cst

Here, :math:`\zeta` is the vortex wake age, measured from its emission
time.

Stretching
''''''''''

If the stretching method is selected, [**WakeRegMethod=1**], the viscous
core radius is modeled by Eq. :eq:`stretch`.

.. math::
   r_c(\zeta,\epsilon) = \sqrt{r_{c0}^2+\int_0^\zeta(1+\epsilon)^{-1}d\zeta}
   :label: stretch

.. math::
   \epsilon = \frac{\Delta l}{l}

Here, :math:`\epsilon` is the vortex-filament strain, and :math:`l` is
the filament length, and :math:`\Delta l` is the change of length
between two time steps. The integral in Eq. :eq:`stretch`
represents strain effects.

Wake Age / Core-Spreading
'''''''''''''''''''''''''

If the wake age method is selected, [], the viscous core radius is
modeled by Eq. :eq:`age`.

.. math::
   r_c(\zeta) = \sqrt{r_{c0}^2+4\alpha\delta\nu \zeta}
   :label: age

where :math:`\alpha=1.25643`, :math:`\nu` is kinematic viscosity, and
:math:`\delta` is a viscous diffusion parameter (typically between
:math:`1` and :math:`1,000`). The parameter :math:`\delta` is provided
in the input file as **CoreSpreadEddyVisc**. Here, the term,
:math:`4\alpha\delta\nu \zeta`, accounts for viscous effects as the wake
propagates downstream. The higher the background turbulence, the more
diffusion of the vorticity with time, and the higher the value of
:math:`\delta` should be. The method is often referred to as the
core-spreading method. It is a way to partially account for viscous
diffusion of the vorticity without solving for the interaction between
the wake vorticity or between the vorticity from the wake and the background
flow. Setting **DiffusionMethod==1** is the same as using the wake age method,
[**WakeRegMethod=2**].

Stretching and Wake Age
'''''''''''''''''''''''

If the stretching and wake-age method is selected [**WakeRegMethod=3**],
the viscous core radius is modeled by
Eq. :eq:`stretchandage`.

.. math::
   r_c(\zeta,\epsilon) = \sqrt{r_{c0}^2 + 4\alpha\delta\nu \zeta + \int_0^\zeta(1+\epsilon)^{-1}d\zeta}
   :label: stretchandage

.. _sec:diffusion:

Diffusion
~~~~~~~~~

The viscous-splitting assumption is used to solve for the convection and
diffusion of the vorticity separately. The diffusion term
:math:`\nu \Delta \vec{\omega}` represents molecular diffusion. This
term allows for viscous connection of vorticity lines. Also, turbulent
flows will diffuse the vorticity in a similar manner based on a
turbulent eddy viscosity.

The parameter **ViscousDiffusion** is used to switch between viscous diffusion
methods.  Currently, only the core-spreading method is implemented. The method
is described in :numref:`sec:corerad` since it is equivalent to the increase of
the regularization parameter with the wake age.

.. _sec:circ:

Lifting-Line Circulation
~~~~~~~~~~~~~~~~~~~~~~~~

The code relies on a lifting-line formulation. Lifting-line methods
effectively lump the loads at each cross-section of the blade onto the
mean line of the blade and do not account directly for the geometry of
each cross-section. In the vorticity-based version of the lifting-line
method, the blade is represented by a line of varying circulation. The
line follows the motion of the blade and is referred to as “bound”
circulation. The bound circulation does not follow the same dynamic
equation as the free vorticity of the wake. Instead, the intensity is
linked to airfoil lift via the Kutta-Joukowski theorem. Spanwise
variation of the bound circulation results in vorticity being emitted
into the the wake, referred to as “trailed vorticity”. Time changes of
the bound circulation are also emitted in the wake, referred to as “shed”
vorticity. Three methods are implemented to determine the bound circulation
strength. They are selected using the input **CircSolvMethod**, and are
presented in the subsequent paragraphs. At the end of a time step, the
circulation of each vortex element is propagated downstream so that vortex
elements with a new intensity can be emitted from the blade at the next time
step.

Cl-Based Iterative Method
^^^^^^^^^^^^^^^^^^^^^^^^^

The Cl-based iterative method is extensively described in the work from van
Garrel and it is only briefly presented here (:cite:`Garrel03_1`). The method
was implemented following the same approach and notations as van Garrel. At
present, it is the preferred method to compute the circulation along the blade
span. It is selected with **CircSolvMethod==1**. In this method, the blade is
discretized into a finite number of segments placed along the lifting line
(i.e., the blade aerodynamic center line), representing the bound circulation,
:math:`\Gamma_b`. The circulation is solved within a nonlinear iterative solver
that makes use of the polar data at each control point located on the lifting
line.

No-flow-through method
^^^^^^^^^^^^^^^^^^^^^^

A Weissinger-L-based representation (:cite:`Weissinger47_1`)
of the lifting surface is also
available (:cite:`Bagai94_1,Gupta06_1,Ribera07_1`). In this
method, the circulation is solved by satisfying a no-flow through
condition at the 3/4-chord points.

Prescribed circulation
^^^^^^^^^^^^^^^^^^^^^^

The final available method prescribes a constant circulation. A user
specified spanwise distribution of circulation is prescribed onto the
blades.

State-Space Representation and Integration with OpenFAST
--------------------------------------------------------

The OLAF module has been integrated into the latest version of OpenFAST
via *AeroDyn15*, following the OpenFAST modularization
framework (:cite:`Jonkman13_1,Sprague15_1`). To follow the
OpenFAST framework, the vortex code is written as a module, and its
formulation comprises state, constraint, and output equations. The data
manipulated by the module include the following vectors: inputs,
:math:`\vec{u}`; states, :math:`\vec{x}`; constrained state,
:math:`\vec{z}`; outputs, :math:`\vec{y}`; and constant parameters,
:math:`\vec{p}`. The vectors are defined as follows:

-  Inputs, :math:`\vec{u}~-` a set of values supplied to the module
   that, along with the states, are needed to calculate future states
   and the system’s output.

-  Outputs, :math:`\vec{y}~-` a set of values calculated and returned
   by the module that depend on the states, inputs, and/or parameters
   through output equations.

-  States, :math:`\vec{x}~-` a set of internal values of the module
   that are influenced by the inputs and used to calculate future state
   values and the output. Continuous states are employed, meaning that
   the states are differentiable in time and characterized by continuous
   time-differential equations.

-  Constraint states, :math:`\vec{z}~-` algebraic variables that are
   calculated using a nonlinear solver, based on values from the current
   time step.

-  Parameters, :math:`\vec{p}~-` a set of internal system values that
   are independent of the states and inputs. The parameters can be fully
   defined at initialization and characterize the system’s state
   equations and output equations.

The parameters of the vortex code include:

-  Fluid characteristics: kinematic viscosity, :math:`\nu`

-  Airfoil characteristics: polar data: (:math:`C_l(\alpha)`,
   :math:`C_d(\alpha)`, :math:`C_m(\alpha)`), and chord :math:`c`

-  Algorithmic methods and parameters for regularization, viscous
   diffusion, discretization, wake geometry, acceleration, and so on.

The inputs of the vortex code are:

-  Position, orientation, translational velocity, and rotational
   velocity of the different nodes of the lifting lines
   (:math:`\vec{r}_{ll}`, :math:`\Lambda_{ll}`,
   :math:`\vec{\dot{r}}_{ll}`, and :math:`\vec{\omega}_{ll}`,
   respectively), gathered into the vector,
   :math:`\vec{x}_{\text{elast},ll}`, for conciseness. These quantities
   are handled using the mesh-mapping functionality and data structure
   of OpenFAST.

-  Undisturbed velocity field at requested locations (lifting-line
   points, :math:`\vec{r}_{ll}`, and a set of locations requested by the
   vortex code, :math:`\vec{r}_r`), written
   :math:`\vec{v}_0=[\vec{v}_{0,ll}, \vec{v}_{0,r}]`. Based on the
   parameters, this velocity field may contain the following influences:
   freestream, shear, veer, turbulence, tower, and nacelle disturbance.
   The locations where the velocity field is requested are typically the
   location of the Lagrangian markers.

The constraint states are:

-  The circulation intensity along the lifting lines,
   :math:`\Gamma_{ll}`.

The continuous states are:

-  The position of the Lagrangian markers, :math:`\vec{r}_m`

-  The vorticity associated with each vortex element,
   :math:`\vec{\omega}_e`. For a projection of the vorticity onto vortex
   segments, this corresponds to the circulation,
   :math:`\vec{\Gamma}_e`, where for each segment,
   :math:`\vec{\Gamma}_e= \Gamma_e \vec{dl}_e =\vec{\omega}_e dV_e`,
   with :math:`\vec{dl}_e` and :math:`dV_e`, the vortex segment length
   and its equivalent vortex volume.

The outputs are  [1]_:

-  The induced velocity at the lifting-line nodes,
   :math:`\vec{v}_{i,ll}`

-  The locations where the undisturbed wind needs to be computed,
   :math:`\vec{r}_{r}` (typically :math:`\vec{r_{r}}=\vec{r}_m`).

State, Constraint, and Output Equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An overview of the main states, constraints, and output equations is
given in this paragraph. More details are provided in
:numref:`sec:FVW`. The constraint equation is used to determine
the circulation distribution along the span of each lifting line. For
the van Garrel method, this circulation is a function of the angle of
attack along the blade and the airfoil coefficients. The angle of attack
at a given lifting-line node is a function of the undisturbed velocity,
:math:`\vec{v}_{0,ll}`, and the velocity induced by the vorticity,
:math:`\vec{v}_{i,ll}`, at that point. Part of the induced velocity is
caused by the vorticity being shed and trailed at the current time step,
which in turn is a function of the circulation distribution along the
lifting line. This constraint equation may be written as:

.. math::
   \vec{Z} = \vec{0} = \vec{\Gamma}_{ll} - \vec{\Gamma}_p(\vec{\alpha}(\vec{x},\vec{u}),\vec{p})

where :math:`\vec{\Gamma}_p` is the function that returns the
circulation along the blade span, based on the distribution of angle of
attacks and the airfoil characteristics. In practice, this nonlinear
equation is solved using an iterative algorithm. The state equation
specifies the time evolution of the vorticity and the convection of the
Lagrangian markers:

.. math::
   \begin{aligned}
       \frac{d \vec{\omega}_e}{dt} &= \left[(\vec{\omega}\cdot\nabla)\vec{v} + \nu\nabla^2 \vec{\omega} \right]_e
           \\
       \frac{d \vec{r}_m}{dt} &= \vec{V}(\vec{r}_m)
    =\vec{V}_0(\vec{r}_m)  + \vec{v}_\omega(\vec{r}_m)
    =\vec{V}_0(\vec{r}_m)  + \vec{V}_\omega(\vec{r}_m, \vec{r}_m, \vec{\omega})
   \end{aligned}
   :label: eq:Convection

where :math:`\vec{v}_\omega` is the velocity induced by the vorticity in
the domain; :math:`\vec{V}_\omega(\vec{r},\vec{r}_m,\vec{\omega})` is
the function that computes this induced velocity at a given point,
:math:`\vec{r}`, based on the location of the Lagrangian markers and the
intensity of the vortex elements; and the subscript, :math:`e`,
indicates that a quantity is applied to an element. The vorticity,
:math:`\vec{\omega}`, is recovered from the vorticity of the vortex
elements by means of discrete convolutions. For vortex-segment
simulations, the viscous-splitting algorithm is used, and the convection
step (Eq. :eq:`eq:Convection`) is the main state
equation being solved for. The vorticity stretching is automatically
accounted for, and the diffusion is performed *a posteriori*. The
velocity function, :math:`\vec{V}_\omega`, uses the Biot-Savart law. The
output equation is:

.. math::
   \begin{aligned}
      \vec{y}_1&=\vec{v}_{i,ll} = \vec{V}_\omega ( \vec{r}_{ll}, \vec{r}_m, \vec{\omega})= \\
      \vec{y}_2&=\vec{r}_{r}
   \end{aligned}

Integration with AeroDyn15
~~~~~~~~~~~~~~~~~~~~~~~~~~

The vortex code has been integrated as a submodule of the aerodynamic
module of OpenFAST, *AeroDyn15*. The data workflow between the different
modules and submodules of OpenFAST is illustrated in
Figure :ref:`FAST-FVW`.

This integration required a restructuring of the *AeroDyn15* module to
isolate the parts of the code related to tower shadow modeling,
induction computation, lifting-line-forces computations, and dynamic
stall. The dynamic stall model will be adapted when used in conjunction
with the vortex code to ensure the effect of shed vorticity is not
accounted for twice. The interface between *AeroDyn15* and the inflow
module, *InflowWind*, was accommodated to include the additionally
requested points by the vortex code.

.. [1]
   The loads on the lifting line are not an output of the vortex code;
   their calculation is handled by a separate submodule of *AeroDyn*.
