.. _OLAF-Input-Files:

Input Files
===========

No lines should be added or removed from the input files, except in tables where
the number of rows is specified.

Units
-----

OLAF uses the International System of Units (e.g., kg, m, s, N). Angles are
assumed to be in degrees unless otherwise specified.

OLAF Primary Input File
-----------------------

The primary OLAF input file defines general free wake options, circulation model
selection and specification, near- and far-wake length, and wake visualization
options. Each section within the file corresponds to an aspect of the OLAF model. For most parameters, the user may
specify the value "default" (with or without quotes), in which case a default
value, defined below, is used by the program.

See :numref:`OLAF-Primary-Input-File` for a sample OLAF primary input file.

General Options
~~~~~~~~~~~~~~~

**IntMethod** [switch] specifies which integration method will be used to
convect the Lagrangian markers. There are four options: 1) fourth-order
Runge-Kutta *[1]*, 2) fourth-order Adams-Bashforth *[2]*, 3) fourth-order
Adams-Bashforth-Moulton *[3]*, and 4) first-order forward Euler *[5]*. The
default option is *[5]*. These methods are specified in :numref:`sec:vortconv`.

**DTfvw** [sec] specifies the time interval at which the module will update the
wake. The time interval must be a multiple of the time step used by
*AeroDyn15*. The blade circulation is updated at each intermediate time
step based on the intermediate blades positions and wind velocities. The default
value is :math:`dt_{aero}`, where :math:`dt_{aero}` is the time step used by
AeroDyn.

**FreeWakeStart** [sec] specifies at what time the wake evolution is classified
as “free." Before this point is reached, the Lagrangian markers are simply
convected with the freestream velocity. After this point, induced velocities are
computed and affect the marker convection. If a time less than or equal to zero
is given, the wake is “free" from the beginning of the simulation.  The default
value is :math:`0`.

**FullCircStart** [sec] specifies at what time the blade circulation reaches its
full strength. If this value is specified to be :math:`>0`, the circulation is
multiplied by a factor of :math:`0` at :math:`t=0` and linearly increasing to a
factor of :math:`1` for :math:`t>\textit{FullCircStart}`. The default
value is :math:`0`.

Circulation Specifications
~~~~~~~~~~~~~~~~~~~~~~~~~~

**CircSolvMethod** [switch] specifies which circulation method is used. There
are three options: 1) :math:`C_l`-based iterative procedure *[1]*, 2) no-flow
through *[2]*, and 3) prescribed *[3]*. The default option is *[1]*. These
methods are described in :numref:`sec:circ`.

**CircSolvConvCrit** [-] specifies the dimensionless convergence criteria used
for solving the circulation. This variable is only used if
*CircSolvMethod* = *[1]*. The default value is
:math:`0.001`, corresponding to :math:`0.1\%` error in the circulation between
two iterations.

**CircSolvRelaxation** [-] specifies the relaxation factor used to solve the
circulation.  This variable is only used if *CircSolvMethod* =
*[1]*. The default value is :math:`0.1`.

**CircSolvMaxIter** [-] specifies the maximum number of iterations used to solve
the circulation. This variable is only used if *CircSolvMethod* = *[1]*. The
default value is :math:`30`.

**PrescribedCircFile** [quoted string] specifies the file containing the
prescribed blade circulation. This option is only used if *CircSolvMethod* =
*[3]*.  The circulation file format is a delimited file with one header line and
two columns. The first column is the dimensionless radial position [r/R]; the
second column is the bound circulation value in [m\ :math:`^2`/s].  The radial
positions do not need to match the AeroDyn node locations. A sample prescribed
circulation file is given in :numref:`Prescribed-Circulation-Input-File`.


Wake Extent and Discretization Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**nNWPanel** [-] specifies the number of FVW time steps (**DTfvw**) for which
the near-wake lattice is computed. In the future, this value will be defined as
an azimuthal span in degrees or a downstream distance in rotor diameter.

**WakeLength** [D] specifies the length, in rotor diameters, of the far wake.
The default value is :math:`8`. [1]_

**FreeWakeLength** [D] specifies the length, in rotor diameters, for which the
turbine wake is convected as “free." If *FreeWakeLength* is greater than
*WakeLength*, then the entire wake is free. Otherwise, the Lagrangian markers
located within the buffer zone delimited by *FreeWakeLength* and *WakeLength*
are convected with the average velocity. The default value is :math:`6`. [2]_

**FWShedVorticity** [flag] specifies whether shed vorticity is included in the
far wake. The default value is *[False]*, specifying that the far wake consists
only of the trailed vorticity from the root and tip vortices.

Wake Regularization and Diffusion Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**DiffusionMethod** [switch] specifies which diffusion method is used to account
for viscous diffusion. There are two options: 1) no diffusion *[0]* and 2) the
core-spreading method *[1]*. The default option is *[0]*.

**RegDetMethod** [switch] specifies which method is used to determine the
regularization parameters. There are two options: 1) manual *[0]* and 2)
optimized *[1]*. The manual option requires the user to specify the parameters
listed in this subsection. The optimized option determines the parameters for
the user.  The default option is *[0]*.

**RegFunction** [switch] specifies the regularization function used to remove
the singularity of the vortex elements, as specified in
:numref:`sec:vortconv`. There are five options: 1) no correction *[0]*,
2) the Rankine method *[1]*, 3) the Lamb-Oseen method *[2]*, 4) the Vatistas
method *[3]*, and 5) the denominator offset method *[4]*. The functions are
given in . The default option is *[3]*.

**WakeRegMethod** [switch] specifies the method of determining viscous core
radius (i.e., the regularization parameter). There are three options: 1)
constant *[1]*, 2) stretching *[2]*, and 3) age *[3]*. The methods are
described in :numref:`sec:corerad`. The default option is *[1]*.

**WakeRegParam** [m] specifies the wake regularization parameter, which is the
regularization value used at the initialization of a vortex element. If the
regularization method is “constant”, this value is used throughout the wake.

**BladeRegParam** [m] specifies the bound vorticity regularization parameter,
which is the regularization value used for the vorticity elements bound to the
blades.

**CoreSpreadEddyVisc** [-] specifies the eddy viscosity parameter
:math:`\delta`.  The parameter is used for the core-spreading method
(*DiffusionMethod* = *[1]*) and the regularization method with age
(*WakeRegMethod* = *[3]*). The variable :math:`\delta` is described in
:numref:`sec:corerad`. The default value is :math:`100`.

Wake Treatment Options
~~~~~~~~~~~~~~~~~~~~~~

**TwrShadowOnWake** [flag] specifies whether the tower potential flow and tower
shadow have an influence on the wake convection. The tower shadow model, when
activated in AeroDyn, always has an influence on the lifting line, hence the
induction and loads on the blade. This option only concerns the wake. The
default option is *[False]*.

**ShearVorticityModel** [switch] specifies whether shear vorticity is modeled in
addition to the sheared inflow prescribed by *InflowWind*. There are two
options: 1) no treatment *[0]* and 2) mirrored vorticity *[1]*. The mirrored
vorticity accounts for the ground effect. Dedicated options to account for the
shear vorticity will be implemented at a later time. The shear velocity profile
is handled by *InflowWind* irrespective of this input. The default option is
*[0]*.


Speedup Options
~~~~~~~~~~~~~~~

**VelocityMethod** [switch] specifies the method used to determine the velocity.
There are two options: 1) Biot-Savart law applied to the vortex segments *[1]*
and 2) tree formulation using a particle representation *[2]*. The default
option is *[1]*.

**TreeBranchFactor** [-] specifies the dimensionless distance, in branch radius,
above which a multipole calculation is used instead of a direct evaluation. This
option is only used in conjunction with the tree code
(*VelocityMethod* = *[2]*).

**PartPerSegment** [-] specifies the number of particles that are used when a
vortex segment is represented by vortex particles. The default value is
:math:`1`.

Output Options
~~~~~~~~~~~~~~

**WrVTK** [flag] specifies if Visualization Toolkit (VTK) visualization files
are to be written out. *WrVTK* = *[0]* does not write out any VTK files. *WrVTK*
= *[1]* outputs a VTK file at every time step. The outputs are written in the
folder, ``vtk_fvw.``   The parameters *WrVTK*, *VTKCoord*, and *VTK_fps* are
independent of the glue code VTK output options.


**VTKBlades** [-] specifies how many blade VTK files are to be written out.
*VTKBlades* :math:`= n` outputs VTK files for :math:`n` blades, with :math:`0`
being an acceptable value. The default value is :math:`1`.

**VTKCoord** [switch] specifies in which coordinate system the VTK files are
written.  There are two options: 1) global coordinate system *[1]* and 2) hub
coordinate system *[2]*. The default option is *[1]*.

**VTK_fps** [:math:`1`/sec] specifies the output frequency of the VTK files. The
provided value is rounded to the nearest allowable multiple of the time step.
The default value is :math:`1/dt_\text{fvw}`. Specifying *VTK_fps* = *[all]*,
is equivalent to using the value :math:`1/dt_\text{aero}`.

AeroDyn15 Input File
--------------------
Input file modifications
~~~~~~~~~~~~~~~~~~~~~~~~

As OLAF is incorporated into the *AeroDyn15* module, a wake computation option
has been added to the *AeroDyn15* input file and a line has been added. These
additions are as follows.

**WakeMod** specifies the type of wake model that is used. *WakeMod* = *[3]* has
been added to allow the user to switch from the traditional BEM method to the
OLAF method.

**FVWFile** [string] specifies the OLAF module file, the path is relative to the
AeroDyn file, unless an absolute path is provided.


Relevant sections
~~~~~~~~~~~~~~~~~
The BEM options (e.g. tip-loss, skew, and dynamic models) are read and discarded
when *WakeMod* = *[3]*. The following sections and parameters remain relevant and
are used by the vortex code:

  - general options (e.g., airfoil and tower modeling);
  - environmental conditions;
  - dynamic stall model options;
  - airfoil and blade information;
  - tower aerodynamics; and
  - outputs.

.. [1]
   At present, this variable is called nFWPanel and specified as the number of far
   wake panels. This will be changed soon.

.. [2]
   At present, this variable is called nFWPanelFree and specified as the number of
   free far wake panels. This will be changed soon.
