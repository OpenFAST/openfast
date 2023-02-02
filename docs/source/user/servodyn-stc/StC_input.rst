.. _StC-Input:

Input Files
===========

The user configures each StC instance with a separate input file. This input
file defines the location of the StC relative to its mounting location, and
defines the properties.  It can also be used with an external forces file to
apply a timeseries load at the location (primarily used for diagnostic
purposes).
 

Units
-----

Structural Control uses the SI system (kg, m, s, N). Angles are assumed to be in
radians unless otherwise specified.


.. raw:: html

   <hr>



.. _StC-Locations:

Structural Control Locations 
----------------------------

The Structural Control input file defines the location and properties of the StC
instance.  The location is relative to the type of StC given in the main
ServoDyn input file (see :numref:`SrvD-StC-Inputs`).  The four location types
are Nacelle, Tower, Blade, and Platform.

The mapping information for the StC will be given in the main OpenFAST summary
file.


Nacelle StC
~~~~~~~~~~~

This StC mounting location is attached relative to the nacelle reference point.
It will track with all nacelle motions (including motions due to yaw, tower
flex, and platform motions).


Tower StC
~~~~~~~~~

This StC mounting location is attached to the tower mesh at the height specified
above the tower base.  This StC attachment will move with the line mesh at that
height. For example, an StC mounted at 85 m on a 90 m tower will move with the
mesh line corresponding to the 85 m height position on the tower center line.


Blade StC
~~~~~~~~~

This StC mounting location is attached to the blade structural center at the
specified distance from the blade root along the z-axis of the blade (IEC
blade coordinate system).  This location will follow all blade deformations and
motions (including blade twist when used with BeamDyn).  This option is
available with both the BeamDyn and ElastoDyn blade representations.

When this option is used, identical StCs will be attached at each of the blades.
The response if each blade mounted StC is tracked separately and is available in
the output channels given in the ServoDyn tab of the
:download:`OutListParameters.xlsx <../../../OtherSupporting/OutListParameters.xlsx>`.


Platform StC
~~~~~~~~~~~~

This StC mounting location is located relative to the platform reference point.
When a rigid body platform is modeled (such as a rigid semi-submersible), it is
attached to the platform reference point.  When a flexible floating body is
modeled, the StC is attached to the SubDyn mesh.


.. raw:: html

   <hr>

.. _StC-Input-File:

Structural Control Input File
-----------------------------

The input file may have an arbitrary number of commented header lines, and
commented lines at any location in the input file.
:download:`(Example Structural Control input file for tuned mass damper on
tower for NREL 5 MW TLP) <ExampleFiles/NRELOffshrBsline5MW_StC.dat>`:

Simulation Control
~~~~~~~~~~~~~~~~~~

**Echo** [flag]

   Echo input data to <RootName>.ech  


StC Degrees of Freedom
----------------------

**StC_DOF_MODE** [switch]

   DOF mode   {0: No StC or TLCD DOF; 1: StC_X_DOF, StC_Y_DOF, and/or StC_Z_DOF
   (three independent StC DOFs); 2: StC_XY_DOF (Omni-Directional StC); 3: TLCD;
   4: Prescribed force/moment time series; 5: Force determined by external DLL}


**StC_X_DOF** [flag]

   DOF on or off for StC X   *[Used only when* **StC_DOF_MODE==1** *]*

**StC_Y_DOF** [flag]

   DOF on or off for StC Y   *[Used only when* **StC_DOF_MODE==1** *]*

**StC_Z_DOF** [flag]

   DOF on or off for StC Z   *[Used only when* **StC_DOF_MODE==1** *]*


StC Location
------------

The location of the StC is relative to the component it is attached to.  This is
specified in the main ServoDyn input file.  See description above.

**StC_P_X** [m]

   At rest X position of StC  

**StC_P_Y** [m]

   At rest Y position of StC  

**StC_P_Z** [m]

   At rest Z position of StC  


StC Initial Conditions
----------------------

*used only when* **StC_DOF_MODE==1 or 2**

**StC_X_DSP** [m]

   StC X initial displacement   *[relative to at rest position]*

**StC_Y_DSP** [m]

   StC Y initial displacement   *[relative to at rest position]*

**StC_Z_DSP** [m]

   StC Z initial displacement   *[relative to at rest position; used only when*
   **StC_DOF_MODE==1** *and* **StC_Z_DOF==TRUE** *]*

**StC_Z_PreLd** [N]

   StC Z spring preload. Either a direct value for the spring preload in
   Newtons,  or **"gravity"** for pre-loading spring to shift the at rest
   position of the StC Z mass when gravity is acting on it using
   :math:`F_{Z_{PreLoad}} = M_Z * G`, or **"none"** to disable spring pre-load.
   See :numref:`SrvD-StCz-PreLoad` for details of implementation.
   *[used only when* **StC_DOF_MODE=1** and **StC_Z_DOF=TRUE** *]*


StC Configuration
-----------------

*used only when* **StC_DOF_MODE==1 or 2**

**StC_X_PSP** [m]

   Positive stop position  -- maximum X mass displacement

**StC_X_NSP** [m]

   Negative stop position  -- minimum X mass displacement

**StC_Y_PSP** [m]

   Positive stop position  -- maximum Y mass displacement

**StC_Y_NSP** [m]

   Negative stop position  -- minimum Y mass displacement

**StC_Z_PSP** [m]

   Positive stop position  -- maximum Z mass displacement *[used only when*
   **StC_DOF_MODE==1** *and* **StC_Z_DOF==TRUE** *]*

**StC_Z_NSP** [m]

   Negative stop position -- minimum Z mass displacement *[used only when*
   **StC_DOF_MODE==1** *and* **StC_Z_DOF==TRUE** *]*

StC Mass, Stiffness, & Damping
------------------------------

*used only when* **StC_DOF_MODE==1 or 2**

**StC_X_M** [kg]

   StC X mass   *[used only when* **StC_DOF_MODE==1** *and* **StC_X_DOF==TRUE**
   *]*

**StC_Y_M** [kg]

   StC Y mass   *[used only when* **StC_DOF_MODE==1** *and* **StC_Y_DOF==TRUE**
   *]*

**StC_Z_M** [kg]

   StC Z mass   *[used only when* **StC_DOF_MODE==1** *and* **StC_Z_DOF==TRUE**
   *]*

**StC_XY_M** [kg]

   StC XY mass   *[used only when* **StC_DOF_MODE==2** *]*

**StC_X_K** [N/m]

   StC X stiffness  

**StC_Y_K** [N/m]

   StC Y stiffness  

**StC_Z_K** [N/m]

   StC Z stiffness   *[used only when* **StC_DOF_MODE==1** *and*
   **StC_Z_DOF==TRUE** *]*

**StC_X_C** [N/(m/s)]

   StC X damping  

**StC_Y_C** [N/(m/s)]

   StC Y damping  

**StC_Z_C** [N/(m/s)]

   StC Z damping   *[used only when* **StC_DOF_MODE==1** *and*
   **StC_Z_DOF==TRUE** *]*

**StC_X_KS** [N/m]

   Stop spring X stiffness  

**StC_Y_KS** [N/m]

   Stop spring Y stiffness  

**StC_Z_KS** [N/m]

   Stop spring Z stiffness   *[used only when* **StC_DOF_MODE==1** *and
   StC_Z_DOF==TRUE]*

**StC_X_CS** [N/(m/s)]

   Stop spring X damping  

**StC_Y_CS** [N/(m/s)]

   Stop spring Y damping  

**StC_Z_CS** [N/(m/s)]

   Stop spring Z damping   *[used only when* **StC_DOF_MODE==1** *and*
   **StC_Z_DOF==TRUE** *]*


StC User-Defined Spring Forces
------------------------------

*used only when* **StC_DOF_MODE==1 or 2**

**Use_F_TBL** [flag]

   Use spring force from user-defined table  

**NKInpSt** [-]

   Number of spring force input stations

The table is expected to contain 6 columns for displacements and equvalent
sprint forces: **X**, **F_X**, **Y**, **F_Y**, **Z**, and **F_Z**.
Displacements are in meters (m) and forces in Newtons (N).

Example spring forces table:

.. container::
   :name: Tab:SpringForce

   .. literalinclude:: ExampleFiles/SpringForce.txt
      :language: none


StructCtrl Control
------------------
*used only when* **StC_DOF_MODE==1 or 2**

**StC_CMODE** [switch]

   Control mode   {0:none; 1: Semi-Active Control Mode; 2: Active Control Mode}.
   When using StC_DOF_MODE==5, StC_CMODE must be 2.

**StC_SA_MODE** [-]

   Semi-Active control mode {1: velocity-based ground hook control; 2: Inverse
   velocity-based ground hook control; 3: displacement-based ground hook control
   4: Phase difference Algorithm with Friction Force 5: Phase difference
   Algorithm with Damping Force}

**StC_X_C_HIGH** [-]

   StC X high damping for ground hook control

**StC_X_C_LOW** [-]

   StC X low damping for ground hook control

**StC_Y_C_HIGH** [-]

   StC Y high damping for ground hook control

**StC_Y_C_LOW** [-]

   StC Y low damping for ground hook control

**StC_Z_C_HIGH** [-]

   StC Z high damping for ground hook control *[used only when*
   **StC_DOF_MODE==1** *and* **StC_Z_DOF==TRUE** *]*

**StC_Z_C_LOW** [-]

   StC Z low damping for ground hook control  *[used only when*
   **StC_DOF_MODE==1** *and* **StC_Z_DOF==TRUE** *]*

**StC_X_C_BRAKE** [-]

   StC X high damping for braking the StC *[currently unused.  set to zero]*

**StC_Y_C_BRAKE** [-]

   StC Y high damping for braking the StC *[currently unused.  set to zero]*

**StC_Z_C_BRAKE** [-]

   StC Z high damping for braking the StC *[used only when* **StC_DOF_MODE==1**
   *and* **StC_Z_DOF==TRUE** *]* *[currently unused.  set to zero]*



TLCD -- Tuned Liquid Column Damper
----------------------------------

*used only when* **StC_DOF_MODE==3**

**L_X** [m]

   X TLCD total length

**B_X** [m]

   X TLCD horizontal length

**area_X** [m^2]

   X TLCD cross-sectional area of vertical column

**area_ratio_X** [-]

   X TLCD cross-sectional area ratio  *[vertical column area divided by
   horizontal column area]*

**headLossCoeff_X** [-]

   X TLCD head loss coeff

**rho_X** [kg/m^3]

   X TLCD liquid density

**L_Y** [m]

   Y TLCD total length  

**B_Y** [m]

   Y TLCD horizontal length  

**area_Y**        [m^2]

   Y TLCD cross-sectional area of vertical column

**area_ratio_Y**  [-] 

   Y TLCD cross-sectional area ratio *[vertical column area divided by
   horizontal column area]*

**headLossCoeff_Y** [-]

   Y TLCD head loss coeff

**rho_Y** [kg/m^3]

   Y TLCD liquid density

Prescribed Time Series
----------------------

A prescribed time series of forces and moments may be applied in place of the
StC damper.  The force and moment may be applied either in a global coordinate
frame, or in a local (following) coordinate frame.  This feature is *used only
when* **StC_DOF_MODE==4**.

**PrescribedForcesCoord** [switch]

   Prescribed forces are in global or local coordinates  {1: global; 2: local}.  
   When using StC_DOF_MODE==5, PrescribedForcesCoord must be 1.

**PrescribedForcesFile** [-]

   Filename for the prescribed forces.  The format expected is 7 columns: time,
   FX, FY, FZ, MX, MY, MZ.  Values will be interpolated from the file between
   the given timestep and value sets.  The input file may have an arbitrary
   number of commented header lines, and commented lines at any location in the
   input file.

Example prescribed time series file :download:`(example prescribed force
timeseries) <ExampleFiles/PrescribedForce.txt>`:

.. container::
   :name: Tab:PrescribedForce

   .. literalinclude:: ExampleFiles/PrescribedForce.txt
      :language: none


