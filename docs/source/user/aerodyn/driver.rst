
.. _ad_drivers:

AeroDyn drivers
===============


Two AeroDyn drivers are available, the baseline AeroDyn driver intended to run multiple simulations of horizontal axis wind turbines, and a general purpose driver 

.. _addm_driver-input-file:




Baseline AeroDyn Driver
-----------------------

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


Multi-rotor AeroDyn Driver
--------------------------


.. figure:: figs/MultiRotor.png
   :width: 80%
           
   Definition of multiple rotors. 


- A turbine is assumed to consist of an optional tower, a hub, and multiple blades. 
- Different frames and origins are defined: the turbine frame, the hub frame, and the blade frames. 
- The tower points are defined in the turbine coordinates
- The hub frame and origin is defined with respect to the turbine coordinates
- The blade frames and origins are defined with respectto the hub coordinates
- The blades are rigidly attached to the hub, and rotate around the x axis of the hub frame. 
- For each blade, the blade frame is such that the zb-axis is along the span, the yb axis is directed towards the "trailing edge" in the absence of pitch, and the xb-axis is directed towards the suction side in the absence of pitch.



**Features:**

- Multiple rotors
- Arbitrary number of blades per rotors
- One tower per rotor 
- Arbitrary rigid body motion of each rotor (tower stays fixed for now). For convenience, we'll add simple prescribed motions as well, such as rotations (and maybe 6-DOF oscillations later)
   

**AeroDyn changes:**

- No changes to AeroDyn input file (except repetition of tower properties).

- Minimize changes to AeroDyn. Internally AeroDyn stores a list of blades. This way we can keep all the storage of NumNodes x TotalNumBlades for now, where TotalNumBlades is the total number of blades for all rotors.

- Introduce "Rotor" types for convenient exchange with the glue code. For each rotor, we store the index mapping from the current rotor blades to the AeroDyn array of blades.
- List of blade is read as usual (TotalNumBlades are read)
- Tower definitions are read for each rotor
- Write outputs that are "per rotor" (e.g. RotorArea) will be postponed to future work and only supported for one rotor.  


**Driver implementation:**

- New driver (named AeroDynMultiRotor?), copy pasted from existing one. 
- New input file to define rotors and motions
- Setup meshes and update them with time (rigid body motion)


**Limitations:**

- Number of points per blades have to be the same for all blade (AeroDyn and FVW limitation). This can be progressively removed, first by storing up to the maximum number of nodes, and later by splitting data in types "per blade".
- Max one tower per rotor
- Number of nodes per tower constant for all rotors
- At first, write outputs might be per blade with no notion of rotor. This can be improved later on. Will be determined during implementation.



Registry changes
~~~~~~~~~~~~~~~~

The array dimensions are written explicitly but they will be allocatable 


**RotorInitInputs** (new)

.. code::

   ReKi     HubPosition             {3}
   ReKi     HubOrientation          {3}{3}
   ReKi     BladeRootPosition       {3}{NumBlade}
   ReKi     BladeRootOrientation    {3}{3}{NumBlade}
   IntKi    NumBlade                                ! NOTE: per rotor 
   IntKi    BladeID                 {NumBlade}  ! (new, index in AeroDyn blade array) 

**InitInput**

.. code::

   CHARACTER(1024)    InputFile 
   Logical            Linearize   
   IntKi              NumRotor  
   ReKi               Gravity    
   CHARACTER(1024)    RootName  
   RotorInitInputs    Rotors      {NumRotor}   ! (new)


**RotorInputs** (new)

.. code::

   MeshType   TowerMotion
   MeshType   HubMotion
   MeshType   BladeRootMotion     {NumBlade}
   MeshType   BladeMotion         {NumBlade}
   ReKi       InflowOnBlade    {3}{NumNodes}{NumBlade}
   ReKi       InflowOnTower    {3}{NumTwrNodes}
   ReKi       UserProp            {NumNodes}{NumBlade}

**Inputs**

.. code::

   RotorInputs   Rotors {NumRotor}     ! (new)
   ReKi          InflowWakeVel {3}{nOLAFPoints} (NOTE: will disappear for future OLAF)


**RotorOutputs** (new)

.. code::

   MeshType    TowerLoad   
   MeshType    BladeLoad    {NumBlade}


**Outputs**

.. code::

   RotorOutputs   Rotor          {NumRotor}   ! (new)
   ReKi           WriteOutput    {:}   


**Misc**:

Store mapping BladeID for each rotor.


Driver Input File
~~~~~~~~~~~~~~~~~


.. code::

    -----  AeroDyn MultiBlade Driver Input File  --------------------------------------
    Comment
    -----  Input Configuration  ----------------------------------------------------
    False          Echo         -  Echo input parameters to "<rootname>.ech"?
    "AeroDyn.dat"  AD_InputFile -  Name of the primary AeroDyn input file
    ----- Turbine Data  -----------------------------------------------------------
    1              NumRotors - Number of rotors
    ====== Rotor 1 ================================================================
    0,0,0          Origin          - x,y,z coordinates of rotor origin (m)
    0,0,0          OrientationInit - successive rotations (theta_x, theta_y, theta_z) defining initial orientation of the rotor frame from the global frame (e.g. roll, tilt, yaw) (deg)
    ----- Rotor 1 Motion -----------------------------------------------------------------
    1              MotionType      - Type of motion prescribed for this rotor {0: fixed, 1: Constant rotation, 2: time varying rotation, 3: arbitrary motion} (flag)
    0,0,0          Omega_r         - Rotational speed of rotor in rotor coordinates (rpm) [used only when MotionType=1]
    0,0,0          UnitOmega_r     - Unitary vector for the rotational speed of rotor in rotor coordinates (rpm) [used only when MotionType=2]
    ""             MotionFileName  - Filename containing rotor motion [used only when MotionType=2 or 3]
    ----- Rotor 1 Tower -----------------------------------------------------------------
    False          HasTower        - True if rotor has a tower (flag)
    0,0,0          TowerBase       - Coordinate of tower base in rotor coordinates [m] [used only wehn HasTower is True]
    0,0,0          TowerTop        - Coordinate of tower top in rotor coordinates [m] [used only wehn HasTower is True]
    ----- Rotor 1 Blades -----------------------------------------------------------------
    3              NumBlades          - Number of blades for current rotor (-)
    0,0,0          BladeOrigins_r     - Orign of blade 1 wrt. rotor origin in rotor coordinates (m)
    0,0,0          BladeOrigins_r     - Orign of blade 2 wrt. rotor origin in rotor coordinates (m)
    0,0,0          BladeOrigins_r     - Orign of blade 3 wrt. rotor origin in rotor coordinates (m)
    0  ,0,0        BladeOrientation_r -  successive rotations (theta_x, theta_y, theta_z) defining initial orientation of the blade frame from the rotor frame such that the "z" is along span, "y" along trailing edge without pitch (deg)
    120,0,0        BladeOrientation_r -  successive rotations (theta_x, theta_y, theta_z) defining initial orientation of the blade frame from the rotor frame such that the "z" is along span, "y" along trailing edge without pitch (deg)
    240,0,0        BladeOrientation_r -  successive rotations (theta_x, theta_y, theta_z) defining initial orientation of the blade frame from the rotor frame such that the "z" is along span, "y" along trailing edge without pitch (deg)
    1.3            BladeHubRad_b      - z-offset in blade coordinates of blade 1 where radial input data start (m)
    1.3            BladeHubRad_b      - z-offset in blade coordinates of blade 1 where radial input data start (m)
    1.3            BladeHubRad_b      - z-offset in blade coordinates of blade 1 where radial input data start (m)
    1              BladeFilenameIndex - Index of blade filename in AeroDyn input file
    2              BladeFilenameIndex - Index of blade filename in AeroDyn input file
    3              BladeFilenameIndex - Index of blade filename in AeroDyn input file
    -----  I/O Settings  -----------------------------------------------------------
      "ES15.8E2"  OutFmt       - Format used for text tabular output, excluding the time channel.  Resulting field should be 10 characters. (quoted string)


Motion File Type 2

.. code::

    #Time_[s], Angle_[rad], RotationalSpeed_[rad/s], RotationalAcceleration_[rad/s]
    0.0000   , 0.0        , 0.12                   , 0.001
    0.1000   , 0.2        , 0.12                   , 0.001
    0.2000   , 0.4        , 0.12                   , 0.001
    0.3000   , 0.6        , 0.12                   , 0.001



Motion File Type 3

.. code::

    # Time_[s] , x_[m], y_[m], z_[m] , e0_[-], e1_[-], e2_[-], e3_[-], xdot_[m/s], ydot_[m/s], zdot_[m/s], omega_x_g_[rad/s], omega_y_g_[rad/s], omega_z_g_[rad/s],  xddot_[m^2/s], yddot_[m^2/s] , zddot_[m^2/s],  alpha_x_g_[rad/s], alpha_y_g_[rad/s], alpha_z_g_[rad/s]

