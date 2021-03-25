
.. _ad_driver:

AeroDyn driver
===============


A standalone AeroDyn driver is provided to perform aerodynamic simulations of rigid turbines 
undergoing rigid body motion (fixed, sinusoidal, or arbitrary). 
The standalone AeroDyn driver code improves on the functionality previously
available in the separate wind turbine rotor-performance tool WT\_Perf.

Examples of applications are:

- Simulation of horizontal axis/ vertical axis wind turbines, kites, quad-rotors.
- Simulation with prescribed time series of wind speed, pitch, yaw, etc.
- Combined case analyses to compute the surfaces of power coefficient (C\ :sub:`P`), thrust coefficient (C\ :sub:`T`), and/or torque coefficient (C\ :sub:`Q`) as a function of tip-speed ratio (TSR) and blade-pitch angle. 
- Simulations with oscillatory motion of the tower base at different amplitudes and frequencies.

More details are provided below.


**Features:**

- Multiple rotors
- Arbitrary number of blades per rotors
- One tower per rotor 
- Arbitrary rigid body motion of each rotor (tower stays fixed for now). For convenience, we'll add simple prescribed motions as well, such as rotations (and maybe 6-DOF oscillations later)
   

**Limitations:**

- Number of points per blades have to be the same for all blade (AeroDyn and FVW limitation). This can be progressively removed, first by storing up to the maximum number of nodes, and later by splitting data in types "per blade".
- Max one tower per rotor
- Number of nodes per tower constant for all rotors
- At first, write outputs might be per blade with no notion of rotor. This can be improved later on. Will be determined during implementation.


Compiling the driver
--------------------

The compilation steps for the AeroDyn driver are similar to the ones of OpenFAST (see :numref:`installation`).  When using the CMake workflow, the driver is compiled automatically when running `make`. To compile only the driver, use `make aerodyn_driver`. The driver will be located in the folder `/build/modules/aerodyn/`. A Visual Studio solution is available for Windows, in the folder  `/vs-build/AeroDyn`.


.. _addm_driver-input-file:




Driver inputs and options
-------------------------

**Main concepts**
The driver supports:

- 2 kinds of turbine definitions: basic (HAWT), and advanced. 
- 2 kinds of inflow definition: basic (uniform power law), and advanced (InflowWind)
- 3 types of analyses: 1) one simulation of multiple turbines, 2) one simulation of one HAWT under time-dependent inputs, 3) combined cases of one HAWT

More details are provided below, where the different sections of the input file are described.


**Input configuration and analysis type**


The driver supports three kinds of analyses, but not all turbine formats and inflow types are available for each analysis type: 

- `AnalysisType=1`: Simulation of one or multiple rotors with basic (HAWT) or arbitrary geometries (HAWT/VAWT, quad-rotor, etc.), with basic or advanced wind inputs, and optional time-varying motion of the tower base, nacelle and individual pitch angles. Arbitrary motion or sinusoidal motion of the tower base are possible.
- `AnalysisType=2`: Simulation of a HAWT with basic time-varying wind, nacelle motion, pitch. 
- `AnalysisType=3`: Combined cases analyses of HAWT with basic steady wind. A table of cases is provided and run sequentially by the driver. 


An example of the input configuration is given below:
.. code::

    ----- Input Configuration -------------------------------------------------------
    False           Echo         - Echo input parameters to "<rootname>.ech"?
            3       AnalysisType - {1: multiple turbines, one simulation, 2: one turbine, one time-dependent simulation, 3: one turbine, combined cases}
           11.0     TMax         - Total run time [used only when AnalysisType/=3] (s)
            0.5     DT           - Simulation time step [used only when AnalysisType/=3] (s)
    "AD.dat"        AeroFile     - Name of the primary AeroDyn input file



**Inflow data**

The inflow can be provided in two ways:

- basic (`CompInflow=0`): uniform wind with a power law shear. The wind is defined using a reference height (`RefHt`), a powerlaw exponent (`PLExp`), and the wind speed at reference height (`HWindSpeed`). In some analyses types, the reference wind speed and power law can be defined as time series.
- advanced (`CompInflow=1`): the InflowWind module is used to compute the inflow, and all available options of InflowWind are then available. The user need to provide the name of the InflowWind input file (`InflowFile`)

An example of inputs is given below:

.. code::

    ----- Inflow Data ---------------------------------------------------------------
              0      CompInflow  - Compute inflow wind velocities (switch) {0=Steady Wind; 1=InflowWind}
    "unused"         InflowFile  - Name of the InflowWind input file [used only when CompInflow=1]
            9.0      HWindSpeed  - Horizontal wind speed   [used only when CompInflow=0 and AnalysisType=1] (m/s)
            140      RefHt       - Reference height for horizontal wind speed [used only when CompInflow=0]  (m)
           0.10      PLExp       - Power law exponent   [used only when CompInflow=0 and AnalysisType=1]                        (-)



**Turbine data**

The user specify the number of turbines as follows:

.. code:: 

    ----- Turbine Data --------------------------------------------------------------
    1               NumTurbines  - Number of turbines (should be 1 for AnalysisType=2 or AnalysisType=3)

As noted in the comment, the number of turbine should be 1 for `AnalysisType=2` and `AnalysisType=3`.
After the number of turbine is defined, the geometry and motion is defined for each turbine. Inputs for each turbine must have the suffix `(i)` where `i` is the turbine number.

An example of configuration with two wind turbines is given in :numref:`fig:MultiRotor`. The figure defines the different frames and origin associated with each turbine: the turbine base frame (t), nacelle frame (n), hub frame (h), and blade frames (b).
Prescribed motions of the turbine base occur at the turbine origin.
Yawing occurs around the :math:`z_n` axis,  the blade rotates about the :math:`x_h` axis, and pitching occurs around the individual :math:`z_b` axes. The definitions of the different frames are standardized when using a basic (HAWT) input format definition, and are arbitrary defined using the advanced input format. More details are given in the next paragraph.

.. figure:: figs/MultiRotor.png
   :width: 80%
   :name: fig:MultiRotor
           
   Definition of multiple rotors. 




**Turbine geometry definition**

Two turbine input formats are supported:

- basic (`BasicHAWTFormat=True`): Basic horizontal axis wind turbine (HAWT) format. In this format, the turbine geometry is entirely determined by the number of blades (`NumBlades`), the hub radius (`HubRad`), the hub height (`HubHt`), the overhang (`Overhang`), the shaft tilt (`ShftTilt`) and the precone (`Precone`). The definition of each parameter follows the OpenFAST convention.
  An example of basic input is given below:

.. code::

    ----- Turbine(1) Geometry -------------------------------------------------------
            True    BasicHAWTFormat(1) - Flag to switch between basic or generic input format {True: next 7 lines are basic inputs, False: Base/Twr/Nac/Hub/Bld geometry and motion must follow}
           0,0,0    BaseOriginInit(1) - Coordinates of turbine base in global coordinates (m)
               3    NumBlades(1)    - Number of blades (-)
              3.    HubRad(1)       - Hub radius (m)
          140.82513 HubHt(1)        - Hub height (m)
              -7    Overhang(1)     - Overhang (m)
              -6    ShftTilt(1)     - Shaft tilt (deg)
              -4    Precone(1)      - Blade precone (deg)


- advanced (`BasicHAWTFormat=False`): The position and orientation of the tower base, nacelle, hub, and individual blades can be arbitrarily defined. This can be used for HAWT and any other turbine concepts. 
  The definition of the different frames are given in :numref:`fig:MultiRotor`.
  The position (`BaseOriginInit`) and orientation (`BaseOriginInit`) of the turbine base frame are defined with respect to the OpenFAST global frame. Orientations are given using the values of three successive rotations (x-y-z Euler angles). If the base undergoes a motion, the orientation of the base frame will consists of the time-varying rotations followed by these initial rotations.

  A flag indicating whether the turbine has a tower is given on the next line (`HasTower`). This flag currently affects the VTK outputs. The user still has to provide tower inputs data in AeroDyn for each turbine (see :numref:`ad_inputs_multirot`).
  The next line indicates which projection AeroDyn is to use in its calculation. It is recommended to use `HAWTProjection=True` for HAWT, which is the default projection used in OpenFAST (projects on the coned-pitched axis). For other rotor concepts, set `HAWTprojection=False`.
  The following lines indicates the position and orientations of the tower, nacelle and hub. 

  The tower and the nacelle are defined with respect to the turbine base (t) origin and frame.
  The tower is assumed to end at the nacelle. 
  The tower station defined in the AeroDyn input file are assumed to be given with respect to the tower origin, unlike OpenFAST which uses ground/MSL as a reference (see :numref:`ad_inputs_multirot`).
  The hub is defined with respect to the nacelle origin and frame (n).

  The definitions of the blades follows, starting with the number of blades `NumBlades`. A rotor with zero blade is supported, and can be used to model a tower.
  If tower shadow/potential is used in AeroDyn, then the isolated tower will disturbed the flow of the vortex wake when OLAF is used.
  When BEM is used, the flow of the blades of a given turbine are disturbed only by the turbine's tower.
  The inputs for turbine `i` and blade `j` are labelled `(i_j)`.
  The origin `BldOrigin_h`) and orientation (`BldOrientation_h`) of each blade are given with respect to the hub origin and frame (h).
  Hub radius inputs (`BldHubRad_Bl`) are provided for convenience. They will effectively offset the blade origin along the :math:`z_b` axis.
  An example of input for an advanced geometry definition is given below. 
  This example corresponds to typical values for a 3-bladed upwind HAWT, with 6 degrees of tilt (-6 in OpenFAST) and -4 degrees of precone (going upstream)

.. code::

    ----- Turbine(1) Geometry -------------------------------------------------------
         False      BasicHAWTFormat(1) - Flag to switch between basic or generic input format {True: next 7 lines are basic inputs, False: Base/Twr/Nac/Hub/Bld geometry and motion must follow}
    0,0,0           BaseOriginInit(1)      - x,y,z coordinates of turbine base origin (m)
    0,0,0           BaseOrientationInit(1) - successive rotations (theta_x, theta_y, theta_z) defining initial orientation of the base frame from the global frame (e.g. roll, tilt, yaw) (deg)
    True            HasTower(1)            - True if turbine has a tower (flag)
    True            HAWTprojection(1)      - True if turbine is a horizontal axis turbine (for AeroDyn projections) (flag)
    0,0,0           TwrOrigin_t(1)         - Coordinate of tower base in base coordinates [used only when HasTower is True] (m)
    0,0,137         NacOrigin_t(1)         - x,y,z coordinates of nacelle origin (and tower top) from base, in base coordinates (m)
    -6.96,0.,3.82   HubOrigin_n(1)         - x,y,z coordinates of hub origin from nacelle origin, in nacelle coordinates (m)
    0,6,0           HubOrientation_n(1)    - successive rotations (theta_x, theta_y, theta_z) defining initial orientation of the hub frame from the nacelle frame (e.g. roll, tilt, yaw). The x axis needs to be aligned with the rotational speed. (deg)
    ----- Turbine(1) Blades -----------------------------------------------------------------
    3               NumBlades(1)          - Number of blades for current rotor (-)
    0,0,0           BldOrigin_h(1_1)      - Orign of blade 1 wrt. hub origin in hub coordinates (m)
    0,0,0           BldOrigin_h(1_2)      - Orign of blade 1 wrt. hub origin in hub coordinates (m)
    0,0,0           BldOrigin_h(1_3)      - Orign of blade 1 wrt. hub origin in hub coordinates (m)
    0  ,-4,0        BldOrientation_h(1_1) - successive rotations (theta_x, theta_y, theta_z) defining initial orientation of the blade frame from the hub frame such that the "z" is along span, "y" along trailing edge without pitch (azimuth, precone, pitch) (deg)
    120,-4,0        BldOrientation_h(1_2) - successive rotations (theta_x, theta_y, theta_z) defining initial orientation of the blade frame from the hub frame such that the "z" is along span, "y" along trailing edge without pitch (azimuth, precone, pitch) (deg)
    240,-4,0        BldOrientation_h(1_3) - successive rotations (theta_x, theta_y, theta_z) defining initial orientation of the blade frame from the hub frame such that the "z" is along span, "y" along trailing edge without pitch (azimuth, precone, pitch) (deg)
    3.0             BldHubRad_bl(1_1)     - z-offset in blade coordinates of blade 1 where radial input data start (m)
    3.0             BldHubRad_bl(1_2)     - z-offset in blade coordinates of blade 2 where radial input data start (m)
    3.0             BldHubRad_bl(1_3)     - z-offset in blade coordinates of blade 3 where radial input data start (m)



**Turbine motion definition**

The turbine motion inputs are different for a basic or advanced definition of the geometry. The definition of the turbine motion are only used when AnalysisType=1, but must always be present in the input file. 


Both definitions define the base motion in the same way.
The motion of the base may be: fixed (`BaseMotionType=0`), sinusoidal (`BaseMotionType=1`) or arbritrary (`BaseMotionType=2`). 
The turbine base motion is applied at each time step before applying the initial position and orientation of the turbine base.
A sinusoidal motion implies that one degree of freedom (`DegreeOfFreedom`) of the turbine base is moving according to a sine function of a given amplitude (`Amplitude`) and frequency (`Frequency`), with zero phase.
The 6 possible degrees of freedom corresponds to translations or rotations of the base frame in global coordinates (g).
An arbitrary motion is specified via a csv file (`BaseMotionFileName`) which contains 19 columns: time, 3 translations (global), three successive rotations (global), 3 translation velocities, 3 rotational velocities (omega, in global), 3 translational acceleration and 3 rotational accelerations (alpha, in global). Example of arbitrary input files are given in :numref:`ad_inputfiles_examples`.


An example of inputs is given below:

.. code::

    ----- Turbine(1) Motion [used only when AnalysisType=1] --------------------------
    1               BaseMotionType(1)      - Type of motion prescribed for this base {0: fixed, 1: Sinusoidal motion, 2: arbitrary motion} (flag)
    1               DegreeOfFreedom(1)     - {1:xg, 2:yg, 3:zg, 4:theta_xg, 5:theta_yg, 6:theta_zg} [used only when BaseMotionType=1] (flag)
    5.0             Amplitude(1)           - Amplitude of sinusoidal motion   [used only when BaseMotionType=1]
    0.1             Frequency(1)           - Frequency of sinusoidal motion   [used only when BaseMotionType=1]
    ""              BaseMotionFileName(1)  - Filename containing arbitrary base motion (19 columns: Time, x, y, z, theta_x, ..., alpha_z)  [used only when BaseMotionType=2]


The different inputs for the basic and advanced geometries are given below:

- basic: The motion of a basic turbine consist of a constant nacelle yaw (`NacYaw`), rotor speed (`RotSpeed`), blade pitch (`BldPitch`). 
  Examples are given below:

.. code::

    0               NacYaw(1)              - Yaw angle (about z_t) of the nacelle (deg)
    7               RotSpeed(1)            - Rotational speed of rotor in rotor coordinates (rpm)
    1               BldPitch(1)            - Blades pitch (deg)

- advanced: when an advanced geometry is provided, the motion sections contains detailled options for fine turbine of the nacelle motion, rotor motion and individual pitch motion. 
  The syntax for each of these motion consist of defining a type (fixed or time-varying), a value for the fixed case, or a file for the time-varying case.
  The input files are CSV files containing time, position, speed and acceletation. Examples of files are given in :numref:`ad_inputfiles_examples`.
  The rotation data in the CSV files are defined in rad and rad/s, whereas they are defined in deg and rpm in the driver input file.
  Examples are given below:
.. code::

    0               NacMotionType(1)      - Type of motion prescribed for the nacelle {0: fixed yaw, 1: time varying yaw angle} (flag)
    0               NacYaw(1)             - Yaw angle (about z_t) of the nacelle [user only when NacMotionType=0] (deg)
    "unused"        NacMotionFileName(1)  - Filename containing yaw motion [used only when NacMotionType=1]
    0               RotMotionType(1)        - Type of motion prescribed for this rotor {0: constant rotation, 1: time varying rotation} (flag)
    6.0             RotSpeed(1)             - Rotational speed of rotor in rotor coordinates [used only when RotorMotionType=0] (rpm)
    "unused"        RotMotionFileName(1)    - Filename containing rotor motion [used only when RotorMotionType=1]
    0               BldMotionType(1)        - Type of pitch motion prescribed for the blades {0: fixed, 1: time varying pitch} (flag)
    0               BldPitch(1_1)           - Blade 1 pitch [used only when BldMotionType=0] (deg)
    0               BldPitch(1_2)           - Blade 2 pitch [used only when BldMotionType=0] (deg)
    0               BldPitch(1_3)           - Blade 3 pitch [used only when BldMotionType=0] (deg)
    "unused" BldMotionFileName(1_1)  - Filename containing blade pitch motion [used only when BldMotionType=1]
    "unused" BldMotionFileName(1_2)  - Filename containing blade pitch motion [used only when BldMotionType=1]
    "unused" BldMotionFileName(1_3)  - Filename containing blade pitch motion [used only when BldMotionType=1]


General considerations for the advanced: 

- A turbine is assumed to consist of an optional tower, a nacelle, a hub, and multiple blades.
- Different frames and origins are defined: the turbine frame (t), the hub frame (h), and the blade frames (b). 
- The tower points are defined in the turbine coordinates
- The hub frame and origin is defined with respect to the turbine coordinates
- The blade frames and origins are defined with respect to the hub coordinates
- The blades are rigidly attached to the hub, and rotate around the x axis of the hub frame. 
- For each blade, the blade frame is such that the zb-axis is along the span, the yb axis is directed towards the "trailing edge" in the absence of pitch, and the xb-axis is directed towards the suction side in the absence of pitch.








  




The blades will rotate rigidly




          3.    HubRad(1)       - Hub radius (m)
      140.82513 HubHt(1)        - Hub height (m)
          -7    Overhang(1)     - Overhang (m)
          -6    ShftTilt(1)     - Shaft tilt (deg)
          -4    Precone(1)      - Blade precone (deg)
  or generic input format {True: next 7 lines are basic inputs, False: Base/Twr/Nac/Hub/Bld geometry and motion must follow}






.. _ad_inputs_multirot:

AeroDyn inputs for multiple turbines
------------------------------------

To minimize the impact of the multipleturbine implementation, the driver currently uses only one AeroDyn input file for all turbines. This means that the AeroDyn options are currently the same for all rotors.


Blades

Tower inputs for each turbines.

The tower station defined in the AeroDyn input file are assumed to be given with respect to the tower origin, unlike OpenFAST which uses ground/MSL as a reference).



.. _ad_inputfiles_examples:

Examples of driver input files
------------------------------








Driver Input Files
~~~~~~~~~~~~~~~~~~

Example of an aerodyn driver for a basic inflow, basic HAWT, and combined case analyses:

.. code::

    ----- AeroDyn Driver Input File -------------------------------------------------
    Driver input file for combined case analyses
    ----- Input Configuration -------------------------------------------------------
    False           Echo         - Echo input parameters to "<rootname>.ech"?
            3       AnalysisType - {1: multiple turbines, one simulation, 2: one turbine, one time-dependent simulation, 3: one turbine, combined cases}
           11.0     TMax         - Total run time [used only when AnalysisType/=3] (s)
            0.5     DT           - Simulation time step [used only when AnalysisType/=3] (s)
    "./OpenFAST_BAR_00_AeroDyn15.dat"        AeroFile - Name of the primary AeroDyn input file
    ----- Inflow Data ---------------------------------------------------------------
              0      CompInflow  - Compute inflow wind velocities (switch) {0=Steady Wind; 1=InflowWind}
    "unused"         InflowFile  - Name of the InflowWind input file [used only when CompInflow=1]
            9.0      HWindSpeed  - Horizontal wind speed   [used only when CompInflow=0 and AnalysisType=1] (m/s)
            140      RefHt       - Reference height for horizontal wind speed [used only when CompInflow=0]  (m)
           0.10      PLExp       - Power law exponent   [used only when CompInflow=0 and AnalysisType=1]                        (-)
    ----- Turbine Data --------------------------------------------------------------
    1               NumTurbines  - Number of turbines
    ----- Turbine(1) ----------------------------------------------------------------
            True    BasicHAWTFormat(1) - Flag to switch between basic or generic input format {True: next 7 lines are basic inputs, False: Base/Twr/Nac/Hub/Bld geometry and motion must follow}
           0,0,0    BaseOriginInit(1) - Coordinate of tower base in base coordinates (m)
               3    NumBlades(1)    - Number of blades (-)
              3.    HubRad(1)       - Hub radius (m)
          140.82513 HubHt(1)        - Hub height (m)
              -7    Overhang(1)     - Overhang (m)
              -6    ShftTilt(1)     - Shaft tilt (deg)
              -4    Precone(1)      - Blade precone (deg)
    ----- Turbine(1) Motion [used only when AnalysisType=1] --------------------------
    0               NacYaw(1)              - Yaw angle (about z_t) of the nacelle (deg)
    7               RotSpeed(1)            - Rotational speed of rotor in rotor coordinates (rpm)
    1               BldPitch(1)            - Blade 1 pitch (deg)
    1               BaseMotionType(1)      - Type of motion prescribed for this base {0: fixed, 1: Sinusoidal motion, 2: arbitrary motion} (flag)
    1               DegreeOfFreedom(1)     - {1:xt, 2:yt, 3:zt, 4:theta_xt, 5:theta_yt, 6:theta_zt} [used only when BaseMotionType=1] (flag)
    5.0             Amplitude(1)           - Amplitude of sinusoidal motion   [used only when BaseMotionType=1]
    0.1             Frequency(1)           - Frequency of sinusoidal motion   [used only when BaseMotionType=1]
    ""              BaseMotionFileName(1)  - Filename containing arbitrary base motion (19 columns: Time, x, y, z, theta_x, ..., alpha_z)  [used only when BaseMotionType=2]
    ----- Time-dependent Analysis [used only when AnalysisType=2] -------------------
    "unused"         TimeAnalysisFileName - Filename containing time series (6 column: Time, HWndSpeed, PLExp, RotSpd, Pitch, Yaw). 
    -----  Combined-Case Analysis [used only when AnalysisType=3] -------------------
             4  NumCases     - Number of cases to run
    HWndSpeed     PLExp     RotSpd       Pitch          Yaw     dT      Tmax  DOF  Amplitude Frequency 
    (m/s)        (-)          (rpm)        (deg)        (deg)   (s)     (s)    (-)     (-)     (Hz)
      08           0.0          6.         0.            0.     1.0     100     0      0        0 
      08           0.0          6.         0.            0.     1.0     100     0      0        0 
      09           0.1          7.         1.            0.     0.5     51      1      5.0      0.1 
      09           0.2          8.         2.            0.     0.51    52      1      2.0      0.2 
    ----- I/O Settings --------------------------------------------------------------
    "ES15.8E2"       OutFmt      - Format used for text tabular output, excluding the time channel.  Resulting field should be 10 characters. (quoted string)
    2                OutFileFmt  - Format for tabular (time-marching) output file (switch) {1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both}
    0                WrVTK       - VTK visualization data output: (switch) {0=none; 1=animation}
    2                VTKHubRad   - HubRadius for VTK visualization (m)
    -1,-1,-1,2,2,2   VTKNacDim   - Nacelle Dimension for VTK visualization x0,y0,z0,Lx,Ly,Lz (m)


Arbitrary base motion file:

.. code::

    time_[s] , x_[m]    , y_[m]    , z_[m]    , theta_x_[rad] , theta_y_[rad] , theta_z_[rad] , xdot_[m/s] , ydot_[m/s] , zdot_[m/s] , omega_x_g_[rad/s] , omega_y_g_[rad/s] , omega_z_g_[rad/s] , xddot_[m^2/s] , yddot_[m^2/s] , zddot_[m^2/s] , alpha_x_g_[rad/s] , alpha_y_g_[rad/s] , alpha_z_g_[rad/s]
    0.000000 , 0.000000 , 0.000000 , 0.000000 , 0.000000      , 0.000000      , 0.000000      , 0.000000   , 0.000000   , 10.053096  , 0.000000          , 0.000000          , 0.000000          , 0.000000      , 0.000000      , -0.000000     , 0.000000          , 0.000000          , 0.000000
    0.100000 , 0.000000 , 0.000000 , 0.963507 , 0.000000      , 0.000000      , 0.000000      , 0.000000   , 0.000000   , 8.809596   , 0.000000          , 0.000000          , 0.000000          , 0.000000      , 0.000000      , -24.344157    , 0.000000          , 0.000000          , 0.000000


Yaw motion file:

.. code::

    time_[s] , yaw_[rad] , yaw_rate_[rad/s] , yaw_acc_[rad/s^2]
    0.000000 , 0.000000  , 0.000000         , 0.000000
    0.100000 , 0.007277  , 0.212647         , 4.029093

Rotor motion file:

.. code::

    time_[s] , azimuth_[rad] , omega_[rad/s] , rotacc_[rad/s^2]
    0.000000 , 0.000000      , 0.000000      , 0.000000
    0.100000 , 0.000000      , 0.000000      , 0.000000

Pitch motion file:

.. code::

    time_[s] , pitch_[rad] , pitch_rate_[rad/s] , pitch_acc_[rad/s^2]
    0.000000 , 0.000000    , 0.000000           , 0.000000
    0.100000 , 0.000000    , 0.000000           , 0.000000
    0.200000 , 0.000000    , 0.000000           , 0.000000


