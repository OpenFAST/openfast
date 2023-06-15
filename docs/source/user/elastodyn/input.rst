.. _ed_input:

Input Files
===========

The user configures the structural model parameters via a primary ElastoDyn
input file, as well as separate input files for the tower and *other stuff that
will be documented here later.*

No lines should be added or removed from the input files.

Units
-----

ElastoDyn uses the SI system (kg, m, s, N). Angles are assumed to be in
radians unless otherwise specified.

ElastoDyn Primary Input File
----------------------------

The primary ElastoDyn input file defines modeling options and geometries for the
OpenFAST structure including the tower, nacelle, drivetrain, and blades (if
BeamDyn is not used).  It also sets the initial conditions for the structure.

Simulation Control
~~~~~~~~~~~~~~~~~~

Set the **Echo** flag to TRUE if you wish to have ElastoDyn echo the
contents of the ElastoDyn primary, airfoil, and blade input files (useful
for debugging errors in the input files). The echo file has the naming
convention of *OutRootFile.ED.ech*. **OutRootFile** is either
specified in the I/O SETTINGS section of the driver input file when
running ElastoDyn standalone, or by the OpenFAST program when running a
coupled simulation.

**Method**

**dT**

Degrees of Freedom
~~~~~~~~~~~~~~~~~~

**FlapDOF1**    - First flapwise blade mode DOF (flag)

**FlapDOF2**    - Second flapwise blade mode DOF (flag)

**EdgeDOF**     - First edgewise blade mode DOF (flag)

**TeetDOF**     - Rotor-teeter DOF (flag) [unused for 3 blades]

**DrTrDOF**     - Drivetrain rotational-flexibility DOF (flag)

**GenDOF**      - Generator DOF (flag)

**YawDOF**      - Yaw DOF (flag)

**TwFADOF1**    - First fore-aft tower bending-mode DOF (flag)

**TwFADOF2**    - Second fore-aft tower bending-mode DOF (flag)

**TwSSDOF1**    - First side-to-side tower bending-mode DOF (flag)

**TwSSDOF2**    - Second side-to-side tower bending-mode DOF (flag)

**PtfmSgDOF**   - Platform horizontal surge translation DOF (flag)

**PtfmSwDOF**   - Platform horizontal sway translation DOF (flag)

**PtfmHvDOF**   - Platform vertical heave translation DOF (flag)

**PtfmRDOF**    - Platform roll tilt rotation DOF (flag)

**PtfmPDOF**    - Platform pitch tilt rotation DOF (flag)

**PtfmYDOF**    - Platform yaw rotation DOF (flag)



Initial Conditions
~~~~~~~~~~~~~~~~~~

**OoPDefl**     - Initial out-of-plane blade-tip displacement (meters)

**IPDefl**      - Initial in-plane blade-tip deflection (meters)

**BlPitch(1)**  - Blade 1 initial pitch (degrees)

**BlPitch(2)**  - Blade 2 initial pitch (degrees)

**BlPitch(3)**  - Blade 3 initial pitch (degrees) [unused for 2 blades]

**TeetDefl**    - Initial or fixed teeter angle (degrees) [unused for 3 blades]

**Azimuth**     - Initial azimuth angle for blade 1 (degrees)

**RotSpeed**    - Initial or fixed rotor speed (rpm)

**NacYaw**      - Initial or fixed nacelle-yaw angle (degrees)

**TTDspFA**     - Initial fore-aft tower-top displacement (meters)

**TTDspSS**     - Initial side-to-side tower-top displacement (meters)

**PtfmSurge**   - Initial or fixed horizontal surge translational displacement of platform (meters)

**PtfmSway**    - Initial or fixed horizontal sway translational displacement of platform (meters)

**PtfmHeave**   - Initial or fixed vertical heave translational displacement of platform (meters)

**PtfmRoll**    - Initial or fixed roll tilt rotational displacement of platform (degrees)

**PtfmPitch**   - Initial or fixed pitch tilt rotational displacement of platform (degrees)

**PtfmYaw**     - Initial or fixed yaw rotational displacement of platform (degrees)

Turbine Configuration
~~~~~~~~~~~~~~~~~~~~~

**NumBl**       - Number of blades (-)

**TipRad**      - The distance from the rotor apex to the blade tip (meters)

**HubRad**      - The distance from the rotor apex to the blade root (meters)

**PreCone(1)**  - Blade 1 cone angle (degrees)

**PreCone(2)**  - Blade 2 cone angle (degrees)

**PreCone(3)**  - Blade 3 cone angle (degrees) [unused for 2 blades]

**HubCM**       - Distance from rotor apex to hub mass [positive downwind] (meters)

**UndSling**    - Undersling length [distance from teeter pin to the rotor apex] (meters) [unused for 3 blades]

**Delta3**      - Delta-3 angle for teetering rotors (degrees) [unused for 3 blades]

**AzimB1Up**    - Azimuth value to use for I/O when blade 1 points up (degrees); for floating MHK turbines, blade 1 will be pointed up (opposite gravity) when `AzimB1Up` = 0; the user can set `AzimB1Up` to 180 degrees to give the same azimuth convention relative to the tower for floating MHK turbines as for fixed MHK turbines

**OverHang**    - Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)

**ShftGagL**    - Distance from rotor apex [3 blades] or teeter pin [2 blades] to shaft strain gages [positive for upwind rotors] (meters)

**ShftTilt**    - Rotor shaft tilt angle (degrees)

**NacCMxn**     - Downwind distance from the tower-top to the nacelle CM (meters)

**NacCMyn**     - Lateral  distance from the tower-top to the nacelle CM (meters)

**NacCMzn**     - Vertical distance from the tower-top to the nacelle CM, typically negative for floating MHK turbines (meters)

**NcIMUxn**     - Downwind distance from the tower-top to the nacelle IMU (meters)

**NcIMUyn**     - Lateral  distance from the tower-top to the nacelle IMU (meters)

**NcIMUzn**     - Vertical distance from the tower-top to the nacelle IMU, typically negative for floating MHK turbines (meters)

**Twr2Shft**    - Vertical distance from the tower-top to the rotor shaft, typically negative for floating MHK turbines (meters)

**TowerHt**     - Height of tower relative to ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] (meters)

**TowerBsHt**   - Height of tower base relative to ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] (meters)

**PtfmCMxt**    - Downwind distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters)

**PtfmCMyt**    - Lateral distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters)

**PtfmCMzt**    - Vertical distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform CM (meters)

**PtfmRefzt**   - Vertical distance from the ground level [onshore], MSL [offshore wind or floating MHK], or seabed [fixed MHK] to the platform reference point (meters)



Mass and Inertia
~~~~~~~~~~~~~~~~

**TipMass(1)**  - Tip-brake mass, blade 1 (kg)

**TipMass(2)**  - Tip-brake mass, blade 2 (kg)

**TipMass(3)**  - Tip-brake mass, blade 3 (kg) [unused for 2 blades]

**HubMass**     - Hub mass (kg)

**HubIner**     - Hub inertia about rotor axis [3 blades] or teeter axis [2 blades] (kg m^2)

**GenIner**     - Generator inertia about HSS (kg m^2)

**NacMass**     - Nacelle mass (kg)

**NacYIner**    - Nacelle inertia about yaw axis (kg m^2)

**YawBrMass**   - Yaw bearing mass (kg)

**PtfmMass**    - Platform mass (kg)

**PtfmRIner**   - Platform inertia for roll tilt rotation about the platform CM (kg m^2)

**PtfmPIner**   - Platform inertia for pitch tilt rotation about the platform CM (kg m^2)

**PtfmYIner**   - Platform inertia for yaw rotation about the platform CM (kg m^2)



Blade
~~~~~

**BldNodes**    - Number of blade nodes (per blade) used for analysis (-)

**BldFile(1)**  - Name of file containing properties for blade 1 (quoted string)

**BldFile(2)**  - Name of file containing properties for blade 2 (quoted string)

**BldFile(3)**  - Name of file containing properties for blade 3 (quoted string) [unused for 2 blades]


Rotor-Teeter
~~~~~~~~~~~~

**TeetMod**     - Rotor-teeter spring/damper model {0: none, 1: standard, 2: user-defined from routine UserTeet} (switch) [unused for 3 blades]

**TeetDmpP**    - Rotor-teeter damper position (degrees) [used only for 2 blades and when TeetMod=1]

**TeetDmp**     - Rotor-teeter damping constant (N-m/(rad/s)) [used only for 2 blades and when TeetMod=1]

**TeetSStP**    - Rotor-teeter soft-stop position (degrees) [used only for 2 blades and when TeetMod=1]

**TeetHStP**    - Rotor-teeter hard-stop position (degrees) [used only for 2 blades and when TeetMod=1]

**TeetSSSp**    - Rotor-teeter soft-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]

**TeetHSSp**    - Rotor-teeter hard-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]



Drivetrain
~~~~~~~~~~

**GBoxEff**     - Gearbox efficiency (%)

**GBRatio**     - Gearbox ratio (-)

**DTTorSpr**    - Drivetrain torsional spring (N-m/rad)

**DTTorDmp**    - Drivetrain torsional damper (N-m/(rad/s))



Furling
~~~~~~~

**Furling**     - Read in additional model properties for furling turbine (flag) [must currently be FALSE)

**FurlFile**    - Name of file containing furling properties (quoted string) [unused when Furling=False]
An example of furling input file is given in :numref:`TF_ed_input-file`.



Tower
~~~~~

**TwrNodes**    - Number of tower nodes used for analysis (-)

**TwrFile**     - Name of file containing tower properties (quoted string)


.. _ED-Outputs:

Outputs
~~~~~~~

**SumPrint** [flag]  Set this value to TRUE if you want ElastoDyn to generate a
summary file with the name **OutFileRoot**.ED.sum*. **OutFileRoot** is specified
by the OpenFAST program when running a coupled simulation.

**OutFile** [switch] is currently unused.  The eventual purpose is to allow 
output from ElastoDyn to be written to a module output file (option 1), or the
main OpenFAST output file (option 2), or both.  At present this switch is
ignored.

**TabDelim** [flag] is currently unused.  Setting this to True will set the
delimeter for text files to the tab character for the ElastoDyn module
**OutFile**.

**OutFmt** [quoted string] is currently unused. ElastoDyn will use this string
as the numerical format specifier for output of floating-point values in its
local output specified by **OutFile**. The length of this string must not exceed
20 characters and must be enclosed in apostrophes or double quotes.  You may not
specify an empty string. To ensure that fixed-width column data align properly
with the column titles, you should ensure that the width of the field is 10
characters. Using an E, EN, or ES specifier will guarantee that you will never
overflow the field because the number is too big, but such numbers are harder to
read. Using an F specifier will give you numbers that are easier to read, but
you may overflow the field. Please refer to any Fortran manual for details for
format specifiers.

**TStart** [s] sets the start time for **OutFile**.  This is currenlty unused.

**DecFact** [-] This parameter sets the decimation factor for output. ElastoDyn
will output data to **OutFile** only once each DecFact integration time steps.
For instance, a value of 5 will cause ElastoDyn to generate output only every fifth
time step. This value must be an integer greater than zero.

**NTwGages** [-] The number of strain-gage locations along the tower indicates
the number of input values on the next line. Valid values are integers from 0 to
5 (inclusive).

**TwrGagNd** [-] The virtual strain-gage locations along the tower are assigned
to the tower analysis nodes specified on this line. Possible values are 1 to
TwrNodes (inclusive), where 1 corresponds to the node closest to the tower base
(but not at the base) and a value of TwrNodes corresponds to the node closest to
the tower top. The exact elevations of each analysis node in the undeflected
tower, relative to the base of the tower, are determined as follows:

   Elev. of node J = TwrRBHt + ( J – 1⁄2 ) • [ ( TowerHt + TwrDraft – TwrRBHt ) / TwrNodes ]
      (for J = 1,2,...,TwrNodes)

You must enter at least NTwGages values on this line.
If NTwGages is 0, this line will be skipped, but you must have a line taking up
space in the input file. You can separate the values with combinations of tabs,
spaces, and commas, but you may use only one comma between numbers.

**NBlGages** [-] specifies the number of strain-gague locations along the blade,
and indicates the number of input values expected in **BldGagNd**. This is only
used when the blade structure is modeled in ElastoDyn.

**BldGagNd** [-] specifies the virtual strain-gage locations along the blade
that should be output. Possible values are 1 to **BldNodes** (inclusive), where
1 corresponds to the node closest to the blade root (but not at the root) and a
value of BldNodes corresponds to the node closest to the blade tip. The node
locations are specified by the ElastoDyn blade input files.  You must enter at
least NBlGages values on this line. If NBlGages is 0, this line will be skipped,
but you must have a line taking up space in the input file. You can separate the
values with combinations of tabs, spaces, and commas, but you may use only one
comma between numbers. This is only used when the blade structure is modeled in
ElastoDyn.


The **OutList** section controls output quantities generated by
ElastoDyn. Enter one or more lines containing quoted strings that in turn
contain one or more output parameter names. Separate output parameter
names by any combination of commas, semicolons, spaces, and/or tabs. If
you prefix a parameter name with a minus sign, “-”, underscore, “_”, or
the characters “m” or “M”, ElastoDyn will multiply the value for that
channel by –1 before writing the data. The parameters are written in the
order they are listed in the input file. ElastoDyn allows you to use
multiple lines so that you can break your list into meaningful groups
and so the lines can be shorter. You may enter comments after the
closing quote on any of the lines. Entering a line with the string “END”
at the beginning of the line or at the beginning of a quoted string
found at the beginning of the line will cause ElastoDyn to quit scanning
for more lines of channel names. Blade and tower node-related quantities
are generated for the requested nodes identified through the
**BldGagNd** and **TwrGagNd** lists above. If ElastoDyn encounters an
unknown/invalid channel name, it warns the users but will remove the
suspect channel from the output file. Please refer to the ElastoDyn tab in the
Excel file :download:`OutListParameters.xlsx <../../../OtherSupporting/OutListParameters.xlsx>`
for a complete list of possible output parameters.

.. _ED-Nodal-Outputs:

.. include:: EDNodalOutputs.rst






.. _TF_ed_input-file:

ElastoDyn furl input file
-------------------------

This section describes the furl input file indicated by the input ``FurlFile`` from the ElastoDyn input file.
OpenFAST will only read this file if the model is designated as a furling machine (when Furling from the primary input
file is set to True).
The input file defines the geometry and stuctural properties of the rotor-furl and tail-furl. 

The rotor-furl and tail-turl coordinate systems and the geometrical inputs are described in
:numref:`ed_rfrl_coordsys` and 
:numref:`ed_tfrl_coordsys`, respectively.

An example of ElastoDyn furl input file is provided below:

.. code::
    
    ---------------------- FAST FURLING FILE ---------------------------------------
    Comment
    ---------------------- FEATURE FLAGS (CONT) ------------------------------------
    False       RFrlDOF     - Rotor-furl DOF (flag)
    True        TFrlDOF     - Tail-furl DOF (flag)
    ---------------------- INITIAL CONDITIONS (CONT) -------------------------------
       0.0      RotFurl     - Initial or fixed rotor-furl angle (deg)
       0.0      TailFurl    - Initial or fixed tail-furl angle (deg)
    ---------------------- TURBINE CONFIGURATION (CONT) ----------------------------
       0.1      Yaw2Shft    - Lateral distance from the yaw axis to the rotor shaft (m)
       0.0      ShftSkew    - Rotor shaft skew angle (deg)
    0., 0., 0.  RFrlCM_n    - Position of the CM of the structure that furls with the rotor [not including rotor] from the tower-top, in nacelle coordinates (m)
    1.7,0.1,0.  BoomCM_n    - Postion of the tail boom CM from the tower top, in nacelle coordinates (m) 
    0., 0., 0.  TFinCM_n    - Position of tail fin CM from the tower top, in nacelle coordinates (m) 
    0., 0., 0.  RFrlPnt_n   - Position of an arbitrary point on the rotor-furl axis from the tower top, in nacelle coordinates (m)
       0.0      RFrlSkew    - Rotor-furl axis skew angle (deg)
       0.0      RFrlTilt    - Rotor-furl axis tilt angle (deg)
    0.3, 0., 0. TFrlPnt_n   - Position of an arbitrary point on the tail-furl axis from the tower top, in nacelle coordinates (m)
     -45.2      TFrlSkew    - Tail-furl axis skew angle (deg)
      78.7      TFrlTilt    - Tail-furl axis tilt angle (deg)
    ---------------------- MASS AND INERTIA (CONT) ---------------------------------
       0.0      RFrlMass    - Mass of structure that furls with the rotor [not including rotor] (kg)
      86.8      BoomMass    - Tail boom mass (kg)
       0.0      TFinMass    - Tail fin mass (kg)
       0.0      RFrlIner    - Inertia of the structure that furls with the rotor about the rotor-furl axis (kg m^2) [not including rotor]
     264.7      TFrlIner    - Tail boom inertia about tail-furl axis (kg m^2)
    ---------------------- ROTOR-FURL ----------------------------------------------
       0        RFrlMod     - Rotor-furl spring/damper model {0: none, 1: standard, 2:user-defined routine} (switch)
       0.0      RFrlSpr     - Rotor-furl spring constant (N-m/rad) [used only when RFrlMod=1]
       0.0      RFrlDmp     - Rotor-furl damping constant (N-m/(rad/s)) [used only when RFrlMod=1]
       0.0      RFrlUSSP    - Rotor-furl up-stop spring position (deg) [used only when RFrlMod=1]
       0.0      RFrlDSSP    - Rotor-furl down-stop spring position (deg) [used only when RFrlMod=1]
       0.0      RFrlUSSpr   - Rotor-furl up-stop spring constant (N-m/rad) [used only when RFrlMod=1]
       0.0      RFrlDSSpr   - Rotor-furl down-stop spring constant (N-m/rad) [used only when RFrlMod=1]
       0.0      RFrlUSDP    - Rotor-furl up-stop damper position (deg) [used only when RFrlMod=1]
       0.0      RFrlDSDP    - Rotor-furl down-stop damper position (deg) [used only when RFrlMod=1]
       0.0      RFrlUSDmp   - Rotor-furl up-stop damping constant (N-m/(rad/s)) [used only when RFrlMod=1]
       0.0      RFrlDSDmp   - Rotor-furl down-stop damping constant (N-m/(rad/s)) [used only when RFrlMod=1]
    ---------------------- TAIL-FURL -----------------------------------------------
       1        TFrlMod     - Tail-furl spring/damper model {0: none, 1: standard, 2:user-defined routine} (switch)
       0.0      TFrlSpr     - Tail-furl spring constant (N-m/rad) [used only when TFrlMod=1]
      10.0      TFrlDmp     - Tail-furl damping constant (N-m/(rad/s)) [used only when TFrlMod=1]
      85.0      TFrlUSSP    - Tail-furl up-stop spring position (deg) [used only when TFrlMod=1]
       3.0      TFrlDSSP    - Tail-furl down-stop spring position (deg) [used only when TFrlMod=1]
       1.0E3    TFrlUSSpr   - Tail-furl up-stop spring constant (N-m/rad) [used only when TFrlMod=1]
       1.7E4    TFrlDSSpr   - Tail-furl down-stop spring constant (N-m/rad) [used only when TFrlMod=1]
      85.0      TFrlUSDP    - Tail-furl up-stop damper position (deg) [used only when TFrlMod=1]
       0.0      TFrlDSDP    - Tail-furl down-stop damper position (deg) [used only when TFrlMod=1]
       1.0E3    TFrlUSDmp   - Tail-furl up-stop damping constant (N-m/(rad/s)) [used only when TFrlMod=1]
     137.0      TFrlDSDmp   - Tail-furl down-stop damping constant (N-m/(rad/s)) [used only when TFrlMod=1]






*Feature Flags*

**RFrlDOF**
The rotor-furl DOF will be enabled when this is True. The initial rotor-furl angle is specified with
RotFurl. If RFrlDOF is disabled, the rotor-furl angle will be fixed at RotFurl. (flag)

**TFrlDOF**
The tail-furl DOF will be enabled when this is True. The initial tail-furl angle is specified with
TailFurl. If TFrlDOF is disabled, the tail-furl angle will be fixed at TailFurl. (flag)


*Initial Conditions*

**RotFurl**
This is the fixed or initial rotor-furl angle. It is positive about the rotor-furl axis as shown in
:numref:`figTFAxes`. The rotor-furl axis is defined through input ``RFrlPnt_n``
RFrlSkew, and RFrlTilt below. This value must be greater than -180 and less than or equal to
180 degrees. (deg)

**TailFurl**
This is the fixed or initial tail-furl angle. It is positive about the tail-furl axis as shown in :numref:`figTFAxes`. 
The tail-furl axis is defined through inputs ``TFrlPnt_n``, ``TFrlSkew``,
and ``TFrlTilt`` below. This value must be greater than -180 and less than or equal to 180 degrees.
(deg)




*Turbine Configuration*

Inputs ``RFrlPnt_n``, ``RFrlSkew``, and ``RFrlTilt`` define the orientation of the rotor-furl axis and associated DOF, ``RFrlDOF``. 
Inputs ``TFrlPnt_n``, ``TFrlSkew``, and ``TFrlTilt`` define the orientation of the tail-furl axis and associated DOF, ``TFrlDOF``. 
See :numref:`figTFAxes`.


**Yaw2Shft**
This is the lateral offset distance from the yaw axis to the intersection of the rotor shaft axis with
the yn-/zn-plane. The distance is measured parallel to the yn-axis. It is positive to the left when
looking downwind as shown in :numref:`figTFFurl`. 
For turbines with rotor-furl, this distance defines the configuration at a furl angle of zero. (m)

**ShftSkew**
This is the skew angle of the rotor shaft in the nominally horizontal plane. Positive skew acts like
positive nacelle yaw as shown in :numref:`figTFFurl`; however, ``ShftSkew`` should only be used to skew the
shaft a few degrees away from the zero-yaw position and must not be used as a replacement for
the yaw angle. This value must be between -15 and 15 degrees (inclusive). 
For turbines with rotor-furl, this angle defines the configuration at a furl angle of zero. (deg)

**RFrlCM_n**
Position of the center of mass of the structure that furls with the rotor
(not including the rotor-reference input ``RFrlMass``) measured from the tower top
and expressed in the nacelle coordinate system.
See :numref:`figTFFurl`. 
For turbines with rotor-furl, this position defines the configuration at a furl angle of zero. (m)

**BoomCM_n**
Position of the tail boom mass center (reference input ``BoomMass``) with respect to the tower top,
expressed in the nacelle coordinate system.
See :numref:`figTFGeom`. 
For turbines with tail-furl, this distance defines the configuration at a furl angle of zero. (m)

**TFinCM_n**
Position of the tail fin mass center (reference input ``TFinMass``) with respect to the top,
expressed in the nacelle coordinate system.
See :numref:`figTFGeom`. 
For turbines with tail-furl, this distance defines the configuration at a furl angle of zero. (m)

**RFrlPnt_n**
Position of an arbitrary point on the rotor-furl axis expressed from the tower top and 
in the nacelle coordinate system.
See :numref:`figTFAxes`. (m)

**RFrlSkew**
This is the skew angle of the rotor-furl axis in the nominally horizontal plane. Positive skew
orients the nominal horizontal projection of the rotor-furl axis about the zn-axis. 
See :numref:`figTFAxes`. 
This value must be greater than -180
and less than or equal to 180 degrees. (deg)

**RFrlTilt**
This is the tilt angle of the rotor-furl axis from the nominally horizontal plane. 
This value must be between -90 and 90 degrees (inclusive). 
See :numref:`figTFAxes`. (deg)

**TFrlPnt_n**
Position from the tower top to an arbitrary point on the tail-furl axis, in nacelle coordinates.
See :numref:`figTFAxes`. (m)

**TFrlSkew**
This is the skew angle of the tail-furl axis in the nominally horizontal plane. 
Positive skew orients the nominal horizontal projection of the tail-furl axis about the zn-axis.  
See :numref:`figTFAxes`. 
This value must be greater than -180 and less than or equal to 180 degrees. 
(deg)

**TFrlTilt**
This is the tilt angle of the tail-furl axis from the nominally horizontal plane. 
See :numref:`figTFAxes`. 
This value must be between -90 and 90 degrees (inclusive).  
(deg)


*Mass and Inertia*

**RFrlMass**
This is the mass of the structure that furls with the rotor (not including the rotor). The center of
this mass is located at the point specified by input ``RFrlCM_n``
relative to the tower-top at a rotor-furl angle of zero. It includes everything that furls with the
rotor excluding the rotor (blades, hub, and tip brakes). This value must not be negative. (kg)

**BoomMass**
This is the mass of the tail boom. The center of the tail boom mass is located at the point specified
by input ``BoomCM_n`` relative to the tower-top at a tail-furl angle
of zero. It includes everything that furls with the tail except the tail fin (see next input). This
value must not be negative. (kg)

**TFinMass**
This is the mass of the tail fin. The center of the tail fin mass is located at the point specified by
input ``TFinCM_n`` relative to the tower-top at a tail-furl angle of
zero. TFinMass and BoomMass combined should include everything that furls with the tail.
This value must not be negative. (kg)

**RFrlIner**
This is the moment of inertia of the structure that furls with the rotor (not including the rotor)
about the rotor-furl axis. It includes all mass contained in ``RFrlMass``. This value must be greater
than: ``RFrlMass*d^2`` where d is the perpendicular distance between rotor-furl axis and C.M. of the structure
that furls with the rotor [not including the rotor]. (kg·m2)

**TFrlIner**
This is the tail boom moment of inertia about the tail-furl axis. It includes all mass contained in
BoomMass. This value must be greater than: ``BoomMass*d^2`` where d is the perpendicular distance between
tail-furl axis and tail boom C.M. (kg·m2)



*Rotor-Furl*


The rotor-furl bearing can be an ideal bearing with
no friction by setting ``RFrlMod`` to 0; by setting
``RFrlMod`` to 1, it also has a standard model that
includes a linear spring and linear damper,
as well as up- and down-stop springs, and up-
and down-stop dampers. 
The formulae are provided in :numref:`ed_rtfrl_theory`.
ElastoDyn models the stop
springs with a linear function of rotor-furl deflection.
The rotor-furl stops start at a specified angle and work
as a linear spring based on the deflection past the stop
angles. The rotor-furl dampers are linear functions of
the furl rate and start at the specified up-stop and
down-stop angles. These dampers are bidirectional,
resisting motion equally in both directions once past
the stop angle.

A user-defined rotor-furl spring and damper model
is also available. To use it, set `RFrlMod` to 2 and
create a subroutine entitled `UserRFrl()` with the
parameters ``RFrlDef``, ``RFrlRate``, ``DirRoot``, ``ZTime``, and
``RFrlMom``:

- ``RFrlDef``: Current rotor-furl angular deflection in radians (input)
- ``RFrlRate``: Current rotor-furl angular rate in rad/sec (input)
- ``ZTime``: Current simulation time in sec (input)
- ``DirRoot``: Simulation root name including the full path to the current working director (input)
- ``RFrlMom``: Rotor-furl moment in N·m (output)

The source file ``ED_UserSubs.f90`` contains a dummy
``UserRFrl()`` routine; replace it with your own and
rebuild ElastoDyn. 


**RFrlMod**
The rotor-furl springs and dampers can be modeled three ways. For a value of 0 for ``RFrlMod``,
there will be no rotor-furl spring nor damper and the moment normally produced will be set to
zero. A ``RFrlMod`` of 1 will invoke simple spring and damper models using the inputs provided
below as appropriate coefficients. If ``RFrlMod`` is set to 2, ElastoDyn will call the routine ``UserRFrl()``
to compute the rotor-furl spring and damper moments. You should replace the dummy routine
supplied with the code with your own, which will need to be linked with the rest of ElastoDyn. Using
values other than 0, 1, or 2 will cause ElastoDyn to abort. (switch)

**RFrlSpr**
The linear rotor-furl spring restoring moment is proportional to the rotor-furl deflection through
this constant. This value must not be negative and is only used when ``RFrlMod`` is set to 1.
(N·m/rad)

**RFrlDmp**
The linear rotor-furl damping moment is proportional to the rotor-furl rate through this constant.
This value must not be negative and is only used when ``RFrlMod`` is set to 1. (N·m/(rad/s))

**RFrlCDmp**
This Coulomb-friction damping moment resists rotor-furl motion, but it is a constant that is not
proportional to the rotor-furl rate. However, if the rotor-furl rate is zero, the damping is zero.
This value must not be negative and is only used when ``RFrlMod`` is set to 1. (N·m)

**RFrlUSSP**
The rotor-furl up-stop spring is effective when the rotor-furl deflection exceeds this value. This
value must be greater than -180 and less than or equal to 180 degrees and is only used when
``RFrlMod`` is set to 1. (deg)

**RFrlDSSP**
The rotor-furl down-stop spring is effective when the rotor-furl deflection exceeds this value. This
value must be greater than -180 and less than or equal to ``RFrlUSSP`` degrees and is only used
when ``RFrlMod`` is set to 1. (deg)

**RFrlUSSpr**
The linear rotor-furl up-stop spring restoring moment is proportional to the rotor-furl up-stop
deflection by this constant and is effective when the rotor-furl deflection exceeds ``RFrlUSSP``.
This value must not be negative and is only used when ``RFrlMod`` is set to 1. (N·m/rad)

**RFrlDSSpr**
The linear rotor-furl down-stop spring restoring moment is proportional to the rotor-furl down-
stop deflection by this constant and is effective when the rotor-furl deflection exceeds ``RFrlDSSP``.
This value must not be negative and is only used when ``RFrlMod`` is set to 1. (N·m/rad)

**RFrlUSDP**
The rotor-furl up-stop damper is effective when the rotor-furl deflection exceeds this value. This
value must be greater than -180 and less than or equal to 180 degrees and is only used when
``RFrlMod`` is set to 1. (deg)

**RFrlDSDP**
The rotor-furl down-stop damper is effective when the rotor-furl deflection exceeds this value.
This value must be greater than -180 and less than or equal to ``RFrlUSDP`` degrees and is only
used when ``RFrlMod`` is set to 1. (deg)

**RFrlUSDmp**
The linear rotor-furl up-stop damping moment is proportional to the rotor-furl rate by this constant
and is effective when the rotor-furl deflection exceeds ``RFrlUSDP``. This value must not be
negative and is only used when ``RFrlMod`` is set to 1. (N·m/(rad/s))

**RFrlDSDmp**
The linear rotor-furl down-stop damping restoring moment is proportional to the rotor-furl rate by
this constant and is effective when the rotor-furl deflection exceeds ``RFrlDSDP``. This value must
not be negative and is only used when ``RFrlMod`` is set to 1. (N·m/(rad/s))





*Tail-Furl*


The tail-furl bearing can be an ideal bearing with
no friction by setting ``TFrlMod`` to 0; by setting
``TFrlMod`` to 1, it also has a standard model that
includes a linear spring and damper ,
as well as up- and down-stop springs, and up-
and down-stop dampers. 
The formulae are provided in :numref:`ed_rtfrl_theory`.
ElastoDyn models the stop
springs with a linear function of tail-furl deflection.
The tail-furl stops start at a specified angle and work as
a linear spring based on the deflection past the stop
angles. The tail-furl dampers are linear functions of
the furl rate and start at the specified up-stop and
down-stop angles. These dampers are bidirectional,
resisting motion equally in both directions once past
the stop angle.

A user-defined tail-furl spring and damper model
is also available. To use it, set ``TFrlMod`` to 2 and
create a subroutine entitled ``UserTFrl()`` with the
arguments ``TFrlDef``, ``TFrlRate``, ``ZTime``, ``DirRoot``, and
``TFrlMom``:

- ``TFrlDef``: Current tail-furl angular deflection in radians (input)
- ``TFrlRate``: Current tail-furl angular rate in rad/sec (input)
- ``ZTime``: Current simulation time in sec (input)
- ``DirRoot``: Simulation root name including the full path to the current working directory (input)
- ``TFrlMom``: Tail-furl moment in N.m (output)

The source file ``ED_UserSubs.f90`` contains a dummy
``UserTFrl()`` routine; replace it with your own and
rebuild ElastoDyn. 


**TFrlMod**
The tail-furl springs and dampers can be modeled three ways. For a value of 0 for ``TFrlMod``, there
will be no tail-furl spring nor damper and the moment normally produced will be set to zero. A
``TFrlMod`` of 1 will invoke simple spring and damper models using the inputs provided below as
appropriate coefficients. If you set ``TFrlMod`` to 2, ElastoDyn will call the routine ``UserTFrl()`` to
compute the tail-furl spring and damper moments. You should replace the dummy routine
supplied with the code with your own, which will need to be linked with the rest of ElastoDyn. Using
values other than 0, 1, or 2 will cause ElastoDyn to abort. (switch)

**TFrlSpr**
The linear tail-furl spring restoring moment is proportional to the tail-furl deflection through this
constant. This value must not be negative and is only used when ``TFrlMod`` is set to 1. (N·m/rad)

**TFrlDmp**
The linear tail-furl damping moment is proportional to the tail-furl rate through this constant. This
value must not be negative and is only used when ``TFrlMod`` is set to 1. (N·m/(rad/s))

**TFrlCDmp**
This Coulomb-friction damping moment resists tail-furl motion, but it is a constant that is not
proportional to the tail-furl rate. However, if the tail-furl rate is zero, the damping is zero. This
value must not be negative and is only used when ``TFrlMod`` is set to 1. (N·m)

**TFrlUSSP**
The tail-furl up-stop spring is effective when the tail-furl deflection exceeds this value. This value
must be greater than -180 and less than or equal to 180 degrees and is only used when ``TFrlMod`` is
set to 1. (deg)

**TFrlDSSP**
The tail-furl down-stop spring is effective when the tail-furl deflection exceeds this value. This
value must be greater than -180 and less than or equal to ``TFrlUSSP`` degrees and is only used
when ``TFrlMod`` is set to 1. (deg)

**TFrlUSSpr**
The linear tail-furl up-stop spring restoring moment is proportional to the tail-furl up-stop
deflection by this constant and is effective when the tail-furl deflection exceeds ``TFrlUSSP``. This
value must not be negative and is only used when ``TFrlMod`` is set to 1. (N·m/rad)

**TFrlDSSpr**
The linear tail-furl down-stop spring restoring moment is proportional to the tail-furl down-stop
deflection by this constant and is effective when the tail-furl deflection exceeds ``TFrlDSSP``. This
value must not be negative and is only used when ``TFrlMod`` is set to 1. (N·m/rad)

**TFrlUSDP**
The tail-furl up-stop damper is effective when the tail-furl deflection exceeds this value. This
value must be greater than -180 and less than or equal to 180 degrees and is only used when
``TFrlMod`` is set to 1. (deg)

**TFrlDSDP**
The tail-furl down-stop damper is effective when the tail-furl deflection exceeds this value. This
value must be greater than -180 and less than or equal to ``TFrlUSDP`` degrees and is only used
when ``TFrlMod`` is set to 1. (deg)

**TFrlUSDmp**
The linear tail-furl up-stop damping moment is proportional to the tail-furl rate by this constant
and is effective when the tail-furl deflection exceeds ``TFrlUSDP``. This value must not be negative
and is only used when ``TFrlMod`` is set to 1. (N·m/(rad/s))

**TFrlDSDmp**
The linear tail-furl down-stop damping restoring moment is proportional to the tail-furl rate by this
constant and is effective when the tail-furl deflection exceeds ``TFrlDSDP``. This value must not be
negative and is only used when ``TFrlMod`` is set to 1. (N·m/(rad/s))





