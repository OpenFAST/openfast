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

**AzimB1Up**    - Azimuth value to use for I/O when blade 1 points up (degrees)

**OverHang**    - Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)

**ShftGagL**    - Distance from rotor apex [3 blades] or teeter pin [2 blades] to shaft strain gages [positive for upwind rotors] (meters)

**ShftTilt**    - Rotor shaft tilt angle (degrees)

**NacCMxn**     - Downwind distance from the tower-top to the nacelle CM (meters)

**NacCMyn**     - Lateral  distance from the tower-top to the nacelle CM (meters)

**NacCMzn**     - Vertical distance from the tower-top to the nacelle CM (meters)

**NcIMUxn**     - Downwind distance from the tower-top to the nacelle IMU (meters)

**NcIMUyn**     - Lateral  distance from the tower-top to the nacelle IMU (meters)

**NcIMUzn**     - Vertical distance from the tower-top to the nacelle IMU (meters)

**Twr2Shft**    - Vertical distance from the tower-top to the rotor shaft (meters)

**TowerHt**     - Height of tower above ground level [onshore] or MSL [offshore] (meters)

**TowerBsHt**   - Height of tower base above ground level [onshore] or MSL [offshore] (meters)

**PtfmCMxt**    - Downwind distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)

**PtfmCMyt**    - Lateral distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)

**PtfmCMzt**    - Vertical distance from the ground level [onshore] or MSL [offshore] to the platform CM (meters)

**PtfmRefzt**   - Vertical distance from the ground level [onshore] or MSL [offshore] to the platform reference point (meters)



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

**TeetCDmp**    - Rotor-teeter rate-independent Coulomb-damping moment (N-m) [used only for 2 blades and when TeetMod=1]

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
For instance, a value of 5 will cause FAST to generate output only every fifth
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
