.. _SrvD-Input:

Input Files
===========

The user configures the servodynamics model parameters via a primary ServoDyn
input file, as well as separate input files for Structural control, and a
controller DLL.  *This information is incomplete and will be documented here
at a later date.*


Units
-----

ServoDyn uses the SI system (kg, m, s, N). Angles are assumed to be in
radians unless otherwise specified.

ServoDyn Primary Input File
----------------------------

The primary ServoDyn input file defines the modeling options for the controller.
This includes some DLL options, and Structural control options (typically a
tuned mass damper system). 


Simulation Control
~~~~~~~~~~~~~~~~~~

**Echo** [flag]

   Echo input data to <RootName>.ech

**DT**   [sec]

   Communication interval for controllers (or "default")


Pitch Control
~~~~~~~~~~~~~

**PCMode** [switch]

   Pitch control mode {0: none, 3: user-defined from routine PitchCntrl, 4:
   user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL}

**TPCOn** [sec]

   Time to enable active pitch control *[unused when* **PCMode==0** *]*

**TPitManS(1)** [sec]

   Time to start override pitch maneuver for blade 1 and end standard pitch
   control

**TPitManS(2)** [sec]

   Time to start override pitch maneuver for blade 2 and end standard pitch
   control

**TPitManS(3)** [sec]

   Time to start override pitch maneuver for blade 3 and end standard pitch
   control *[unused for 2 blades]*

**PitManRat(1)** [deg/s]

   Pitch rate at which override pitch maneuver heads toward final pitch angle
   for blade 1

**PitManRat(2)** [deg/s]

   Pitch rate at which override pitch maneuver heads toward final pitch angle
   for blade 2

**PitManRat(3)** [deg/s]

   Pitch rate at which override pitch maneuver heads toward final pitch angle
   for blade 3 *[unused for 2 blades]*

**BlPitchF(1)** [deg]

   Blade 1 final pitch for pitch maneuvers

**BlPitchF(2)** [deg]

   Blade 2 final pitch for pitch maneuvers

**BlPitchF(3)** [deg]

   Blade 3 final pitch for pitch maneuvers *[unused for 2 blades]*


Generator and Torque Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**VSContrl** [switch]

   Variable-speed control mode {0: none, 1: simple VS, 3: user-defined from
   routine UserVSCont, 4: user-defined from Simulink/Labview, 5: user-defined from
   Bladed-style DLL}

**GenModel** [switch]

   Generator model {1: simple, 2: Thevenin, 3: user-defined from routine
   UserGen} *[used only when* **VSContrl==0** *]*

**GenEff**   [\%]

   Generator efficiency *[ignored by the Thevenin and user-defined generator
   models]*

**GenTiStr** [flag]

   Method to start the generator {T: timed using TimGenOn, F: generator speed
   using SpdGenOn}

**GenTiStp** [Flag]

   Method to stop the generator {T: timed using TimGenOf, F: when generator
   power = 0}

**SpdGenOn** [rpm]

   Generator speed to turn on the generator for a startup (HSS speed) *[used
   only when* **GenTiStri==False** *]*

**TimGenOn** [sec]

   Time to turn on the generator for a startup *[used only when*
   **GenTiStr==True** *]*

**TimGenOf** [sec]

   Time to turn off the generator *[used only when* **GenTiStp==True** *]*


Simple Variable-speed Torque Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**VS_RtGnSp** [rpm]

   Rated generator speed for simple variable-speed generator control (HSS side)
   *[used only when* **VSContrl==1** *]*

**VS_RtTq**   [N-m]

   Rated generator torque/constant generator torque in Region 3 for simple
   variable-speed generator control (HSS side) *[used only when* **VSContrl==1**
   *]*

**VS_Rgn2K**  [N-m/rpm^2]

   Generator torque constant in Region 2 for simple variable-speed generator
   control (HSS side) *[used only when* **VSContrl==1** *]*

**VS_SlPc**   [\%]

   Rated generator slip percentage in Region 2 1/2 for simple variable-speed
   generator control *[used only when* **VSContrl==1** *]*


Simple Induction Generator
~~~~~~~~~~~~~~~~~~~~~~~~~~

**SIG_SlPc**     [\%]

   Rated generator slip percentage *[used only when* **VSContrl==0** *and*
   **GenModel==1** *]*

**SIG_SySp**     [rpm]

   Synchronous (zero-torque) generator speed *[used only when* **VSContrl==0**
   *and* **GenModel==1** *]*

**SIG_RtTq**     [N-m]

   Rated torque *[used only when* **VSContrl==0** *and* **GenModel==1** *]*

**SIG_PORt**     [-]

   Pull-out ratio (Tpullout/Trated) *[used only when* **VSContrl==0** *and*
   **GenModel==1** *]*


Thevenin-Equivalent Induction Generator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**TEC_Freq**     [Hz]

   Line frequency [50 or 60] *[used only when* **VSContrl==0** *and*
   **GenModel==2** *]*

**TEC_NPol**     [-]

   Number of poles [even integer > 0] *[used only when* **VSContrl==0** *and*
   **GenModel==2** *]*

**TEC_SRes**     [ohms]

   Stator resistance *[used only when* **VSContrl==0** *and* **GenModel==2** *]*

**TEC_RRes**     [ohms]

   Rotor resistance *[used only when* **VSContrl==0** *and* **GenModel==2** *]*

**TEC_VLL**      [volts]

   Line-to-line RMS voltage *[used only when* **VSContrl==0** *and*
   **GenModel==2** *]*

**TEC_SLR**      [ohms]

   Stator leakage reactance *[used only when* **VSContrl==0** *and*
   **GenModel==2** *]*

**TEC_RLR**      [ohms]

   Rotor leakage reactance *[used only when* **VSContrl==0** *and*
   **GenModel==2** *]*

**TEC_MR**       [ohms]

   Magnetizing reactance *[used only when* **VSContrl==0** *and* **GenModel==2**
   *]*


High-speed Shaft Brake
~~~~~~~~~~~~~~~~~~~~~~

**HSSBrMode**     [switch]

   HSS brake model {0: none, 1: simple, 3: user-defined from routine UserHSSBr,
   4: user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL}

**THSSBrDp**      [sec]

   Time to initiate deployment of the HSS brake

**HSSBrDT**       [sec]

   Time for HSS-brake to reach full deployment once initiated *[used only when*
   **HSSBrMode==1** *]*

**HSSBrTqF**      [N-m]

   Fully deployed HSS-brake torque


Nacelle-yaw Control
~~~~~~~~~~~~~~~~~~~

**YCMode**        [switch]

   Yaw control mode {0: none, 3: user-defined from routine UserYawCont, 4:
   user-defined from Simulink/Labview, 5: user-defined from Bladed-style DLL}

**TYCOn**         [sec]

   Time to enable active yaw control *[unused when* **YCMode==0** *]*

**YawNeut**       [deg]

   Neutral yaw position--yaw spring force is zero at this yaw

**YawSpr**        [N-m/rad]

   Nacelle-yaw spring constant

**YawDamp**       [N-m/(rad/s)]

   Nacelle-yaw damping constant

**TYawManS**      [sec]

   Time to start override yaw maneuver and end standard yaw control

**YawManRat**     [deg/s]

   Yaw maneuver rate (in absolute value)

**NacYawF**       [deg]

   Final yaw angle for override yaw maneuvers


.. _SrvD-AfC-inputs:

Aerodynamic Flow Control
~~~~~~~~~~~~~~~~~~~~~~~~
 
**AfCmode**       [switch]

   Airfoil control mode {0: none, 1: sine wave cycle, 4: user-defined from
   Simulink/Labview, 5: user-defined from Bladed-style DLL}

**AfC_Mean**      [-]

   Mean level for cosine cycling or steady value *[used only with*
   **AfCmode==1** *]*

**AfC_Amp**       [-]

   Amplitude for cosine cycling of flap signal *[used only with*
   **AfCmode==1** *]*

**AfC_Phase**     [deg]

   Phase relative to the blade azimuth (0 is vertical) for cosine cycling of
   flap signal *[used only with* **AfCmode==1** *]*

When **AfCmode==1**, the signal for the airfoil flow control is set by the
expression *AfC_Mean + p%AfC_Amp*cos( Azimuth + AfC_phase)* where the azimuth
is the azimuth of that particular blade (azimuth=0 is considered vertical).


.. _SrvD-CableControl-inputs:

Cable Control
~~~~~~~~~~~~~

Control of cable elements specified in either the MoorDyn or SubDyn modules can
be controlled through ServoDyn by a Bladed-style controller.  Each cable
receives a pair of controller channels, one for the requested cable length
change (DeltaL), and one for the cable length rate of change (DeltaLdot).  The
channel assignments are requested by the modules with the cable elements
(MoorDyn and/or SubDyn at present), and mapped to the appropriate control
channel.  A summary of which module requested the channels is available in the
summary file output from ServoDyn.  Up to 100 channel groups may be requested
when linking to a DLL, or 20 channel groups when linking to Simulink.

**CCmode**        [switch]

   Cable control mode {0: none, 4: user-defined from Simulink/Labview, 5:
   user-defined from Bladed-style DLL}.

   Each cable control channel group consists of a channel for DeltaL (requested
   cable length change) and a channel for DeltaLdot (cable length change
   rate) from the controller/Simulink interface.


.. _SrvD-StC-inputs:

Structural Control
~~~~~~~~~~~~~~~~~~

See :numref:`StC-Locations` for descriptions of the mounting locations for each
of the following options.

**NumBStC**      [integer]

   Number of blade structural controllers

**BStCfiles**      [-]

   Name of the files for blade structural controllers (quoted strings on one
   line) *[unused when* **NumBStC==0** *]*

**NumNStC**      [integer]

   Number of nacelle structural controllers

**NStCfiles**      [-]

   Name of the files for nacelle structural controllers (quoted strings on one
   line) *[unused when* **NumNStC==0** *]*

**NumTStC**      [integer]

   Number of tower structural controllers 

**TStCfiles**      [-]

   Names of the file for tower structural control damping (quoted strings on one
   line) *[unused when* **NumTStC==0** *]*

**NumSStC**   [integer]

   Number of substructure structural controllers

**SStCfiles**   [-]

   Name of the files for substructure structural controllers (quoted strings on one
   line) *[unused when* **NumSStC==0** *]*


Bladed Controller Interface
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**DLL_FileName**  [-]

   Name/location of the dynamic library {.dll [Windows] or .so [Linux]} in the
   Bladed-DLL format *[used only with Bladed Interface]*

**DLL_InFile**    [-]

   Name of input file sent to the DLL *[used only with Bladed Interface]*

**DLL_ProcName**  [-]

   Name of procedure in DLL to be called *[case sensitive; used only with DLL
   Interface]*

**DLL_DT**        [sec]

   Communication interval for dynamic library (or "default") *[used only with
   Bladed Interface]*

**DLL_Ramp**      [flag]

   Whether a linear ramp should be used between DLL_DT time steps [introduces
   time shift when true] *[used only with Bladed Interface]*

**BPCutoff**      [Hz]

   Cutoff frequency for low-pass filter on blade pitch from DLL *[used only with
   Bladed Interface]*

**NacYaw_North**  [deg]

   Reference yaw angle of the nacelle when the upwind end points due North
   *[used only with Bladed Interface]*

**Ptch_Cntrl**    [switch]

   Record 28: Use individual pitch control {0: collective pitch; 1: individual
   pitch control} *[used only with Bladed Interface]*

**Ptch_SetPnt**   [deg]

   Record  5: Below-rated pitch angle set-point *[used only with Bladed
   Interface]*

**Ptch_Min**      [deg]

   Record  6: Minimum pitch angle *[used only with Bladed Interface]*

**Ptch_Max**      [deg]

   Record  7: Maximum pitch angle *[used only with Bladed Interface]*

**PtchRate_Min**  [deg/s]

   Record  8: Minimum pitch rate (most negative value allowed) *[used only with
   Bladed Interface]*

**PtchRate_Max**  [deg/s]

   Record  9: Maximum pitch rate  *[used only with Bladed Interface]*

**Gain_OM**       [N-m/(rad/s)^2]

   Record 16: Optimal mode gain *[used only with Bladed Interface]*

**GenSpd_MinOM**  [rpm]

   Record 17: Minimum generator speed *[used only with Bladed Interface]*

**GenSpd_MaxOM**  [rpm]

   Record 18: Optimal mode maximum speed *[used only with Bladed Interface]*

**GenSpd_Dem**    [rpm]

   Record 19: Demanded generator speed above rated *[used only with Bladed
   Interface]*

**GenTrq_Dem**    [N-m]

   Record 22: Demanded generator torque above rated *[used only with Bladed
   Interface]*

**GenPwr_Dem**    [W]

   Record 13: Demanded power *[used only with Bladed Interface]*


Bladed Interface Torque-Speed Look-up table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**DLL_NumTrq**    [-]

   Record 26: No. of points in torque-speed
   look-up table {0 = none and use the optimal mode parameters; nonzero = ignore
   the optimal mode PARAMETERs by setting Record 16 to 0.0} *[used only with Bladed
   Interface]*
   The following 2 column table format is expected:

   +------------+------------+
   | GenSpd_TLU | GenTrq_TLU |
   |   (rpm)    |   (N-m)    |
   +------------+------------+


.. _SrvD-Outputs:

Output
~~~~~~

**SumPrint**      [flag]

   Print summary data to <RootName>.sum.  This file contains a summary of the
   inputs, and will give a detailed list of the communication channels with a
   Bladed-style controller when used.  This information may be helpful in
   debugging a controller, or verifying how ServoDyn is configured.

**OutFile**       [-]

   Switch to determine where output will be placed: {1: in module output file
   only; 2: in glue code output file only; 3: both} *(currently unused)*

**TabDelim**      [flag]

   Use tab delimiters in text tabular output file? *(currently unused)*

**OutFmt**        [-]

   Format used for text tabular output (except time).  Resulting field should be
   10 characters. (quoted string) *(currently unused)*

**TStart**        [sec]

   Time to begin tabular output *(currently unused)*

**OutList** section controls output quantities generated by
ServoDyn. Enter one or more lines containing quoted strings that in turn
contain one or more output parameter names. Separate output parameter
names by any combination of commas, semicolons, spaces, and/or tabs. If
you prefix a parameter name with a minus sign, “-”, underscore, “_”, or
the characters “m” or “M”, ServoDyn will multiply the value for that
channel by –1 before writing the data. The parameters are written in the
order they are listed in the input file. ServoDyn allows you to use
multiple lines so that you can break your list into meaningful groups
and so the lines can be shorter. You may enter comments after the
closing quote on any of the lines. Entering a line with the string “END”
at the beginning of the line or at the beginning of a quoted string
found at the beginning of the line will cause ServoDyn to quit scanning
for more lines of channel names.  If ServoDyn encounters an
unknown/invalid channel name, it warns the users but will remove the
suspect channel from the output file. Please refer to the ServoDyn tab in the
Excel file :download:`OutListParameters.xlsx <../../../OtherSupporting/OutListParameters.xlsx>`
for a complete list of possible output parameters.


