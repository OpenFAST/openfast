
Implementation
==============


Datatypes
----------

**InitInput**:

.. code::

  CHARACTER(1024) InputFile   - - - "Name of the input file" -
  CHARACTER(1024) RootName    - - - "RootName for writing output files"


**InitOutput**:

The module returns init outputs needed for AeroDyn and ServoDyn, respecting the conventions of ElastoDyn.

.. code::

  CHARACTER(ChanLen) WriteOutputHdr {:}   - - "Names of the output-to-file channels" -
  CHARACTER(ChanLen) WriteOutputUnt {:}   - - "Units of the output-to-file channels" -
  IntKi              NumBl              - - - "Number of blades on the turbine" -
  ReKi               BlPitch            - - - "Initial blade pitch angle" rad
  ReKi               TowerHeight        - - - "Tower Height" m
  ReKi               HubHt              - - - "Height of the hub" m
  ReKi               PlatformPos {6}      - - "Initial platform position (6 DOFs)"
  ReKi               HubRad             - - - "Preconed hub radius (distance from the rotor apex to the blade root)" m
  ReKi               BladeLength        - - - "Blade length (for AeroDyn)" m
  ReKi               RotSpeed           - - - "Initial or fixed rotor speed" rad/s
  LOGICAL            isFixed_GenDOF     - - - "Whether the generator is fixed or free" -




**Inputs**:

Note: the yaw rate is only used to setup the velocities on the meshes.

.. code::

  ReKi     AerTrq         -     - -   "Aerodynamic torque" N-m
  ReKi     HSSBrTrqC      -     - -   "High speed side brake torque" N-m
  ReKi     GenTrq         -     - -   "Electrical generator torque on HSS" N-m
  ReKi     BlPitchCom    {:}    - 2pi "Commanded blade pitch angles" rad
  ReKi     Yaw            -     - -   "Yaw angle"  rad
  ReKi     YawRate        -     - -   "Yaw rate"  rad/s
  

**Outputs**:


The module returns outputs needed for AeroDyn and ServoDyn, respecting the conventions of ElastoDyn and using the proper units.
Towerline, nacelle and hub are needed for AeroDyn. Platform mesh is for possible future implementation with HydroDyn.

Note: LSSTipPxa is obtained from the state, but modulo 2 pi

Note: additional scalar outputs may be needed by ServoDyn.

.. code::

  MeshType BladeRootMotion {:} - - "For AeroDyn/BeamDyn: motions at the blade roots"
  MeshType HubPtMotion       - - - "For AeroDyn and Lidar(InflowWind): motions of the hub"
  MeshType NacelleMotion     - - - "For AeroDyn14 & ServoDyn/TMD: motions of the nacelle."
  MeshType TowerLn2Mesh      - - - "Tower line2 mesh with positions/orientations/velocities/accelerations"	-
  MeshType PlatformPtMesh    - - - "Platform reference point positions/orientations/velocities/accelerations" -
  ReKi     LSSTipPxa       - - 2pi "Rotor azimuth angle (position)" radians
  ReKi     RotSpeed          - - - "Rotor azimuth angular speed" rad/s
  ReKi     WriteOutput     {:} - - "Data to be written to an output file: see WriteOutputHdr for names of each variable" "see WriteOutputUnt"
  ReKi     HSS_Spd           - - - "High-speed shaft (HSS) speed" rad/s


**States**:

.. code::

 R8Ki   QT   {1} - - "Current estimate of Q (displacement) for each degree of freedom" -
 R8Ki   QDT  {1} - - "Current estimate of QD (velocity)    for each degree of freedom"


**Misc**:

None anticipated at the moment.

**Parameters**:

.. code::


  R8Ki        DeltaT          - - - "Time step for module time integration"  s
  IntKi       IntMethod       - - - "Integration method {1: RK4, 2: AB4, or 3: ABM4}"  -
  ReKi        J_DT            - - - "Drivetrain inertia (blades+hub+shaft+generator) kgm^2
  ReKi        PtfmPitch       - - - "Static platform tilt angle"	rad
  MeshMapType mapPtf2Twr      - - - "Mesh mapping from Ptfm to Tower line"	-
  MeshMapType mapTwr2Nac      - - - "Mesh mapping from Tower to Nacelle"	-
  MeshMapType mapNac2Hub      - - - "Mesh mapping from Nacelle to Hub"	-
  LOGICAL     isFixed_GenDOF  - - - "whether the generator is fixed or free" -

..
  ReKi        GBoxEff         - - - "Gear box efficiency"	-

**InputFileInput**

.. code::

  R8Ki               DeltaT      -   -   -   "Time step for module time integration"  s
  IntKi              IntMethod   -   -   -   "Integration method {1: RK4, 2: AB4, or 3: ABM4}"  -
  LOGICAL            GenDOF      -   -   -   "whether the generator is fixed or free" -
  R8Ki               Azimuth     -   -   -   "Initial azimuth angle for blade 1"  deg
  ReKi               BlPitch     -   -   -   "Initial blade pitch angles" radians
  ReKi               RotSpeed    -   -   -   "Initial or fixed rotor speed" RPM
  ReKi               PtfmPitch   -   -   -   "Initial platform position (6 DOFs)"
  IntKi              NumBl       -   -   -   "Number of blades on the turbine" -
  ReKi               TipRad      -   -   -   "Preconed blade-tip radius (distance from the rotor apex to the blade tip)"  m
  ReKi               HubRad      -   -   -   "Preconed hub radius (distance from the rotor apex to the blade root)" m
  ReKi               PreCone     -   -   -   "Rotor precone angles" deg
  ReKi               OverHang    -   -   -   "Distance from yaw axis to rotor apex or teeter pin"  m
  ReKi               ShftTilt    -   -   -   "Rotor shaft tilt angle"  deg
  ReKi               Twr2Shft    -   -   -   "Vertical distance from the tower-top to the rotor shaft"  m
  ReKi               TowerHt     -   -   -   "Height of tower above ground level [onshore] or MSL [offshore]"  m
  ReKi               RotIner     -   -   -   "Hub inertia about teeter axis (2-blader) or rotor axis (3-blader)"  "kg m^2"
  ReKi               GenIner     -   -   -   "Generator inertia about HSS" "kg m^2"
  LOGICAL            SumPrint    -   -   -   "Print summary data to <RootName>.sum" -
  ReKi               GBRatio     -   -   -   "Gearbox ratio" -

..
  ReKi               GBoxEff     -   -   -   "Gearbox efficiency" %


Workflow of main routines
-------------------------

Init
~~~~

- Read input file
- Transfer InputFile data to parameters
- Set reference positions of meshes 
    - Set `PlatformPtMesh` (at (0,0,0))
    - Set `TowerL2Mesh` using two nodes at `PlatformPtMesh` and  at the `TowerHt`
    - Set `NacellePtMotion`, based on `NacYaw` and last point of `TowerL2Mesh`
    - Set `HubPtMotion`     based on geometry (`Twr2Shft`, `Tilt`, `OverHang`, and zero azimuth)
    - Set `BladeRootMotion`, distributing the blades azimuthally based on the number of blades.
- Define mesh mappings:
    - Set mesh mapping between `PlatformPtMesh` and `TowerL2Mesh` 
    - Set mesh mapping between `TowerL2Mesh` and `NacellePtMotion`
    - Set mesh mapping between `NacellePtMotion` and `HubPtMotion`
    - Set mesh mapping between `HubPtMotion` and `BladeRootMotion`

- Return quantities needed by AeroDyn and ServoDyn


UpdateStates
~~~~~~~~~~~~

If the generator degrees of freedom is on (`isFixed_GenDOF`) , this routine calls one of the time integration method, based on `IntMethod`, each calling the function `CalcConStateDerivative`.
Otherwise, compute states based on Eq. :eq:`sed_stateEqGenDOF`.

CalcConStateDerivative
~~~~~~~~~~~~~~~~~~~~~~

Returns derivative of states according to Eq. :eq:`sed_stateEq`.



CalcOutput
~~~~~~~~~~

- Set relative motions and successively update meshes using mapping:
    - Set relative motion of `NacellePtMotion` based on yaw and yawrate
    - Set relative motion of `HubPtMotion` using  :math:`\psi` &  :math:`\dot{\psi}`
    - Set relative motion of `BladeRootMotion` based on `BlPitchCom`
- Compute `Outputs` and `WriteOutputs` ("Azimuth", "RotSpeed", "RotAcc", "GenSpeed", "GenAcc") from states and inputs. 
