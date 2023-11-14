.. _sd_input-files:

Input Files
===========

The user specifies the substructure model parameters, including its
geometry and properties, via a primary SubDyn input file. When used in
stand-alone mode, an additional driver input file is required. This
driver file specifies inputs normally provided to SubDyn by FAST,
including motions of the TP reference point.

No lines should be added or removed from the input files, except in
tables where the number of rows is specified.

Additional input files containing soil-structure information (*SSIfile*)
can be provided by the user specifying their paths in the main SubDyn
input file under the section titled *BASE REACTION JOINTS*.

Units
-----

SubDyn uses the SI system (kg, m, s, N). Angles are assumed to be in
radians unless otherwise specified.

.. _sd_driver-input-file:

SubDyn Driver Input File
-------------------------

The driver input file is only needed for the stand-alone version of
SubDyn and contains inputs that are normally set by FAST, and that are
necessary to control the simulation for uncoupled models. It is possible
to provide per-time-step inputs to SubDyn, even in stand-alone mode, by
tying the driver file to an additional input file containing
time-histories of the TP motion (displacements, velocities, and
accelerations). A sample SubDyn driver input file is given in 
:numref:`sd_appendix_B`.

Users can set the **Echo** flag in this file to TRUE so that
*SubDyn\_win32.exe* echoes the contents of the driver input file (useful
for debugging errors in the driver file). The echo file has the naming
convention of **OutRootName.dvr.ech**. **OutRootName** is specified
in the SUBDYN section of the driver input file (see below).

Environmental conditions
~~~~~~~~~~~~~~~~~~~~~~~~

Set the gravity constant using the **Gravity** parameter. SubDyn
expects a magnitude, so in SI units this would be set to 9.80665
:math:`\frac{m}{s^{2}}` for standard gravity. **WtrDpth** specifies
the water depth (depth of the seabed), based on the reference MSL, and
must be a value greater than zero.


SubDyn module inputs
~~~~~~~~~~~~~~~~~~~~

**SDInputFile** is the file name of the primary SubDyn input file.
This name should be in quotations and can contain an absolute path or a
relative path. All SubDyn-generated output files will be prefixed with
**OutRootName**. If this parameter includes a file path, the output
will be generated in that folder. If this output is left empty,
the driver filename is used (without the extension) is used.
**NSteps** specifies the number of
simulation time steps, and **TimeStep** specifies the time between
steps. Next, the user must specify the location of the TP reference
point **TP\_RefPoint** (in the global reference system). This is
normally set by FAST through the ElastoDyn input file, and it is the
so-called *platform* reference point location. When coupled to FAST, the
*platform* reference point location is identified by only one (*Z*)
coordinate. The interface joints, defined in SubDyn’s main input file,
are rigidly connected to this reference point. To utilize the same
geometry definition within SubDyn’s main input file, while still
allowing for different substructure orientations about the vertical, the
user can set **SubRotateZ** to a prescribed angle in degrees with
respect to the global *Z*-axis. The entire substructure will be rotated
by that angle. (This feature is only available in stand-alone mode.)


Input motion 
~~~~~~~~~~~~

Setting **InputsMod** = 0 sets all TP reference-point input motions to
zero for all time steps. Setting **InputsMod** = 1 allows the user to
provide steady (fixed) inputs for the TP motion in the STEADY INPUTS
section of the file—\ **uTPInSteady**, **uDotTPInSteady**, and
**uDotDotTPInSteady** following the same convention as Table 1
(without time). Setting **InputsMod** = 2 allows the user to input a
time-series file whose name is specified via the **InputsFile**
parameter. The time-series input file is a text-formatted file. This
file has no header lines, **NSteps** rows, and each *i*\ :sup:`th` row
has the first column showing time as *t* = ( *i* – 1 )\*\ **TimeStep**
(the data will not be interpolated to other times). The remainder of
each row is made of white-space-separated columns of floating point
values representing the necessary motion inputs as shown in Table 1. All
motions are specified in the global, inertial-frame coordinate system.
SubDyn does not check for physical consistency between the displacement,
velocity, and acceleration motions specified for the TP reference point
in the driver file.

Table 1. TP Reference Point Inputs Time-Series Data File Contents

+-----------------+-------------------------------------------------------------------------------------------------------+------------------------------------------+
| Column Number   | Input                                                                                                 | Units                                    |
+=================+=======================================================================================================+==========================================+
| 1               | Time step value                                                                                       |  `s`                                     |
+-----------------+-------------------------------------------------------------------------------------------------------+------------------------------------------+
| 2-4             | TP reference point translational displacements along *X*, *Y*, and *Z*                                |  `m`                                     |
+-----------------+-------------------------------------------------------------------------------------------------------+------------------------------------------+
| 5-7             | TP reference point rotational displacements about *X*, *Y*, and *Z* (small angle assumptions apply)   | `rad/s`                                  |
+-----------------+-------------------------------------------------------------------------------------------------------+------------------------------------------+
| 8-10            | TP reference point translational velocities along *X*, *Y*, and *Z*                                   | `m/s`                                    |
+-----------------+-------------------------------------------------------------------------------------------------------+------------------------------------------+
| 11-13           | TP reference point rotational velocities about *X*, *Y*, and *Z*                                      | `rad/s`                                  |
+-----------------+-------------------------------------------------------------------------------------------------------+------------------------------------------+
| 14-16           | TP reference point translational accelerations along *X*, *Y*, and *Z*                                | `m/s^2`                                  |
+-----------------+-------------------------------------------------------------------------------------------------------+------------------------------------------+
| 17-19           | TP reference point rotational accelerations about *X*, *Y*, and *Z*                                   | `rad/s^2`                                |
+-----------------+-------------------------------------------------------------------------------------------------------+------------------------------------------+


Applied loads
~~~~~~~~~~~~~
The next section of the input file provides options to apply loads at given joints of the structure.
**nAppliedLoads** [-] specifies the number of applied loads listed in the subsequent table.
The user can specify a combination of steady loads and unsteady loads (both are added together).
The loads are in the global coordinate sytem.
The steady loads are given as columns of the table
(Fx, Fy, Fz, Mx, My, Mz), whereas the unsteady loads are provided in a CSV file.
The CSV filename is provided in the last entry of the table. 
If the filename is empty, the unsteady loads are not read.
An example of applied loads table is given below:

.. code::

   ---------------------- LOADS --------------------------------------------------------------------
   1    nAppliedLoads  - Number of applied loads at given nodes
   ALJointID    Fx     Fy    Fz     Mx     My     Mz   UnsteadyFile
      (-)       (N)    (N)   (N)   (Nm)   (Nm)   (Nm)     (-)
      15       100      0     0     0       0      0      ""
      23        0       0     0     0       0      0      "Force_TS.csv"

In the above example, a steady applied force of 100N is applied at the joint with ID=15 of the structure,
and an unsteady load is applied to joint 23. The time series of unsteady loads is a CSV file with
7 columns (Time, Fx, Fy, Fz, Mx, My, Mz) and one line of header. The time vector needs to be increasing, 
but does not need to be linear or cover the full range of the simulation. Interpolation is done in between
time stamps, and the first are last values are used for times smaller and larger than the simulation time range respectively.
An example of time series is shown below:

.. code::

   #Time_[s] , Fx_[N] , Fy_[N] , Fz_[N] , Mx_[Nm] , My_[Nm] , Mz_[Nm]
   0.0       , 0.0    , 0.0    , 0.0    , 0.0     , 0.0     , 0.0
   10.0      , 100.0  , 0.0    , 0.0    , 0.0     , 0.0     , 0.0
   11.0      , 0.0    , 0.0    , 0.0    , 0.0     , 0.0     , 0.0




.. _sd_main-input-file:

SubDyn Primary Input File
-------------------------
The SubDyn input file defines the substructure geometry, integration and
simulation options, finite-element parameters, and output channels. The
geometry of members is defined by joint coordinates of the undisplaced
substructure in the global reference system (inertial-frame coordinate
system), with the origin at the intersection of the undeflected tower
centerline with MSL or ground level for land-based structures. A member
connects two joints; multiple members can use a common joint. The
hydrodynamic and gravity loads are applied at the nodes, which are the
resultant of member refinement into multiple (**NDiv** input) elements
(nodes are located at the ends of each element), as calculated by the
module. Member properties include outer diameter, thickness, material
density, and Young’s and shear moduli. Member properties are specified
at the joints; if properties change from one joint to the other, they
will be linearly interpolated for the inner nodes. Unlike the geometric
properties, the material properties are not allowed to change within a
single member.

Future releases will allow for members of different cross-sections,
i.e., noncircular members. For this reason, the input file has
(currently unused) sections dedicated to the identification of direction
cosines that in the future will allow the module to identify the correct
orientation of noncircular members. The current release only accepts
tubular (circular) members.

The file is organized into several functional sections. Each section
corresponds to an aspect of the SubDyn model and substructure.

If this manual refers to an ID in a table entry, it is an integer
identifier for the table entry and must be unique for a given table
entry.

A sample SubDyn primary input file is given in :numref:`sd_appendix_A`.

The input file begins with two lines of header information, which is for
the user but is not used by the software.


Simulation Control Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Users can set the **Echo** flag to TRUE to have SubDyn echo the
contents of the SubDyn input file (useful for debugging errors in the
input file). The echo file has the naming convention of
**OutRootName.SD.ech**. **OutRootName** is either specified in the
SUBDYN section of the driver input file when running SubDyn standalone,
or by FAST, when running a coupled simulation, from FAST’s main input
file.

**SDdeltaT** specifies the fixed time step of the integration in
seconds. The keyword ‘DEFAULT’ may be used to indicate that the module
should employ the time step prescribed by the driver code
(FAST/standalone driver program).

**IntMethod** specifies the integration algorithm to use. There are
four options: 1) Runge-Kutta 4\ :sup:`th`-order explicit (RK4); 2)
Adams-Bashforth 4\ :sup:`th`-order explicit predictor (AB4); 3)
Adams-Bashforth-Moulton 4\ :sup:`th`-order explicit predictor-corrector
(ABM4); 4) Adams-Moulton implicit 2\ :sup:`nd`-order (AM2). See Section
on how to properly select this and the previous parameter values.

**SttcSolve** is a flag that specifies whether the static improvement method 
(SIM, see :numref:`SD_SIM`)
shall be employed. Through this method, all (higher frequency) modes
that are not considered by the C-B reduction are treated
quasi-statically. This treatment helps
minimize the number of retained modes needed to capture effects such as
static gravity and buoyancy loads, and high-frequency loads transferred
from the turbine. Recommended to set to True.


**GuyanLoadCorrection** is a flag to specify whether the extra moment due to 
the lever arm from the Guyan deflection of the structure is to be added to the loads
passed to SubDyn, and, whether the FEM representation should be expressed in the rotating 
frame in the floating case (the rotation is induced by the rigid body Guyan modes).
See section :numref:`SD_Loads` for details. Recommended to set to True.


FEA and Craig-Bampton Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**FEMMod** specifies one of the following options for finite-element
formulation: 1) Euler-Bernoulli; 3) Timoshenko. Tapered formulations (2
and 4) have yet to be implemented and will be available in a future
release.

**NDiv** specifies the number of elements per member. Analysis nodes
are located at the ends of elements and the number of analysis nodes per
member equals **NDiv** + 1. **NDiv** is applied uniformly to all
members regardless of the member’s length, hence it could result in
small elements in some members and long elements in other members.
Increasing the number of elements per member may increase accuracy, with
the trade-off of increased memory usage and computation time. We
recommend using **NDiv** > 1 when modeling tapered members.

**CBMod** is a flag that specifies whether or not the C-B reduction
should be carried out by the module. If FALSE, then the full
finite-element model is retained and **Nmodes** is ignored.

**Nmodes** sets the number of internal C-B modal DOFs to retain in the
C-B reduction. **Nmodes** = 0 corresponds to a Guyan (static)
reduction. **Nmodes** is ignored if **CBMod** is set to FALSE,
meaning the full finite-element model is retained by keeping all modes
(i.e. a modal analysis is still done, and all the modes are used as DOFs)  .


**JDampings** specifies value(s) of damping coefficients as a
percentage of critical damping for the retained C-B modes. Distinct
damping coefficients for each retained mode should be listed on the same
line, separated by white space. If the number of **JDampings** is less
than the number of retained modes, the last value will be replicated for
all the remaining modes. (see :numref:`SD_DampingSpecifications`)

**GuyanDampMod** Guyan damping [0=none, 1=Rayleigh Damping, 2= user specified 6x6 matrix] (see :numref:`SD_DampingSpecifications`)


**RayleighDamp** Mass and stiffness proportional damping  coefficients (:math:`(\alpha,\beta)` Rayleigh damping) [only if GuyanDampMod=1]
Guyan damping matrix (6x6) [only if GuyanDamgMod=2] (see :numref:`SD_DampingSpecifications`)


**Guyan damping matrix**:
The 6 lines following this input line consits of the 6x6 coefficients of the damping matrix to be applied at the interface. (see :numref:`SD_DampingSpecifications`)


For more information on these parameters and guidelines on how to set
them, see Sections :numref:`sd_modeling-considerations` and :numref:`subdyn-theory`.

Structure Joints
~~~~~~~~~~~~~~~~

The finite-element model is based on a substructure composed of joints
interconnected by members. **NJoints** is the user-specified number of
joints, and determines the number of rows in the subsequent table.
Because a member connects two joints, **NJoints** must be greater than
or equal to two. Each joint listed in the table is identified by a
unique integer, **JointID**; each integer between one and
**NJoints** must be present in the table, but they need not be
sequential. The (*X*,\ *Y*,\ *Z*) coordinate of each joint is specified
in the substructure (SS) coordinate system, which coincides with the
global inertial-frame coordinate system via **JointXss**,
**JointYss**, and **JointZss**, respectively. This version of SubDyn
does not consider overlap when multiple members meet at a common joint,
therefore, it tends to overestimate the total substructure mass. Member
overlap and node offset calculations will be considered in a future
release of SubDyn.
The fifth column specifies the **JointType** (see :numref:`SD_FEM`):

- Cantilever joints (*JointType=1*)

- Universal joint (*JointType=2*)

- Pin joint (*JointType=3*)

- Ball joint (*JointType=4*)

The three following columns specify the vector coordinates of the direction around which rotation is free for a pin joints.
The last column, **JointStiff** specify a value of additional stiffness to be added to the "free" rotational DOFs of Ball, Pin and Universal joints.


Note for HydroDyn coupling: modeling a fixed-bottom substructure
embedded into the seabed (e.g., through piles or suction buckets)
requires that the lowest member joint(s) in HydroDyn lie(s) below the
water depth. Placing a joint at or above the water depth will result in
static and dynamic pressure loads being applied at the joint. When
SubDyn is coupled to FAST, the joints and members need not match between
HydroDyn and SubDyn—FAST’s mesh-mapping utility handles transfer of
motion and loads across meshes in a physically relevant manner (Sprague
et al. 2014), but consistency between the joints and members in HydroDyn
and SubDyn is advised.   


An example of joint table is given below

.. code::

    3   NJoints  - Number of joints (-)
    JointID JointXss JointYss  JointZss JointType JointDirX JointDirY JointDirZ JointStiff 
      (-)      (m)      (m)       (m)     (-)        (-)       (-)       (-)     (Nm/rad) 
      101      0.0      0.0      50.0      1         0.0       0.0       0.0       0.0    
      111      0.0      0.0      10.0      2         0.0       1.0       0.0     100.0    
      102      0.0      0.0     -45.0      1         0.0       0.0       0.0       0.0    


Base Reaction Joints
~~~~~~~~~~~~~~~~~~~~~

SubDyn requires the user to specify the boundary joints. **NReact**
should be set equal to the number of joints (defined earlier) at the
bottom of the structure (i.e., seabed) that are fully constrained;
**NReact** also determines the number of rows in the subsequent table.
In SubDyn, **NReact** must be greater than or equal to one. Each joint
listed in the table is identified by a unique integer, **RJointID**,
which must correspond to the **JointID** value found in the STRUCTURE
JOINTS table. The flags **RctTDXss**, **RctTDYss**, **RctTDZss**,
**RctRDXss**, **RctRDYss**, **RctRDZss** indicate the fixity value
for the three translations (TD) and three rotations (RD) in the SS
coordinate system (global inertial-frame coordinate system). One denotes
fixed and zero denotes free (instead of TRUE/FALSE). **SSIfile**
points to the relative path and filename for an SSI information file.
This version of SubDyn can, in fact, handle partially restrained joints
by setting one or more DOF flags to 0 and providing the appropriate
stiffness and mass matrix elements for that DOF via the **SSIfile**.
If a DOF flag is set to 1, then the node DOF is considered restrained
and the associated matrix elements potentially provided in the
**SSIfile** will be ignored.


An example of base reaction and interface table is given below

.. code::

    ------------------- BASE REACTION JOINTS
      1   NReact      - Number of Joints with reaction forces
    RJointID RctTDXss RctTDYss RctTDZss RctRDXss RctRDYss RctRDZss  SSIfile
      (-)     (flag)   (flag)   (flag)   (flag)   (flag)   (flag)   (string)
      61         1        1        1        1        1        1	    "SSI.txt"
    ------------------- INTERFACE JOINTS
      1   NInterf     - Number of interface joints locked to the Transition Piece (TP)
    IJointID ItfTDXss ItfTDYss ItfTDZss ItfRDXss ItfRDYss ItfRDZss 
      (-)     (flag)   (flag)   (flag)   (flag)   (flag)   (flag)
      24         1        1        1        1        1        1


Interface Joints
~~~~~~~~~~~~~~~~

SubDyn requires the user to specify the interface joints. **NInterf**
should be set equal to the number of joints at the top of the structure
(i.e., TP); **NInterf** also determines the number of rows in the
subsequent table. In SubDyn, **NInterf** must be greater than or equal
to one. Note that these joints will be assumed to be rigidly connected
to the platform reference point of ElastoDyn (see FAST documentation)
when coupled to FAST, or to the TP reference point if SubDyn is run in
stand-alone mode. Each joint listed in the table is identified by a
unique integer, **IJointID**, which must correspond to the *JointID*
value found in the STRUCTURE JOINTS table. The flags **ItfTDXss**,
**ItfTDYss**, **ItfTDZss**, **ItfRDXss**, **ItfRDYss**,
**ItfRDZss** indicate the fixity value for the three translations (TD)
and three rotations (RD) in the SS coordinate system (global
inertial-frame coordinate system). One denotes fixed and zero denotes
free (instead of TRUE/FALSE). This version of SubDyn cannot handle
partially restrained joints, so all flags must be set to one; different
degrees of fixity will be considered in a future release.

Members
~~~~~~~

**NMembers** is the user-specified number of members and determines
the number of rows in the subsequent table. Each member listed in the
table is identified by a unique integer, **MemberID**. Each integer
between one and **NMembers** must be present in the table, but they
need not be sequential. For each member distinguished by **MemberID**,
**MJointID1** specifies the starting joint and **MJointID2**
specifies the ending joint, corresponding to an identifier
(**JointID**) from the STRUCTURE JOINTS table. Likewise,
**MPropSetID1** corresponds to the identifier **PropSetID** from the
MEMBER X-SECTION PROPERTY table (discussed next) for starting
cross-section properties and **MPropSetID2** specifies the identifier
for ending cross-section properties, allowing for tapered members.
The sixth column specify the member type  **MType**.
A member is one of the three following types (see :numref:`SD_FEM`):

- Beams (*MType=1*), Euler-Bernoulli (*FEMMod=1*) or Timoshenko (*FEMMod=3*)

- Pretension cables (*MType=2*)

- Rigid link (*MType=3*)

**COSMID** refers to the IDs of the members' cosine matrices for noncircular
members; the current release uses SubDyn's default direction cosine convention
if it's not present or when COSMID values are -1.


An example of member table is given below

.. code::

     2   NMembers    - Number of frame members
  MemberID   MJointID1   MJointID2   MPropSetID1   MPropSetID2  MType   COSMID
    (-)         (-)         (-)          (-)           (-)        (-)      (-)
     10        101         102            2             2          1
     11        102         103            2             2          1




Member Cross-Section Properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Members in SubDyn are assumed to be straight, circular, possibly
tapered, and hollow cylinders. Future releases will allow for generic
cross-sections to be employed. These special cross-section members will
be defined in the second of two tables in the input file (Member
X-Section Property data 2/2), which is currently ignored.

For the circular cross-section members, properties needed by SubDyn are
material Young’s modulus, **YoungE**, shear modulus, **ShearG**, and
density, **MatDens**, member outer diameter, **XsecD**, and member
thickness, **XsecT**. Users will need to create an entry in the first
table within this section of the input file distinguished by
**PropSetID**, for each unique combination of these five properties.
The member property-set table contains **NPropSets** rows. The member
property sets are referred to by their **PropSetID** in the MEMBERS
table, as described in Section . Note, however, that although diameter
and thickness will be linearly interpolated within an individual member,
SubDyn will not allow *material* properties to change within an
individual member.

The second table in this section of the input file (not to be used in
this release) will have **NXPropSets** rows (assumed to be zero for
this release), and have additional entries when compared to the previous
table, including: cross-sectional area (**XsecA**), cross-sectional
shear area along the local principal axes *x* and *y* (**XsecAsx**,
**XsecAsy**), cross-sectional area second moment of inertia about *x*
and *y* (**XsecJxx**, **XsecJyy**), and cross-sectional area polar
moment of inertia (**XsecJ0**). The member cosine matrix section (see
Section ) will help determine the correct orientation of the members
within the assembly.





Cable Properties
~~~~~~~~~~~~~~~~


Members that are specified as pretension cables (**MType=2**), 
have their properties defined in the cable properties table. 
The table lists for each cable property: the property ID (**PropSetID**), the cable tension stiffness (**EA**), 
the material density (**MatDens**), the pretension force (**T0**), and the control channel (**CtrlChannel**).
The control channel is only used if ServoDyn provides dedicated control signal, in which case
the cable tension (given in terms of a length change :math:`\Delta l`) 
is dynamically changed (see :numref:`SD_ControlCable`).
The FEM representation of pretension cable is given in :numref:`SD_PretensionCable`.

An example of cable properties table is given below:

.. code::

    -------------------------- CABLE PROPERTIES  -------------------------------------
                 2   NCablePropSets   - Number of cable cable properties
    PropSetID   EA     MatDens    T0    CtrlChannel
      (-)      (N)     (kg/m)    (N)      (-)
       11      210E7   7850.0    2E7       1 
       10      210E7   7850.0    1E7       0 


Rigid link Properties
~~~~~~~~~~~~~~~~~~~~~

Members that are specified as rigid links (**MType=3**), 
have their properties defined in the rigid link properties table. 
The table lists the material density (**MatDens**) for each rigid link property.
The FEM representation of rigid links is given in :numref:`SD_RigidLinks`.

An example of rigid link properties table is given below

.. code::

   ----------------------- RIGID LINK PROPERTIES ------------------------------------
                1   NRigidPropSets - Number of rigid link properties
   PropSetID   MatDens   
     (-)       (kg/m)    
      12       7850.0
       3       7000.0
       
















Member Cosine Matrices COSM (i,j)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**NCOSMs** rows, one for each unique member orientation set, will need
to be provided. Each row of the table will list the nine entries of the
direction cosine matrices (COSM11, COSM12,…COSM33) for matrix elements.
Each row is a vector in the global coordinate system for principal axes 
in the x, y and z directions respectively. These vectors need to be 
specified with an extremely high level of precision for results to be
equivalent to an internal calculation.

Joint Additional Concentrated Masses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SubDyn can accept **NCmass** lumped masses/inertias defined at the
joints. The subsequent table will have **NCmass** rows, in which for
each joint distinguished by **CMJointID** (corresponding to an
identifier, **JointID**, from the STRUCTURE JOINTS table), **JMass**
specifies the lumped mass value, and **JMXX**, **JMYY**, **JMZZ**
specify the mass second moments of inertia with respect to the SS
coordinate system (not the element system).
Latest version of SubDyn accept 6 additional columns 
(**JMXY**, **JMXZ**, **JMYZ**, **MCGX**, **MCGY**, **MCGZ**) 
to specify off-diagonal terms.

The additional mass matrix added to the node is computed in the SS system as follows:

.. math::

      M_\text{add}=
      \begin{bmatrix}
      m    & 0    & 0    & 0                    & z m                    & -y m          \\
      0    & m    & 0    & -z m                 & 0                      & x m           \\
      0    & 0    & m    & y m                  & -x m                   & 0             \\
      0    & -z m & y m  & J_{xx} + m (y^2+z^2) & J_{xy} - m x y         & J_{xz}  - m x z  \\
      z m  & 0    & -x m & J_{xy} - m x y       & J_{yy} + m (x^2+z^2)   & J_{yz}  - m y z  \\
      -y m & x m  & 0    & J_{xz} - m x z       & J_{yz} - m y z         & J_{zz}  + m (x^2+y^2)\\
      \end{bmatrix}

with :math:`m` the parameter **JMass**, and :math:`x,y,z`, the CG offsets.


An example of concentrated mass table is given below:

.. code::

         2  NCmass - Number of joints with concentrated masses; (SS coord system)
    CMJointID  JMass    JMXX    JMYY    JMZZ   JMXY     JMXZ   JMYZ   MCGX  MCGY MCGZ   
      (-)       (kg)  (kgm^2) (kgm^2) (kgm^2) (kgm^2) (kgm^2) (kgm^2)  (m)  (m)  (m)
       1        4090     0       0       0       0        0       0      0    0    0
       3        4.2e6    0       0     3.3e9     0        0       0      0    0    0


Output: Summary and Outfile
~~~~~~~~~~~~~~~~~~~~~~~~~~~
In this section of the input file, the user sets flags and switches for
the desired output behavior.

Specifying **SumPrint** = TRUE causes SubDyn to generate a summary file
with name **OutRootName**.SD.sum*. **OutRootName** is either
specified in the SUBDYN section of the driver input file when running
SubDyn in stand-alone mode, or in the FAST input file when running a
coupled simulation. See Section 4.2 for summary file details.

The following two inputs specified whether mode shapes should be written
to disk.  **OutCBModes** is a flag that controls the output of the Guyan
and Craig-Bampton modes. Similarly, **OutFEMModes**, controls the output
of the FEM modes (full sytem with constraints prior to the CB-reduction).
For now, only the first 30 FEM modes are written to disk, but all CB modes
selected by the users are written. 
For both inputs, the following options are available: 0, no ouput, 1, outputs
in JSON format. The JSON files contain nodes coordinates, connectivity between the nodes, 
displacements for each modes and nodes, and frequencies for each modes.
The reading of these files should be straightforward using Matlab or Python using a JSON format parser. 
The files can be opened to visualize the modes using the tool viz3danim
(see the `live version <https://ebranlard.github.io/viz3Danim/>`_
, or its `github repository <https://github.com/ebranlard/viz3danim>`_).

Currently, **OutCOSM** is ignored. In future releases,
specifying **OutCOSM** = TRUE will cause SubDyn to include direction
cosine matrices (undeflected) in the summary file for only those members
requested in the list of output channels.

Specifying **OutAll** = TRUE causes SubDyn to output forces and
moments at all of the joints (not internal nodes). That is, the static
(elastic) and dynamic (inertia) components of the three forces and three
moments at the end node of each member connected to a given joint are
output for all joints. These outputs are included within the
**OutRootName**.SD.out* output file in addition to those directly
specified through the output channels section below.

If **OutSwtch** is set to one, outputs are sent to a file with the
name **OutRootName**.SD.out*. If **OutSwtch** is set to two, outputs
are sent to the calling program (FAST) for writing in its main output
file (not available in stand-alone mode). If **OutSwtch** is set to
three, both file outputs occur. In stand-alone mode, setting
**OutSwtch** to two results in no output file being produced.

If **TabDelim** is set to TRUE and **OutSwtch** is set to one, the
output file **OutRootName**.SD.out* will be tab-delimited.

With **OutDec** set to an integer value greater than one, the output
file data rate will be decimated, and only every **OutDec**-th value
will be written to the file. This applies only to SubDyn’s output file
(**OutRootName**.SD.out*)—not FAST’s.

The **OutFmt** and **OutSFmt** parameters control the formatting of
SubDyn’s output file for the output data and the channel headers,
respectively. SubDyn currently does not check the validity of these
format strings. They need to be valid Fortran format strings.
**OutSFmt** is used for the column header and **OutFmt** is used for
the channel data. Therefore, in order for the headers and channel data
to align properly, the width specification should match. For example:

| "ES11.4" OutFmt
| "A11" OutSFmt.


.. _SD_Member_Output:

Member Output List
~~~~~~~~~~~~~~~~~~

SubDyn can output load and kinematic quantities at up to nine locations
for up to nine different members, for a total of 81 possible local
member output locations. **NMOutputs** specifies the number of members
that output is requested for. The user must create a table entry for
each requested member. Within a row of this table, **MemberID** is the
ID specified in the MEMBERS table, and **NOutCnt** specifies how many
nodes along the member will generate output. **NodeCnt** specifies
those node numbers (a separate entry on the same line for each node) for
output as an integer index from the start-joint (node 1) to the
end-joint (node **NDiv** + 1) of the member. The outputs specified in
the SDOutList section determines which quantities are actually output at
these locations.

Output Channels- SDOutList Section
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section specifies which quantities are output by SubDyn. Enter one
or more lines containing quoted strings that in turn contain one or more
output parameter names. Separate output parameter names by any
combination of commas, semicolons, spaces, and/or tabs. If a parameter
name is prefixed with a minus sign, “-”, underscore, “\_”, or the
characters “m” or “M”, SubDyn will multiply the value for that channel
by –1 before writing the data. The parameters are written in the order
they are listed in the input file. SubDyn allows the use of multiple
lines so that users can break their lists into meaningful groups and so
the lines can be shorter. Comments may also be entered after the closing
quote on any of the lines. Entering a line with the string “END” at the
beginning of the line or at the beginning of a quoted string found at
the beginning of the line will cause SubDyn to quit scanning for more
lines of channel names. Modal kinematics and member-node-, base-, and
interface-related kinematic and load quantities can be selected.
Member-node-related data follow the organization described in Section .
If SubDyn encounters an unknown/invalid channel name, it prints an error
message and halts execution. Please refer to :numref:`sd_appendix_C` for a complete
list of possible output parameters and their names.

.. _sd_ssi_inputfile:

SSI Input File
--------------

Individual SSI files (*SSIfiles*) can be provided for each restrained
node, therefore the maximum number of SSIfiles is **NReact**. In an
SSIfile, up to 21 elements for the SSI mass matrix and up to 21 SSI
stiffness matrix elements can be provided. The mass and stiffness
elements account for both pile and soil effects. No additional damping
can be provided at this point.

The order of the elements is not important, because each element value
is accompanied by a string label that identifies the actual element. The
stiffness matrix accepted labels are: 'Kxx', 'Kxy', 'Kyy', 'Kxz', 'Kyz’,
'Kzz’, 'Kxtx', 'Kytx', 'Kztx', 'Ktxtx', 'Kxty', 'Kyty','Kzty’, 'Ktxty',
'Ktyty', ‘Kxtz', 'Kytz', 'Kztz', 'Ktxtz', 'Ktytz', 'Ktztz'.

If any matrix element is not provided it will be set to infinity (i.e.,
machine ‘huge’) by default.

For the mass matrix the accepted labels are:
'Mxx','Mxy','Myy','Mxz','Myz', 'Mzz','Mxtx','Mytx','Mztx', 'Mtxtx',
'Mxty', 'Myty', 'Mzty', 'Mtxty', 'Mtyty', 'Mxtz', 'Mytz', 'Mztz',
'Mtxtz', 'Mtytz', 'Mtztz'. If any matrix element is not provided it will
be set to 0 by default. The labels contain ‘K’ or ‘M’ to specify
stiffness or mass matrix elements, and then the directions they apply
to, e.g., ‘Kxy’ refers to the force along x due to a unit displacement
along y; the ‘t’ refers to the rotation about one of the ‘x’,’y’, or ’z’
axes in the global coordinate system.

Units are in SI system (N/m; N/m/rad; Nm/rad, Kg, kgm, kgm2).

Note that by selecting fixities of 1 in the various DOFs of the
restrained nodes, the columns and rows associated with those DOFs will
be removed, therefore the associated matrix elements will be ignored.

A sample SubDyn SSI input file is given in :numref:`sd_appendix_C`.
