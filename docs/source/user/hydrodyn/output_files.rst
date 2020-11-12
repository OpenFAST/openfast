.. _hd-output:

Output Files
============
HydroDyn produces four types of output files: an echo file, a
wave-elevations file, a summary file, and a time-series results file.
The following sections detail the purpose and contents of these files.

Echo Files
~~~~~~~~~~
If you set the **Echo** flag to TRUE in the HydroDyn driver file or the
HydroDyn primary input file, the contents of those files will be echoed
to a file with the naming conventions, **OutRootName**\ *.dvr.ech* for
the driver input file and **OutRootName**\ *.HD.ech* for the HydroDyn
primary input file. **OutRootName** is either specified in the HYDRODYN
section of the driver input file, or by the FAST program. The echo files
are helpful for debugging your input files. The contents of an echo file
will be truncated if HydroDyn encounters an error while parsing an input
file. The error usually corresponds to the line after the last
successfully echoed line.

Wave-Elevations File
~~~~~~~~~~~~~~~~~~~~
Setting **WaveElevSeriesFlag** in the driver file to TRUE enables the
outputting of a grid of wave elevations to a text-based file with the
name ``OutRootName.WaveElev.out``. The grid consists of
**WaveElevNX** by **WaveElevNY** wave elevations (centered at *X* = 0,
*Y* = 0) with a **dX** and **dY** spacing in the global inertial-frame
coordinate system. These wave elevations are distinct and output
separately from the wave elevations determined by **NWaveElev** in the
HydroDyn primary input file, such that the total number of wave
elevation outputs is **NWaveElev** + ( **WaveElevNX** × **WaveElevNY**
). The wave-elevation output file ``OutRootName.WaveElev.out``
contains the total wave elevation, which is the sum of the first- and
second-order terms (when the second-order wave kinematics are optionally
enabled).

.. _hd-summary-file:

Summary File
~~~~~~~~~~~~
HydroDyn generates a summary file with the naming convention,
**OutRootName**\ *.HD.sum* if the **HDSum** parameter is set to TRUE.
This file summarizes key information about your hydrodynamics model,
including buoyancy, substructure volumes, marine growth weight, the
simulation mesh and its properties, first-order wave frequency
components, and the radiation kernel.

When the text refers to an index, it is referring to a given row in a
table. The indexing starts at 1 and increases consecutively down the
rows.

WAMIT-model volume and buoyancy information
-------------------------------------------
This section summarizes the buoyancy of the potential-flow-model
platform in its undisplaced configuration. For a hybrid
potential-flow/strip-theory model, these buoyancy values must be added
to any strip-theory member buoyancy reported in the subsequent sections
to obtain the total buoyancy of the platform.

Substructure Volume Calculations
--------------------------------
This section contains a summary of the total substructure volume, the
submerged volume, volume of any marine growth, and fluid-filled
(flooded/ballasted) volume for the substructure in its undisplaced
configuration. Except for the fluid-filled volume value, the reported
volumes are only for members that have the **PropPot** flag set to
FALSE. The flooded/ballasted volume applies to any fluid-filled member,
regardless of its **PropPot** flag.

Integrated Buoyancy Loads
-------------------------
This section details the buoyancy loads of the undisplaced substructure
when summed about the WRP (0,0,0). The external buoyancy includes the
effects of marine growth, and only applies to members whose **PropPot**
flag is set to FALSE. The internal buoyancy is the negative effect on
buoyancy due to flooding or ballasting and is independent of the
**PropPot** flag.

Integrated Marine Growth Weights
--------------------------------
This section details the marine growth weight loads of the undisplaced
substructure when summed about the WRP (0,0,0).

Simulation Node Table
---------------------
This table details the undisplaced nodal information and properties for
all internal analysis nodes used by the HydroDyn model. The node index
is provided in the first column. The second column maps the node to the
input joint index (not to be confused with the **JointID**). If a value
of -1 is found in this column, the node is an interior node and results
from an input member being split somewhere along its length due to the
requirements of the **MDivSize** parameter in the primary input file
members table.
The third column indicates if this node is part of a Super
Member (**JointOvrlp** = 1). The next column tells you the corresponding
input member index (not to be confused with the **MemberID**). **Nxi**,
**Nyi**, and **Nzi**, provide the (*X*,\ *Y*,\ *Z*) coordinates in the
global inertial-frame coordinate system. **InpMbrDist** provides the
normalized distance to the node from the start of the input member.
**R** is the outer radius of the member at the node (excluding marine
growth), and **t** is the member wall thickness at the node. **dRdZ** is
the taper of the member at the node, **tMG** is the marine growth
thickness, and **MGDens** is the marine growth density. **PropPot**
indicates whether the element attached to this node is modeled using
potential-flow theory. If **FilledFlag** is TRUE, then **FillDens**
gives the filled fluid density and **FillFSLoc** indicates the
free-surface height (*Z*-coordinate). **Cd**, **Ca**, **Cp**, **AxCa**,
**AxCp**, **JAxCd**, **JAxCa**, and **JAxCp** are the viscous-drag,
added-mass, dynamic-pressure, axial added-mass, axial dynamic-pressure,
end-effect axial viscous-drag, end-effect axial added-mass, and
end-effect axial dynamic-pressure coefficients, respectively. **NConn**
gives the number of elements connected to node, and **Connection List**
is the list of element indexes attached to the node.

.. TODO 7.5.2 is the theory section which does not yet exist.
.. See Section 7.5.2 for the member splitting rules used by HydroDyn.

Simulation Element Table
------------------------
This section details the undisplaced simulation elements and their
associated properties. A suffix of 1 or 2 in a column heading refers to
the element’s starting or ending node, respectively. The first column is
the element index. **node1** and **node2** refer to the node index found
in the node table of the previous section. Next are the element
**Length** and exterior **Volume**. This exterior volume calculation
includes any effects of marine growth. **MGVolume** provides the volume
contribution due to marine growth. The cross-sectional properties of
outer radius (excluding marine growth), marine growth thickness, and
wall thickness for each node are given by **R1**, **tMG1**, **t1**,
**R2**, **tMG2**, and **t2**, respectively. **MGDens1** and **MGDens2**
are the marine growth density at node 1 and 2. **PropPot** indicates if
the element is modeled using potential-flow theory. If the element is
fluid-filled (has flooding or ballasting), **FilledFlag** is set to
**T** for TRUE. **FillDensity** and **FillFSLoc** are the filled fluid
density and the free-surface location’s *Z*-coordinate in the global
inertial-frame coordinate system. **FillMass** is calculated by
multiplying the **FillDensity** value by the element’s interior volume.
Finally, the element hydrodynamic coefficients are listed. These are the
same coefficients listed in the node table (above).

Summary of User-Requested Outputs
---------------------------------
The summary file includes information about all requested member and
joint output channels.

Member Outputs
++++++++++++++
The first column lists the data channel’s string label, as entered in
the OUTPUT CHANNELS section of the HydroDyn input file. **Xi**, **Yi**,
**Zi**, provide the output’s undisplaced spatial location in the global
inertial-frame coordinate system. The next column, **InpMbrIndx**, tells
you the corresponding input member index (not to be confused with the
**MemberID**). Next are the coordinates of the starting (**StartXi**,
**StartYi**, **StartZi**) and ending (**EndXi**, **EndYi**, **EndZi**)
nodes of the element containing this output location. **Loc** is the
normalized distance from the starting node of this element.

Joint Outputs
+++++++++++++
The first column lists the data channel’s string label, as entered in
the OUTPUT CHANNELS section of the HydroDyn input file. **Xi**, **Yi**,
**Zi**, provide the output’s undisplaced spatial location in the global
inertial-frame coordinate system. **InpJointID** specifies the
**JointID** for the output as given in the MEMBER JOINTS table of the
HydroDyn input file.

The Wave Number and Complex Values of the Wave Elevations as a Function of Frequency
------------------------------------------------------------------------------------
This section provides the frequency-domain description (in terms of a
Discrete Fourier Transform or DFT) of the first-order wave elevation at
(0,0) on the free surface, but is not written when **WaveMod** = 0 or 6.
The first column, **m**, identifies the index of each wave frequency
component. The finite-depth wave number, frequency, and direction of the
wave component are given by **k**, **Omega**, and **Direction**,
respectively. The last two columns provide the real
(**REAL(DFT{WaveElev})**) and imaginary (**IMAG(DFT{WaveElev})**)
components of the DFT of the first-order wave elevation. The DFT
produces includes both the negative- and positive-frequency components.
The negative-frequency components are complex conjugates of the positive
frequency components because the time-domain wave elevation is
real-valued. The relationships between the negative- and
positive-frequency components of the DFT are given by
:math:`k\left( - \omega \right) = - k\left( \omega \right)` and
:math:`H\left( - \omega \right) = {H\left( \omega \right)}^{*}`, where
*H* is the DFT of the wave elevation and *\** denotes the complex
conjugate.

Radiation Memory Effect Convolution Kernel
------------------------------------------

In the potential-flow solution based on frequency-to-time-domain
transforms, HydroDyn computes the radiation kernel used by the
convolution method for calculating the radiation memory effect through
the cosine transform of the 6x6 frequency-dependent hydrodynamic damping
matrix from the radiation problem. The resulting time-domain radiation
kernel (radiation impulse-response function)—which is a 6x6
time-dependent matrix—is provided in this section. **n** and **t** give
the time-step index and time, which are followed by the elements
(**K11**, **K12**, etc.) of the radiation kernel associated with that
time. Because the frequency-dependent hydrodynamic damping matrix is
symmetric, so is the radiation kernel; thus, only the diagonal and
upper-triangular portion of the matrix are provided. The radiation
kernel should decay to zero after a short amount of time, which should
aid in selecting an appropriate value of **RdtnTMax**.

Results File
~~~~~~~~~~~~

The HydroDyn time-series results are written to a text-based file with
the naming convention ``OutRootName.HD.out`` when **OutSwtch** is
set to either 1 or 3. If HydroDyn is coupled to FAST and **OutSwtch** is
set to 2 or 3, then FAST will generate a master results file that
includes the HydroDyn results. The results are in table format, where
each column is a data channel (the first column always being the
simulation time), and each row corresponds to a simulation output time
step. The data channels are specified in the OUTPUT CHANNELS section of
the HydroDyn primary input file. The column format of the
HydroDyn-generated file is specified using the **OutFmt** and
**OutSFmt** parameter of the primary input file.
