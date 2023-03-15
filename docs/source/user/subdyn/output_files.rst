.. _sd_output-files:
   
Output Files
============

SubDyn produces three types of output files: an echo file, a summary
file, and a time-series results file. The following sections detail the
purpose and contents of these files.

Echo File
---------

If the user sets the **Echo** flag to TRUE in the SubDyn driver file
or the primary SubDyn input file, the contents of those files will be
echoed to a file with the naming conventions, **OutRootName.dvr.ech**
for the driver input file and **OutRootName.SD.ech** for the primary
SubDyn input file. **OutRootName** is either specified in the SUBDYN
section of the driver input file, or in the FAST input file. The echo
files are helpful for debugging the input files. The contents of an echo
file will be truncated if SubDyn encounters an error while parsing an
input file. The error usually corresponds to the line after the last
successfully echoed line.

.. _sd_sum-file:

Summary File
------------

SubDyn generates a summary file with the naming convention,
**OutRootName.SD.sum** if the **SumPrint** parameter is set to TRUE.
This file summarizes key information about the substructure model,
including:

-  Undisplaced node geometry: a list of all of the (**NNodes**) nodes
   and the *X*,\ *Y*,\ *Z* coordinates in the global SS coordinate
   system. Note that **NNodes** may be greater or equal to
   **NJoints**, depending on **NDiv** (primary input file
   parameters).

-  Element connectivity and properties at end nodes: a list of all
   (**NElems**) elements, the start and end nodes (**Node\_I**,
   **Node\_J**) and the ID of the property set (**Prop\_I**,
   **Prop\_J**) at the start and end nodes. **NElems** may be
   greater or equal to **NMembers**, depending on **NDiv** (primary
   input file parameters).

-  Property sets. If tapered members are used, additional property sets
   may be included beyond those specified in the main input file, based
   on interpolated diameter and thickness values. Headers and their
   meanings are identical to those described in Section .

-  Reaction DOFs and interface DOFs and their associated fixity; the
   actual indices of the DOFs (**DOF\_ID**) associated with reaction
   and interface nodes are listed together with the (1/0) flag to
   distinguish the fixity level.

-  Concentrated mass schedule. This is an echo of the equivalent section
   in the primary input file. Refer to Section .

-  Member schedule including connectivity to joints, nodes, and their
   masses. A table lists all of the members by identifier
   (**MemberID**), with their start and end nodes (**Joint1\_ID**,
   **Joint2\_ID**), associated mass (**Mass**), and list of node
   identifiers along the length of the members.

-  Direction cosine matrices for the members. Each row (columns 2-10)
   corresponds to the direction cosine matrix entries (**DC(1,1)**
   through **DC(3,3)**) for the member whose identifier is listed in
   the first column. The direction cosine matrices specify the
   transformation from the global reference to the local coordinate
   system for each member.

-  Sorted eigenfrequencies [in Hertz (Hz)] for the full substructural
   system (neglecting a possible coupling to ElastoDyn through FAST),
   assuming the TP reference point is a free end. There are a total of
   **NDOFs** eigenfrequencies and eigenvectors.

-  Sorted eigenfrequencies (in Hz) for the C-B reduced system, assuming
   the TP reference point is a fixed end. There are a total of
   **Nmodes** C-B reduced eigenfrequencies and eigenvectors.

-  Full substructural system eigenvectors. Each column represents an
   eigenvector associated with the corresponding eigenfrequency
   identified previously in the file.

-  C-B reduced system eigenvectors (**PhiM** matrix). Each column
   represents an eigenvector associated with the corresponding
   eigenfrequency identified previously in the file.

-  **PhiR** matrix or displacements of the internal nodes caused by
   unit rigid body motions of the interface DOFs (see Section ). Each
   column of the matrix represents the internal DOF displacements for a
   given unit rigid-body motion along an interface DOF for each base and
   interface joint.

-  Substructure equivalent stiffness and mass matrices referred to the
   TP reference point (**KBBt** and **MBBt**), based on a Guyan
   reduction. These are useful to calculate effects of substructure
   flexibility while calculating tower eigenmodes for ElastoDyn.

-  Rigid-body-equivalent mass matrix relative to global origin
   (**MRB**); a 6x6 mass matrix.

-  Substructure total (dry) mass.

-  Substructure center of mass coordinates in the global coordinate
   system.

The various sections of the summary file and variables are
self-explanatory and easily identifiable in the file.

Results File
------------

The SubDyn time-series results are written to a text-based file with the
naming convention **OutRootName.SD.out** when **OutSwtch** is set to
either one or three. If SubDyn is coupled to FAST and **OutSwtch** is
set to two or three, then FAST will generate a master results file that
includes the SubDyn results. The results in **OutRootName.SD.out** are
in table format, where each column is a data channel (the first column
always being the simulation time), and each row corresponds to a
simulation time step. The data channels are specified in the SDOutList
section of the input file. The column format of the SubDyn-generated
file is specified using the **OutFmt** and **OutSFmt** parameters of
the input file.

