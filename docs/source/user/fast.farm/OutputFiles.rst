.. _FF:Output:

Output Files
============

FAST.Farm produces five types of output files: an echo file, a summary
file, visualization output files, a time-series results file, and
OpenFAST-related files. The following sections detail the purpose and
contents of these files.

Echo File
---------

If **Echo** = TRUE in the FAST.Farm primary input file, the contents of
the file will be echoed to a file with the naming convention
<*RootName.ech*, where <*RootName*> is as defined in
:numref:`FF:AmbWindVTK`. The echo file is helpful for
debugging the primary input file. The contents of an echo file will be
truncated if FAST.Farm encounters an error while parsing the primary
input file. The error usually corresponds to the line after the last
successfully echoed line.

.. _FF:Output:Sum:

Summary File
------------

If **SumPrint** = TRUE in the FAST.Farm primary input file, FAST.Farm
will generate a summary file with the naming convention of
<*RootName.sum*>. This file summarizes key information about the wind
farm model, including the wind turbine locations and OpenFAST primary
input files; wake dynamics finite-difference grid and parameters; time
steps of the various model components; and the name units and order of
the outputs that have been selected.

.. _FF:Output:Vis:

Visualization Output Files
--------------------------

If **WrDisWind** = TRUE in the FAST.Farm primary input file, FAST.Farm
will generate full 3D low- and high-resolution disturbed wind data
output files, i.e., the ambient wind and wake interactions across the
wind farm for visualization. The VTK data format and spatial resolutions
(number of grid points, origin, and spacing) of these output files
matches those of the corresponding low- and high-resolution ambient wind
data used by the FAST.Farm simulation. The VTK files are written to a
directory named *vtk_ff* where the FAST.Farm primary file is stored. The
naming conventions of these output files are
*<RootName>.Low.Dis.<n*\ :sub:`low`\ *>.vtk* and
*<RootName>.HighT<n*\ :sub:`t`\ *>.Dis.<n*\ :sub:`high`\ *>.vtk* for
the low- and high-resolution disturbed wind data files, respectively,
where *<n*\ :sub:`t`\ *>*, *<n*\ :sub:`low`\ *>*, and
*<n*\ :sub:`high` are as defined in
:numref:`FF:AmbWindVTK`, but with leading zeros.

- Likewise, if **NOutDisWindXY**, **NOutDisWindYZ**, or
  **NOutDisWindXZ** are set to be greater than zero in the FAST.Farm
  primary input file, FAST.Farm will generate low-resolution disturbed
  wind data (including wakes) output files that are 2D slices of the
  full low-resolution domain. The 2D slices are parallel to the *X-Y*,
  *Y-Z*, and/or *X-Z* planes of the global inertial-frame coordinate
  system, respectively. The VTK files are written to a directory named
  *vtk_ff* where the FAST.Farm primary file is stored. The naming
  conventions of these output files are
  *<RootName>.Low.DisXY<n*\ :sub:`Out`\ *>.<n*\ :sub:`low`\ *>.vtk*,
  *<RootName>.Low.DisYZ<n*\ :sub:`Out`\ *>.<n*\ :sub:`low`\ *>.vtk*,
  and
- *<RootName>.Low.DisXZ<n*\ :sub:`Out`\ *>.<n*\ :sub:`low`\ *>.vtk*
  for the *X-Y*, *Y-Z*, and *X-Z* slices, respectively, where
  *<n*\ :sub:`Out`\ *>* is as defined in
  :numref:`FF:AmbWindVTK`, but with leading zeros.

The time step (inverse of the frame rate) of all disturbed wind data
files is set by input parameter **WrDisDT** in the FAST.Farm primary
input file. Note that the full high-resolution disturbed wind data
output files are not output at a frame rate of :math:`1/`\ **DT_High**,
but are only output every **WrDisDT** seconds.

Each visualization output file follows the same VTK format used for the
ambient wind data files for the high-fidelity precursor simulations. See
:numref:`FF:AmbWindIfW` for details on the file format.

Visualizing the ambient wind and wake interactions can be useful for
interpreting results and debugging problems. However, FAST.Farm will
generate many files per output option when **WrDisWind** = TRUE and/or
when **NOutDisWindXY**, **NOutDisWindYZ**, and/or **NOutDisWindXZ** are
set greater than zero. This file generation will slow down FAST.Farm and
take up a lot of disk space, especially when generating full low- and
high-resolution disturbed wind data files. Therefore, disabling
visualization is recommended when running many FAST.Farm simulations.

.. _FF:Output:Time:

Time-Series Results File
------------------------

The FAST.Farm time-series results are written to an ASCII text-based
file with the naming convention <*RootName.out*. The results are in
table format, where each column is a data channel with the first column
always being the simulation time; each row corresponds to a simulation
output time step. The data channels are specified in the **OutList**
section of the *OUTPUT* section of the FAST.Farm primary input file. The
column format of the FAST.Farm-generated file is specified using the
**OutFmt** parameter of the FAST.Farm primary input file.

OpenFAST Output Files
---------------------

In addition to the FAST.Farm-generated output files, the OpenFAST model
of each wind turbine may also generate output files. The various output
files that OpenFAST may generate (at both the driver/glue code and
module levels, as specified within the OpenFAST input files) are
described in the OpenFAST documentation and include summary (*.sum*)
files, time-series results (ASCI *.out* or binary *.outb*) files,
visualization (*.vtk*) files, etc. FAST.Farm simulations will generate
these same files, but with the the path/rootname changed to *<RootName
of WT_FASTInFile>.T<n*\ :sub:`t`\ *>*.
