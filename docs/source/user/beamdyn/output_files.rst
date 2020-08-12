.. _bd-output-files:
   
Output Files
============

BeamDyn produces three types of output files, depending on the options
selected: an echo file, a summary file, and a time-series results file.
The following sections detail the purpose and contents of these files.

Echo File
---------

If the user sets the ``Echo`` flag to TRUE in the BeamDyn primary
input file, the contents of this file will be echoed to a file with the
naming convention ``InputFile.ech``. The echo file is helpful for
debugging the input files. The contents of an echo file will be
truncated if BeamDyn encounters an error while parsing an input file.
The error usually corresponds to the line after the last successfully
echoed line.

.. _sum-file:

Summary File
------------

In stand-alone mode, BeamDyn generates a summary file with the naming
convention, ``InputFile.sum`` if the ``SumPrint`` parameter is set
to TRUE. When coupled to FAST, the summary file is named
``InputFile.BD.sum``. This file summarizes key information about the
simulation, including:

-  Blade mass.

-  Blade length.

-  Blade center of mass.

-  Initial global position vector in BD coordinate system.

-  Initial global rotation tensor in BD coordinate system.

-  Analysis type.

-  Numerical damping coefficients.

-  Time step size.

-  Maximum number of iterations in the Newton-Raphson solution.

-  Convergence parameter in the stopping criterion.

-  Factorization frequency in the Newton-Raphson solution.

-  Numerical integration (quadrature) method.

-  FE mesh refinement factor used in trapezoidal quadrature.

-  Number of elements.

-  Number of FE nodes.

-  Initial position vectors of FE nodes in BD coordinate system.

-  Initial rotation vectors of FE nodes in BD coordinate system.

-  Quadrature point position vectors in BD coordinate system. For Gauss
   quadrature, the physical coordinates of Gauss points are listed. For
   trapezoidal quadrature, the physical coordinates of the quadrature
   points are listed.

-  Sectional stiffness and mass matrices at quadrature points in local
   blade reference coordinate system. These are the data being used in
   calculations at quadrature points and they can be different from the
   section in Blade Input File since BeamDyn linearly interpolates the
   sectional properties into quadrature points based on need.

-  Initial displacement vectors of FE nodes in BD coordinate system.

-  Initial rotational displacement vectors of FE nodes in BD coordinate
   system.

-  Initial translational velocity vectors of FE nodes in BD coordinate
   system.

-  Initial angular velocity vectors of FE nodes in BD coordinate system.

-  Requested output information.

All of these quantities are output in this file in the BD coordinate
system, the one being used internally in BeamDyn calculations. The
initial blade reference coordinate system, denoted by a subscript
:math:`r0` that follows the IEC standard, is related to the internal BD
coordinate system by :numref:`IECBD` in :numref:`beamdyn-theory`.

Results File
------------

The BeamDyn time-series results are written to a text-based file with
the naming convention ``DriverInputFile.out`` where
``DriverInputFile`` is the name of the driver input file when BeamDyn
is run in the stand-alone mode. If BeamDyn is coupled to FAST, then FAST
will generate a master results file that includes the BeamDyn results.
The results in ``DriverInputFile.out`` are in table format, where each
column is a data channel (the first column always being the simulation
time), and each row corresponds to a simulation time step. The data
channel are specified in the OUTPUT section of the primary input file.
The column format of the BeamDyn-generated file is specified using the
``OutFmt`` parameters of the primary input file.

