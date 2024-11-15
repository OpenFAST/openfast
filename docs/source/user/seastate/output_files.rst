.. _sea-output:

Output Files
============
SeaState produces five types of output files: echo files, a wave-elevation 
file, wave-kinematics files, a summary file, and a time-series file. 
The following sections detail the purpose and contents of these files.

Echo Files
~~~~~~~~~~
If you set the **Echo** flag to TRUE in the SeaState driver file or the
SeaState primary input file, the contents of those files will be echoed
to a file with the naming conventions, **OutRootName**\ *.dvr.ech* for
the driver input file and **OutRootName**\ *.SeaSt.ech* for the SeaState
primary input file. **OutRootName** is either specified in the SEASTATE
section of the driver input file, or by the OpenFAST program. The echo files
are helpful for debugging your input files. The contents of an echo file
will be truncated if SeaState encounters an error while parsing an input
file. The error usually corresponds to the line after the last
successfully echoed line.

Wave-Elevation File
~~~~~~~~~~~~~~~~~~~~
Setting **WaveElevSeriesFlag** in the driver file to TRUE enables the
output of the wave elevation across the entire SeaState wave grid to 
a text-based file with the name ``OutRootName.WaveElev.out``. 
The output wave grid is controlled by the SPATIAL DISCRETIZATION 
section of the primary SeaState input file. These wave elevations are 
distinct and output separately from the wave elevations determined by 
**NWaveElev** in the SeaState primary input file, which are instead 
recorded in the time-series output file. The wave-elevation output file 
``OutRootName.WaveElev.out`` contains the total wave elevation, which is 
the sum of the first- and second-order terms (when the second-order wave 
kinematics are optionally enabled). The wave-elevation file described 
here is useful for visualization purposes.

Wave-Kinematics Files
~~~~~~~~~~~~~~~~~~~~~
**WrWvKinMod** controls the wave kinematics output from the SeaState driver. 
The driver will output the wave-elevation time series at the global origin 
(0,0) in a separate *.Elev* file if **WrWvKinMod** = 1. This file also serves 
as a valid **WvKinFile** for **WaveMod** = 5 (externally generated wave-elevation 
time series) in the primary SeaState input file. This is a separate wave-
elevation output independent from the wave-field output obtained with 
**WaveElevSeriesFlag** = TRUE and the wave-elevation time series obtained with 
**NWaveElev** in the SeaState primary input file. 

If **WrWvKinMod** = 2, SeaState will output the full wave kinematics (velocity, 
acceleration, dynamic pressure, and wave elevation) at all wave grid points in 
eight output files with the extensions *.Vxi*, *.Vyi*, *.Vzi*, *.Axi*, *.Ayi*, 
*.Azi*, *.DynP*, and *.Elev*. The velocity and acceleration outputs are all in the 
global earth-fixed coordinate system. These files are also valid as 
**WvKinFile** for **WaveMod** = 6 (externally generated full wave-kinematics 
time series) and can be used as templates if users would like to build 
their own input files for **WaveMod** = 6.

.. _sea-summary-file:

Summary File
~~~~~~~~~~~~
SeaState generates a summary file with the naming convention,
**OutRootName**\ *.SeaSt.sum* if the **SeaStSum** parameter is set 
to TRUE. This file summarizes key information about your sea-state 
model, including water density, water depth, still-water level, 
the wave grid, first-order wave frequency components, and any 
user-requested time-series outputs.

Summary of User-Requested Outputs
---------------------------------
The summary file includes information about all user-requested 
wave-elevation and wave-kinematics time-series output channels.

Wave-Kinematics Output Locations
++++++++++++++++++++++++++++++++
The first column lists the index of the wave-kinematics output locations, 
as entered in the OUTPUT CHANNELS section of the SeaState input file. 
**Xi**, **Yi**, and **Zi** provide the spatial coordinates of the wave-kinematics 
output locations in the global inertial-frame coordinate system. 

Wave-Elevation Output Locations
+++++++++++++++++++++++++++++++
The first column lists the index of the wave-elevation output locations, 
as entered in the OUTPUT CHANNELS section of the SeaState input file. 
**Xi** and **Yi** provide the spatial coordinates of the wave-elevation 
output locations in the global inertial-frame coordinate system. 

Requested Output Channels
+++++++++++++++++++++++++
The string labels of all requested time-series output channels, as entered 
in the OUTPUT CHANNELS section of the SeaState input file, are summarized here.

Wave Frequency Components
-------------------------
This section provides the frequency-domain description (in terms of a
Discrete Fourier Transform or DFT) of the first-order wave elevation at
(0,0) on the free surface, but is not written when **WaveMod** = 0 or 6.
The first column, **index**, identifies the index of each wave frequency
component. The finite-depth wave number, angular frequency, and direction of the
wave components are given by **k**, **Omega**, and **Direction**,
respectively. The last two columns provide the real
(**REAL(DFT{WaveElev})**) and imaginary (**IMAG(DFT{WaveElev})**)
parts of the first-order DFT wave amplitudes. The DFT
produces both negative- and positive-frequency components.
The negative-frequency components are complex conjugates of the 
positive-frequency components because the time-domain wave elevation is
real-valued. The relationships between the negative- and
positive-frequency components of the DFT are given by
:math:`k\left( - \omega \right) = - k\left( \omega \right)` and
:math:`H\left( - \omega \right) = {H\left( \omega \right)}^{*}`, where
*H* is the DFT of the wave elevation and *\** denotes the complex
conjugate.


Results File
~~~~~~~~~~~~

The SeaState time-series results are written to a text-based file with
the naming convention ``OutRootName.SeaSt.out`` when **OutSwtch** is
set to either 1 or 3. If SeaState is coupled to OpenFAST and **OutSwtch** is
set to 2 or 3, then OpenFAST will generate a master results file that
includes the SeaState results. The results are in table format, where
each column is a data channel (the first column is always the
simulation time), and each row corresponds to a simulation output time
step. The data channels are specified in the OUTPUT CHANNELS section of
the SeaState primary input file. The column format of the
SeaState-generated file is specified using the **OutFmt** and
**OutSFmt** parameter of the primary input file.
