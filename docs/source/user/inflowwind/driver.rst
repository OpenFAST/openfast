InflowWind Driver
=================
Example input files are
included in :numref:`ifw_appendix`.

Command-line syntax for InflowWind driver:

::

    InlowWind_Driver <filename> [options]

          where:  <filename>     --   Name of driver input file to use
        options:  /ifw           --   treat <filename> as name of InflowWind input file (no driver input file)

        The following options will override values in the driver input file:
                  /DT[#]         --   timestep
                  /TStart[#]     --   start time
                  /TSteps[#]     --   number of timesteps
                  /xrange[#:#]   --   range of x (#'s are reals)
                  /yrange[#:#]   --   range of y
                  /zrange[#:#]   --   range in z (ground = 0.0)
                  /Dx[#]         --   spacing in x
                  /Dy[#]         --   spacing in y
                  /Dz[#]         --   spacing in z
                  /points[FILE]  --   calculates at x,y,z coordinates specified in a white space delimited FILE
                  /v             --   verbose output
                  /vv            --   very verbose output
                  /hawc          --   convert wind file specified in InflowWind to HAWC format
                  /bladed        --   convert wind file specified in InflowWind to Bladed format
                  /vtk           --   convert wind file specified in InflowWind to VTK format
                  /help          --   print this help menu and exit

::

    Notes:
    -   Unspecified ranges and resolutions default to what is in the file.
    -   If no XRange is specified, assumed to be only at X=0
    -   Options are not case sensitive.

The `InflowWind Manual <https://wind.nrel.gov/nwtc/docs/InflowWind_Manual.pdf>`__
contains a description of file formats that it can read.

Specifying the InflowWind Input File
------------------------------------

The InflowWind driver input file requires that an InflowWind input file
be specified within it. See an example InflowWind input
file in :numref:`ifw_appendix`.

Within the InflowWind input file, if the wind file being specified is
Bladed native format (``WindType = 7``), please also see 
:numref:`ifw_native_bladed`.

Wind-file output formats
------------------------

The InflowWind driver is capable of writing the wind data read from the
input wind file into wind files of various formats.

HAWC2
~~~~~

This format generates the following files:

- three binary files, one for each component:
  ``<RootName>-HAWC.u``, ``<RootName>-HAWC.v``, and ``<RootName>-HAWC.w``

- a text summary file in the style of HAWC2 input files:
  ``<RootName>-HAWC.sum``

In the conversion script, the u component will have the (approximate)
mean removed at each height. The mean value that was removed is
displayed as comments in the text summary file. The turbulence is not
scaled, so it will have the same scaling as the original file.

Bladed
~~~~~~

This format generates a packed binary file and a text summary file.

This output format is in the Bladed-style format that TurbSim generates. That
means that **the shear is included** in the file.

VTK
~~~

This format creates files in a subdirectory called ``vtk``. There is one
vtk file for each time in the full-field data structure, and the entire
Y-Z grid is printed in each file. This format can be used to visualize
the wind field using a viewer such as ParaView.

Converting uniform wind to full-field wind format
-------------------------------------------------

When converting from a uniform wind file to a full-field wind format,
the following assumptions are used: - The advection speed is the
time-averaged horizontal wind speed in the uniform wind file (it does
not include the gust speed). - The constant time-step used in the output
file is the smallest difference between any two entries in the
hub-height file. - The maximum time in the uniform wind file will be
used as the maximum time in the FF binary file. - The grid is generated
with 5-m resolution in the lateral (Y) and horizontal (Z) directions. -
The size of the grid is based on the ``RefLength`` parameter in the
InflowWind input file. The converter adds approximately 10% to the grid
width, with the exact size determined by achieving the desired grid
resolution. The grid is centered in the lateral direction; it extends
vertically above ``RefHt`` by the same distance as the grid width, and
extends below ``RefHt`` to the ground (or within one grid point of the
ground).

Note that there is a potential time shift between the uniform and
full-field wind files, equal to the time it takes to travel the distance
of half the grid width. When using the resulting full-field files, care
must be taken that the aeroelastic code does not treat it as periodic.

