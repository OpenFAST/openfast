.. _ifw_appendix:

Appendix
========

.. _ifw_input_files:

InflowWind Input Files
----------------------

In this appendix we describe the InflowWind input-file structure and provide examples.

1) InflowWind Driver Input File 
:download:`(driver input file example) <examples/inflowwind_driver_example.inp>`: 

The driver input file is needed only for the standalone version of InflowWind. It contains inputs regarding the InflowWind file, interpolation parameters, and the desired output files. 
The InflowWind driver can also be run without this file by using command-line arguments instead.

2) InflowWind Primary Input File 
:download:`(primary input file example) <examples/inflowwind_example.dat>`: 

The primary InflowWind input file defines the inflow that is generated or read from other files. The InflowWind file contains sections for each type of wind-file format.

3) Native Bladed Scaling File
:download:`(primary input file example) <examples/inflowwind_bladedscaling_example.dat>`: 

This file includes lines that determine how to scale the non-dimensional full-field
turbulence files from Bladed. 

4) Uniform Wind Data File
:download:`(uniform wind input file example) <examples/inflowwind_uniform_example.dat>`: 

This file includes lines that define uniform (deterministic) wind data files.


.. _ifw_output_channels:

InflowWind List of Output Channels
----------------------------------

This is a list of all possible output parameters for the InflowWind module. 
See the InflowWind tab of the 
:download:`(OutListParameters.xlsx file) <../../../OtherSupporting/OutListParameters.xlsx>`: 

