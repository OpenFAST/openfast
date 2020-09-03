.. _ifw_input:

InflowWind Input Files
======================

.. _ifw_native_bladed:

Native Bladed wind file support in InflowWind
---------------------------------------------

The ability to read native Bladed wind files (without scaling) has been added to InflowWind. 
To use this feature, the ``WindType`` must be set to ``7`` on line 5 of the primary
InflowWind input file. An example of this file is given inAn example of this Native Bladed scaling file is included in 
:numref:`ifw_appendix`.

::

   7                  WindType       - switch for wind file type (1=steady; 2=uniform; 3=binary TurbSim FF; 4=binary Bladed-style FF; 5=HAWC format; 6=User defined; 7=Bladed native)

In the section for ``WindType = 4``, the name of an intermediate Bladed
wind file should be given (including the file extension). The tower file
flag is ignored.

::

    ================== Parameters for Binary Bladed-style Full-Field files   [used only for WindType = 4] =========
    "tw06_80hh_s200.BladedWind.ipt"  FilenameRoot   - Name of the Full field wind file to use (.wnd, .sum)
     F                               TowerFile      - Have tower file (.twr) (flag)


The intermediate Bladed wind scaling file must contain the following information, which can be retrieved
directly from the Bladed project simulation file from the
``MSTART WINDSEL`` and ``MSTART WINDV`` sections. Additionally, the file
may include an ``XOFFSET`` line, which allows the wind to be shifted by
a given distance. If not included, ``XOFFSET`` is assumed to be 0.
An example of this Native Bladed scaling file is included in 
:numref:`ifw_appendix`.


::

    UBAR  12
    REFHT  90
    TI  0.033333
    TI_V  0.026667
    TI_W  0.016667
    WDIR  0
    FLINC  .139626222222222
    WINDF    "../tw06_80hh_s200.wnd"
    WSHEAR    .2
    XOFFSET  0

In the above file, the names correspond to the following:

+--------+-----------------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------+
| Line   | Variable Name   | Units   | Description                                                                                                                                     |
+========+=================+=========+=================================================================================================================================================+
| 1      | ``UBAR``        | (m/s)   | Mean wind speed                                                                                                                                 |
+--------+-----------------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------+
| 2      | ``REFHT``       | (m)     | Reference height (turbine hub height)                                                                                                           |
+--------+-----------------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------+
| 3      | ``TI``          | (-)     | Turbulence intensity in longitudinal (mean wind flow) direction                                                                                 |
+--------+-----------------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------+
| 4      | ``TI_V``        | (-)     | Turbulence intensity in horizontal direction (orthogonal to mean flow direction)                                                                |
+--------+-----------------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------+
| 5      | ``TI_W``        | (-)     | Turbulence intensity in vertical direction (orthogonal to mean flow direction)                                                                  |
+--------+-----------------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------+
| 6      | ``WDIR``        | (rad)   | Wind direction (meteorological rotation direction)                                                                                              |
+--------+-----------------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------+
| 7      | ``FLINC``       | (rad)   | Upflow angle (positive is up)                                                                                                                   |
+--------+-----------------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------+
| 8      | ``WINDF``       | (-)     | Name of native Bladed wind file (absolute or relative path, 200 character limit)                                                                |
+--------+-----------------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------+
| 9      | ``WSHEAR``      | (-)     | Power law wind shear exponent                                                                                                                   |
+--------+-----------------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------+
| 10     | ``XOFFSET``     | (m)     | Turbulence box offset in the X direction (how far ahead of the turbine the turbulence box starts). If missing, this value is assumed to be 0.   |
+--------+-----------------+---------+-------------------------------------------------------------------------------------------------------------------------------------------------+

Limitations: - Wind file is centered on hub height ("Best fit for rotor
and tower" not implemented) - Always allow wind file to wrap around
(unchecked box not implemented) - Only power-law wind profile is
implemented (not logarithmic, none, or user-defined)

