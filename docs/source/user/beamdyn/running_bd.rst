.. _running-beamdyn:

Running BeamDyn
===============

This section discusses how to obtain and execute BeamDyn from a personal
computer. Both the stand-alone version and the FAST-coupled version of
the software are considered.

Downloading the BeamDyn Software
--------------------------------

There are two forms of the BeamDyn software to choose from: stand-alone
and coupled to the FAST simulator. Although the user may not necessarily
need both forms, he/she would likely need to be familiar with and run
the stand-alone model if building a model of the blade from scratch. The
stand-alone driver is also helpful for model troubleshooting, even if
the goal is to conduct aero-hydro-servo-elastic simulations of
onshore/offshore wind turbines within FAST.

Running BeamDyn
---------------

Running the Stand-Alone BeamDyn Program
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The stand-alone BeamDyn program, ``BeamDyn_Driver.exe``, simulates static
and dynamic responses of the user's input model, without coupling to
FAST. Unlike the coupled version, the stand-alone software requires the
use of a driver file in addition to the primary and blade BeamDyn input
files. This driver file specifies inputs normally provided to BeamDyn by
FAST, including motions of the blade root and externally applied loads.
Both the BeamDyn summary file and the results output file are available
when using the stand-alone BeamDyn (see Section :ref:`bd-output-files` for
more information regarding the BeamDyn output files).

Run the stand-alone BeamDyn software from a DOS command prompt by
typing, for example:

.. code-block:: bash

    >BeamDyn_Driver.exe Dvr_5MW_Dynamic.inp

where, ``Dvr_5MW_Dynamic.inp`` is the name of the BeamDyn driver input
file, as described in Section :ref:`driver-input-file`.

Running BeamDyn Coupled to FAST
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the coupled FAST software from a command prompt by typing, for
example:

.. code-block:: bash

    >openfast.exe Test26.fst

where ``Test26.fst`` is the name of the primary FAST input file. This
input file has a feature switch to enable or disable the BeamDyn
capabilities within FAST, and a corresponding reference to the BeamDyn
input file. See the documentation supplied with FAST for further
information.
