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
stand-alone version is also helpful for model troubleshooting, even if
the goal is to conduct aero-hydro-servo-elastic simulations of
onshore/offshore wind turbines within FAST.

Stand-Alone BeamDyn Archive
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users can download the stand-alone BeamDyn archive from our Web server
at https://nwtc.nrel.gov/BeamDyn. The file has a name similar to
``BD_v1.00.00a.exe``, but may have a different version number. The user
can then download the self-extracting archive (*.exe*) to expand the
archive into a folder he/she specifies.

The archive contains the ``bin``, ``CertTest``, ``Compiling``,
``Docs``, and ``Source`` folders. The ``bin`` folder includes the
main executable file, ``BeamDyn_Driver.exe``, which is used to execute
the stand-alone BeamDyn program. The ``CertTest`` folder contains a
collection of sample BeamDyn input files and driver input files that can
be used as templates for the user’s own models. This document may be
found in the ``Docs`` folder. The ``Compiling`` folder contains files
for compiling the stand-alone ``BeamDyn_v1.00.00.exe`` file with either
Visual Studio or gFortran. The Fortran source code is located in the
``Source`` folder.

FAST Archive
~~~~~~~~~~~~

Download the FAST archive, which includes BeamDyn, from our Web server
at https://nwtc.nrel.gov/FAST8. The file has a name similar to
``FAST_v8.12.00.exe``, but may have a different version number. Run the
downloaded self-extracting archive (``.exe``) to expand the archive into a
user-specified folder. The FAST executable file is located in the
archive’s ``bin`` folder. An example model using the NREL 5-MW
reference turbine is located in the ``CertTest`` folder.

Running BeamDyn
---------------

Running the Stand-Alone BeamDyn Program
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The stand-alone BeamDyn program, ``BeamDyn_Driver.exe``, simulates static
and dynamic responses of the user’s input model, without coupling to
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

Run the coupled FAST software from a DOS command prompt by typing, for
example:

.. code-block:: bash

    >FAST_Win32.exe Test26.fst

where ``Test26.fst`` is the name of the primary FAST input file. This
input file has a feature switch to enable or disable the BeamDyn
capabilities within FAST, and a corresponding reference to the BeamDyn
input file. See the documentation supplied with FAST for further
information.
