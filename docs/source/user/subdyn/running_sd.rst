.. _running-subdyn:

Running SubDyn
===============

This section discusses how to obtain and execute SubDyn from a personal
computer. Both the stand-alone version and the FAST-coupled version of
the software are considered.

Downloading the SubDyn Software
--------------------------------

There are two forms of the SubDyn software to choose from: stand alone
and coupled to the FAST simulator. Although the user may not necessarily
need both forms, he/she would likely need to be familiar with and run
the stand-alone model if building a model of the substructure from
scratch. The stand-alone version is also helpful for model
troubleshooting and may benefit users who are interested in conducting
aero-hydro-servo-elastic simulations of an offshore wind turbine. 

Users can refer to the OpenFAST installation to download and compile SubDyn. 


Running SubDyn
---------------

Running the Stand-Alone SubDyn Program
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The stand-alone SubDyn program, *SubDyn\_win32.exe*, simulates
substructure dynamic responses of the user’s input model, without
coupling to FAST. Unlike the coupled version, the stand-alone software
requires the use of a driver file in addition to the primary SubDyn
input file. This driver file specifies inputs normally provided to
SubDyn by FAST, including motions of the TP reference point. Both the
SubDyn summary file and the results output file are available when using
the stand-alone SubDyn (see Section 4 for more information regarding the
SubDyn output files).

Run the standalone SubDyn software from a DOS command prompt by typing,
for example:

.. code-block:: bash
	
    >SubDyn_win32.exe MyDriverFile.dvr

where, *MyDriverFile.dvr* is the name of the SubDyn driver file, as
described in :numref:`sd_main-input-file`. The SubDyn primary input file is described in
Section :numref:`sd_driver-input-file`.


Running SubDyn Coupled to FAST  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run the coupled FAST software from a DOS command prompt by typing, for
example:

.. code-block:: bash

    >FAST_Win32.exe Test21.fst

where, *Test21.fst* is the name of the primary FAST input file. This
input file has a feature switch to enable or disable the SubDyn
capabilities within FAST, and a corresponding reference to the SubDyn
input file. See the documentation supplied with FAST for further
information.
