.. _WaveTank:

WaveTank
========

The WaveTank glue-code is an experimental code for coupling hardware-in-the-loop
MHK models in a wavetank to software simulating the MHK turbine loads that
cannot be physically modeled in the wave tank.  The  *OpenFAST* modules
*SeaState*, *AeroDyn*, *MoorDyn*, and *InflowWind* are statically linked into a
single dynamic library (``cmake`` target ``wavetanktesting_c_binding``) with a
c-binding based interface.  This library can be called from *LabView* or another
code.

Inputs to the library include the time and motions, including the velocities and
accelerations, located at a single reference poitn at each time step.  The
resulting forces and moments are returned to the calling code.

Restrictions
~~~~~~~~~~~~
The current setup WaveTank library has several restrictions:

- rigid structure including platform, tower, and nacelle
- no yaw DOF
- rigid rotor
- constant rotor RPM for entire simulation
- no option for controller interfacing at present
- visualization limitted to *AeroDyn* and *SeaState*
- Current implementation only supports floating MHK turbines (``MHK = 2``). Other modes are present but not fully implemented.




Input File
~~~~~~~~~~


.. toctree::
   :maxdepth: 2
   
   wavetank_input.rst
