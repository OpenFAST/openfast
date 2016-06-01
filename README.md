# FAST Registry for Automatic Code Generation
A utility for generating code for module data types in the [FAST Modularization Framework](https://nwtc.nrel.gov/FAST-Developers "FAST Developers")

**Authors**: John Michalakes and [Bonnie Jonkman](mailto:bonnie.jonkman@nrel.gov), NREL

The FAST Registry allows the developer to specify the data types for a module once and in a single 
location, automating the time consuming and error-prone task of generating the code.  The tables 
in the Registry text input files serve as a data dictionary for improving understandability and maintainability of 
the code.

The FAST Registry is borrowed from a mechanism that was originally developed at NCAR for 
the [Weather Research and Forecast (WRF) model software](http://www.mmm.ucar.edu/wrf/WG2/software_2.0/registry_schaffer.pdf).
It has been modified for the FAST Modularization Framework to create *ModuleName*_Types.f90 files along with any necessary C source code or header files associated
with the *ModuleName*_Types.f90 files.


For more information and syntax, please refer to the [NWTC Programmer's Handbook](https://nwtc.nrel.gov/system/files/ProgrammingHandbook_Mod20130717.pdf).
