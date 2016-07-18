# NWTC Subroutine Library
A library of general-use Fortran 2003 routines used in many NWTC computer-aided engineering tools

**Authors**: [Bonnie Jonkman](mailto:bonnie.jonkman@nrel.gov) and Marshall Buhl, NREL

The NWTC Subroutine Library consists of several modules, each contained in its own source file. The library is written in standard Fortran 2003, 
with the exception of some of the routines and data in the Sys\*.f90 files. We currently have several compiler-specific files available to choose 
from, including SysIVF.f90 for the Intel Visual Fortran compiler for Windows and SysGnuLinux.f90 for the GNU Fortran compiler for Linux. 
If you want to port our programs to a different compiler, you should only have to modify one of the Sys\*.f90 files to get it to compile.

For more information, please refer to documentation on the [NWTC Subroutine Library web site](https://nwtc.nrel.gov/NWTC_Library "NWTC_Library").
Source-code documentation can be found here: [Doxygen documentation] (http://wind.nrel.gov/nwtc/doxygen/nwtc_library/html/index.html).

