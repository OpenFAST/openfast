# NWTC Library Module
The legacy version of this module and additional documentation are available
at the [NWTC Software Portal](https://nwtc.nrel.gov/NWTC_Library/).

## Overview
Over the years, the researchers at the NWTC have written many general-purpose
routines that we use in many of our software. In the past, we would copy those
needed for any given program from the source of another code. We decided that
doing it that way was just too much work; especially when we wanted to change
one of the routines and have it affect all the packages. So, we decided to
create a central set of routines and include their source files in our Fortran
projects. Thus, the NWTC Subroutine Library was born.

The NWTC Subroutine Library consists of several modules, each contained in its
own source file. The library is written in standard Fortran 2003, with the
exception of some of the routines and data in the Sys*.f90 files. We currently
have several compiler-specific files available to choose from, including
- SysIVF.f90 for the Intel Visual Fortran compiler for Windows
- SysGnuLinux.f90 for the GNU Fortran compiler for Linux

If you want to port our programs to a different compiler, you should only have
to modify one of the Sys*.f90 files to get it to compile.
