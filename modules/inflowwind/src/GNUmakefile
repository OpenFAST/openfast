AMREX_HOME ?= ../amrex

DEBUG	= TRUE

DIM	= 3

COMP    = gcc

USE_MPI   = TRUE
USE_OMP   = FALSE

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

F90EXE_sources += main.F90

CEXE_sources += read_amrex_subdomain.cpp
F90EXE_sources += read_amrex_subdomain_module.F90

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/F_Interfaces/Base/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
