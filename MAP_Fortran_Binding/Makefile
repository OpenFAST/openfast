# =========== MAP Makefile        =============================================
#	     
#  Copyright Sept. 2012
#  National Renewable Energy Laboratory, Golden, Colorado, USA
#  Author: Marco D. Masciola
# 	
#  MAP is free software: you can redistribute it and/or modify it under the terms 
#  of the GNU General Public License as published by the Free Software Foundation, 
#  either version 3 of the License, or (at your option) any later version.
#
#  MAP is distributed in the hope that it will be useful, but WITHOUT ANY 
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
#  FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License along 
#  with MAP. If not, see:
#  
#  < http://www.gnu.org/licenses/>
#  ============================================================================


REGISTRY = /media/sf_Virtual_Machine/MAP/MSQS/FAST_Dev/Registry/Source/registry.exe
NWTC_LIB = /media/sf_Virtual_Machine/MAP/MSQS/FAST_Dev/NWTC_Lib_trunk/source

PREC        = SingPrec
SYS_FILE    = SysGnuLinux

FC       = gfortran
FCFLAGS  = -O0 -fbacktrace -fbounds-check -x f95-cpp-input -ffree-line-length-none 
OPT      = -fmax-errors=20

OBJ = $(PREC).o         \
	$(SYS_FILE).o   \
	NWTC_IO.o       \
	ModMesh_Types.o \
	ModMesh.o       \
	NWTC_Num.o      \
	NWTC_Library.o  \
	MAP_C_Types.o   \
	MAP_Types.o     \
	MAP.o 	        \
	MAP_Glue.o


all : $(OBJ)
	$(FC) $(OPT) -o prog $(OBJ) -L. /media/sf_Virtual_Machine/MAP/MSQS/MSQS_Program/MAP_PETSc/MAP.so


MAP_Glue.o : MAP_Glue.f90
	$(FC) $(OPT) $(FCFLAGS) -c MAP_Glue.f90

$(PREC).o : $(NWTC_LIB)/$(PREC).f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/$(PREC).f90 

$(SYS_FILE).o : $(NWTC_LIB)/$(SYS_FILE).f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/$(SYS_FILE).f90 

NWTC_IO.o : $(NWTC_LIB)/NWTC_IO.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/NWTC_IO.f90 

ModMesh_Types.o : $(NWTC_LIB)/ModMesh_Types.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/ModMesh_Types.f90 

ModMesh.o : $(NWTC_LIB)/ModMesh.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/ModMesh.f90 

NWTC_Num.o : $(NWTC_LIB)/NWTC_Num.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/NWTC_Num.f90 

NWTC_Library.o : $(NWTC_LIB)/NWTC_Library.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/NWTC_Library.f90 

MAP_C_Types.o : MAP_C_Types.f90
	gfortran $(OPT) $(FCFLAGS) -c MAP_C_Types.f90

MAP_Types.o : MAP_Types.f90
	$(FC) $(OPT) $(FCFLAGS) -c MAP_Types.f90

MAP.o : MAP.f90
	$(FC) $(OPT) $(FCFLAGS) -c MAP.f90





clean:		
	rm -f *.o prog *.*~ *.*# Makefile~ *.mod *.valgrind

memcheck:
	valgrind --log-file=valgrind.out --leak-check=yes --show-reachable=yes --track-origins=yes prog

map_registry:
	$(REGISTRY) -ccode MAP_Registry.txt




#all : MAP_C_Types.o MAP_Types.o MAP_Glue.o .o 
#	gfortran $(OPT) -o prog MAP_C_Types.o MAP_Types.o MAP.o Driver.o -L. /media/sf_Virtual_Machine/MAP/MSQS/MSQS_Program/MAP_PETSc/MAP.so
#
#MAP_Glue.o : MAP_Glue.f90
#	gfortran $(OPT) $(FCFLAGS) -c MAP_Glue.f90
#
#MAP_Types.o : MAP_Types.f90
#	gfortran $(OPT) $(FCFLAGS) -c MAP_Types.f90
#
#MAP_C_Types.o : MAP_C_Types.f90
#	gfortran $(OPT) $(FCFLAGS) -c MAP_C_Types.f90
#
#MAP.o : MAP.f90
#	gfortran $(OPT) $(FCFLAGS) -c MAP.f90