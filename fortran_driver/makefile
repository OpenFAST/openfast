# =========== MAP Makefile        =============================================
#	     
#  Copyright 2014
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

NWTC_LIB = ./nwtc_lib/src
NWTC_NUM = ./nwtc_lib/src/NetLib

PREC        = SingPrec
SYS_FILE    = SysGnuLinux

FC       = gfortran
FCFLAGS  = -g -O0 -fbacktrace -fbounds-check -x f95-cpp-input -ffree-line-length-none
OPT      = -fmax-errors=20

OBJ = $(PREC).o              \
	NWTC_Base.o          \
	$(SYS_FILE).o        \
	NWTC_Library_Types.o \
	NWTC_IO.o            \
	NWTC_Num.o           \
	NWTC_LAPACK.o        \
	ModMesh_Types.o      \
	ModMesh.o            \
	ModMesh_Mapping.o    \
	NWTC_Library.o       \
	MAP_Types.o          \
	map.o 	             \
	map_glue.o


all : $(OBJ)
	$(FC) $(OPT) -o f_driver $(OBJ) -L. ../src/libmap-1.20.00.so -lblas -llapack

map_glue.o : map_glue.f90
	$(FC) $(OPT) $(FCFLAGS) -c map_glue.f90

$(PREC).o : $(NWTC_LIB)/$(PREC).f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/$(PREC).f90 

$(SYS_FILE).o : $(NWTC_LIB)/$(SYS_FILE).f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/$(SYS_FILE).f90 

NWTC_Base.o : $(NWTC_LIB)/NWTC_Base.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/NWTC_Base.f90 

NWTC_LAPACK.o : $(NWTC_NUM)/NWTC_LAPACK.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_NUM)/NWTC_LAPACK.f90 

NWTC_IO.o : $(NWTC_LIB)/NWTC_IO.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/NWTC_IO.f90 

ModMesh_Types.o : $(NWTC_LIB)/ModMesh_Types.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/ModMesh_Types.f90 

ModMesh.o : $(NWTC_LIB)/ModMesh.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/ModMesh.f90 

ModMesh_Mapping.o : $(NWTC_LIB)/ModMesh_Mapping.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/ModMesh_Mapping.f90 

NWTC_Num.o : $(NWTC_LIB)/NWTC_Num.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/NWTC_Num.f90 

NWTC_Library.o : $(NWTC_LIB)/NWTC_Library.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/NWTC_Library.f90 

NWTC_Library_Types.o : $(NWTC_LIB)/NWTC_Library_Types.f90
	$(FC) $(OPT) $(FCFLAGS) -c $(NWTC_LIB)/NWTC_Library_Types.f90 

MAP_C_Types.o : ../registry/MAP_C_Types.f90
	gfortran $(OPT) $(FCFLAGS) -c ../registry/MAP_C_Types.f90

MAP_Types.o : ../registry/MAP_Types.f90
	$(FC) $(OPT) $(FCFLAGS) -c ../registry/MAP_Types.f90

map.o : map.f90
	$(FC) $(OPT) $(FCFLAGS) -c map.f90

clean:		
	rm -f *.o map_driver *.*~ *.*# Makefile~ *.mod *.valgrind baseline.sum*

memcheck:
	valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --track-origins=yes ./f_driver
