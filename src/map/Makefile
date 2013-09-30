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

# Python header files
PYTHON_INC = -I/usr/include/python2.7
#PYTHON_INC = 

# Boost Python include files and library
#BOOST_INC = /home/marco/Boost/boost_1_50_0/
#BOOST_LIB = /home/marco/Boost/boost_1_50_0/stage/lib
# Use this for Ubuntu 12.04 LTS
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib

# compiler options
COMPILER  = g++
OPT       = -g -fPIC -O2 -Wall -D WITH_PYTHON -D DEBUG -Wno-uninitialized -Wno-unused-variable
#OPT       = -g -fPIC -O2 -Wall -D DEBUG -Wno-uninitialized -Wno-unused-variable
LDFLAGS   = -fPIC
LIBS      = -lm 

include $(PETSC_DIR)/conf/variables
include $(PETSC_DIR)/conf/rules

OBJ = MAP_Message.o \
	MAP_InitInputType.o \
	VarType.o \
	Node.o \
	Element.o \
	FunctionEvals.o \
	MAP_BaseType.o \
	UserData.o \
	MAP_OtherStateType.o \
	Numerics.o \
	NWTCFunctions.o \
	PythonBinding.o \
	MAP_Types.o \
	MAP_FortranBinding.o \
	MAP_PythonBinding.o 

ALL: all

#all: $(OBJ)
#	-$(CLINKER) -shared $(OPT) -Wl,-no-undefined $(LDFLAGS) -lpython2.7 $(LIBS) $(OBJ) -o MAP.so $(PETSC_SNES_LIB) libboost_python-gcc45-sd-1_50.a liblapack.a
# Use this for Ubuntu 12.04 LTS
all: $(OBJ)
	-$(CLINKER) -shared -Wl,--export-dynamic $(OBJ) -L$(BOOST_LIB) -lboost_python -L/usr/lib/python2.7/config -lpython2.7 $(LIBS) -o MAP.so $(PETSC_SNES_LIB)

#  ============================================================================
#      Object Files     
#  ============================================================================

Element.o: Element.cpp
	$(COMPILER) $(OPT) -c Element.cpp

FunctionEvals.o: FunctionEvals.cpp
	$(COMPILER) $(OPT) -c FunctionEvals.cpp

MAP_InitInputType.o: MAP_InitInputType.cpp
	$(COMPILER) $(OPT) -c MAP_InitInputType.cpp

MAP_Message.o: MAP_Message.cpp
	$(COMPILER) $(OPT) -c MAP_Message.cpp 

Node.o: Node.cpp
	$(COMPILER) $(OPT) -c Node.cpp

VarType.o: VarType.cpp
	$(COMPILER) $(OPT) -c VarType.cpp

MAP_BaseType.o: MAP_BaseType.cpp
	$(COMPILER) $(OPT) $(PYTHON_INC) -c MAP_BaseType.cpp

MAP_OtherStateType.o: MAP_OtherStateType.cpp
	$(COMPILER) $(OPT) $(PYTHON_INC) $(PETSC_CC_INCLUDES) -c MAP_OtherStateType.cpp

Numerics.o: Numerics.cpp
	$(COMPILER) $(OPT) $(PYTHON_INC) $(PETSC_CC_INCLUDES) -c Numerics.cpp 

NWTCFunctions.o: NWTCFunctions.cpp
	$(COMPILER) $(OPT) $(PYTHON_INC) $(PETSC_CC_INCLUDES) -c NWTCFunctions.cpp

PythonBinding.o: PythonBinding.cpp
	$(COMPILER) $(OPT) $(PYTHON_INC) -I$(BOOST_INC) $(PETSC_CC_INCLUDES) -c PythonBinding.cpp

UserData.o: UserData.cpp
	$(COMPILER) $(OPT) $(PYTHON_INC) -I$(BOOST_INC) $(PETSC_CC_INCLUDES) -c UserData.cpp

MAP_Types.o: MAP_Types.c
	-gcc        $(OPT) $(PYTHON_INC) -I$(BOOST_INC) $(PETSC_CC_INCLUDES) -c MAP_Types.c

MAP_FortranBinding.o: MAP_FortranBinding.cpp
	$(COMPILER) $(OPT) $(PYTHON_INC) -I$(BOOST_INC) $(PETSC_CC_INCLUDES) -c MAP_FortranBinding.cpp

MAP_PythonBinding.o: MAP_PythonBinding.cpp
	$(COMPILER) $(OPT) $(PYTHON_INC) -I$(BOOST_INC) $(PETSC_CC_INCLUDES) -c MAP_PythonBinding.cpp


#  ============================================================================
#      Clean
#  ============================================================================
clean::		
	rm -f *.o *.so *.*~


#  ============================================================================
#      Memory Leak Checker  
#  ============================================================================
memcheck::
	valgrind --leak-check=yes python ./MAP_start.py bexco_low