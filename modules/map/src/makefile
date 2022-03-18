#   Copyright (C) 2014 mdm                                     
#   map[dot]plus[dot]plus[dot]help[at]gmail                                
#                                                              
# Licensed to the Apache Software Foundation (ASF) under one   
# or more contributor license agreements.  See the NOTICE file 
# distributed with this work for additional information        
# regarding copyright ownership.  The ASF licenses this file   
# to you under the Apache License, Version 2.0 (the            
# "License"); you may not use this file except in compliance   
# with the License.  You may obtain a copy of the License at   
#                                                              
#   http://www.apache.org/licenses/LICENSE-2.0                 
#                                                              
# Unless required by applicable law or agreed to in writing,   
# software distributed under the License is distributed on an  
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       
# KIND, either express or implied.  See the License for the    
# specific language governing permissions and limitations            
# under the License.                                             

MAP_VER:=1.20.10

BITS = -m64
PLATFORM = $(shell uname -s)
VPATH = cminpack:simclist:bstring:lapack


ifeq ($(OS),Windows_NT)
  DEL_CMD   = del
  LIB_EXT   = dll
  CFLAGS    = $(BITS) -g -std=c99 -DMAP_DLL_EXPORTS -DCMINPACK_NO_DLL  -DNDEBUG -D_WINDOWS -D_USRDLL -D_MINGW
  LDFLAGS   = $(BITS) -g -shared -Wl,--export-all-symbols
  LIB_FLAGS :=
  CC_TOOLS  = gcc
else
  PLATFORM = $(shell uname -s)
  DEL_CMD   = rm -rf
  GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)

  ifeq ($(PLATFORM),Darwin)
    LIB_EXT   = dylib
    CC_TOOLS  = clang
    CFLAGS    = -g -O1 -fno-omit-frame-pointer -fPIC -D DEBUG -Icminpack -Isimclist
    LDFLAGS   = -g -dynamiclib
    LIB_FLAGS := -lm
  else ifeq ($(PLATFORM),Linux)
    LIB_EXT   = so
    CC_TOOLS  = gcc
    #CFLAGS    = -g -O1 -fsanitize=address -fno-omit-frame-pointer -fPIC -D DEBUG -DGITVERSION=\"$(GIT_VERSION)\" -Icminpack -Isimclist
    #CFLAGS    = -g -fPIC -std=c99 -Wuninitialized -Wall -D DEBUG -DGITVERSION=\"$(GIT_VERSION)\"
    CFLAGS    = $(BITS) -g -fPIC -std=c99 -D DEBUG -DGITVERSION=\"$(GIT_VERSION)\" -D WITH_LAPACK
    LDFLAGS   = $(BITS) -g -shared -llapacke
    LIB_FLAGS := -lm -llapacke
  endif
endif

DEBUG     = 
OBJ       = lmder.o \
	dpmpar.o \
	lmpar.o \
	qrsolv.o \
	enorm.o \
	enorm_u.o \
	qrfac.o \
	simclist.o \
	bstrlib.o \
	bstraux.o \
	freedata.o \
	mapinit.o \
	maperror.o \
	lineroutines.o \
	numeric.o \
	outputstream.o \
	mapapi.o \
	lmroutines.o \
	jacobian.o \
	residual.o \


all : $(OBJ)
	$(CC_TOOLS) $(LDFLAGS) -o libmap-$(MAP_VER).$(LIB_EXT) $(DEBUG) $(OBJ) $(LIB_FLAGS) 

.c.o :
	$(CC_TOOLS) -c $(CFLAGS) $<

clean:
	$(DEL_CMD) *.$(LIB_EXT) *.o *~

memcheck:	
	valgrind --tool=memcheck --leak-check=yes python -E -tt ./main.py	
