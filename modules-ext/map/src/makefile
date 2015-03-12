#   Copyright (C) 2014 mdm                                     
#   marco[dot]masciola[at]gmail                                
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


BITS = -m64
PLATFORM = $(shell uname -s)
# CMINPACK_DIR = cminpack
# BSTRING_DIR = bstring
# SIMCLIST_DIR = simclist
# LAPACK_DIR = lapack

VPATH = cminpack:simclist:bstring:lapack


ifeq ($(PLATFORM),Darwin)
  CC_TOOLS  = clang
  CFLAGS    = -g -O1 -fsanitize=address -fno-omit-frame-pointer -fPIC -D DEBUG -Icminpack -Isimclist
  LDFLAGS   = -g -fsanitize=address -fno-omit-frame-pointer -dynamiclib 
else ifeq ($(PLATFORM),Linux)
  CC_TOOLS  = gcc
  #CFLAGS    = -g -O1 -fsanitize=address -fno-omit-frame-pointer -fPIC -D DEBUG -DGITVERSION=\"$(GIT_VERSION)\" -Icminpack -Isimclist
  #CFLAGS    = -g -fPIC -std=c99 -Wuninitialized -Wall -D DEBUG -DGITVERSION=\"$(GIT_VERSION)\" 
  CFLAGS    = $(BITS) -g -fPIC -std=c99 -D DEBUG -DGITVERSION=\"$(GIT_VERSION)\" -D WITH_LAPACK
  LDFLAGS   = $(BITS) -g -shared -llapacke
endif





GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)
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
	$(CC_TOOLS) $(LDFLAGS) -o libmap-1.10.00.so $(DEBUG) $(OBJ) -lm -llapacke

.c.o :
	$(CC_TOOLS) -c $(CFLAGS) $<

clean:
	rm -rf *.so *.o *~

memcheck:	
	valgrind --tool=memcheck --leak-check=yes python -E -tt ./main.py	
