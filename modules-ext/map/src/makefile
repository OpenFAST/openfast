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


#BITS = -m64
#PLATFORM = $(shell uname -s)
BITS = -m32
PLATFORM = Windows

# CMINPACK_DIR = cminpack
# BSTRING_DIR = bstring
# SIMCLIST_DIR = simclist
# LAPACK_DIR = lapack

VPATH = cminpack:simclist:bstring:lapack

ifeq ($(PLATFORM),Windows)
  DEL_CMD   = del
  LIB_EXT   = dll

  CC_TOOLS  = gcc
  CFLAGS    = $(BITS) -g -std=c99 -DMAP_DLL_EXPORTS -DCMINPACK_NO_DLL  -DNDEBUG -D_WINDOWS -D_USRDLL -D_MINGW
  LDFLAGS   = $(BITS) -g -shared -Wl,--export-all-symbols
else
  LIB_EXT   = so
  DEL_CMD   = rm -rf
  GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always)

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
	$(CC_TOOLS) $(LDFLAGS) -o libmap-1.10.00.$(LIB_EXT) $(DEBUG) $(OBJ)
# -lm -llapacke

.c.o :
	$(CC_TOOLS) -c $(CFLAGS) $<

clean:
	$(DEL_CMD) *.$(LIB_EXT) *.o *~

memcheck:
	valgrind --tool=memcheck --leak-check=yes python -E -tt ./main.py
