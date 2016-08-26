#================================================================================#
# This makefile created by B. Jonkman.                                           #
# (c) 2013-2015 National Renewable Energy Laboratory                             #
#                                                                                #
# It has been tested on                                                          #
#     Windows 7 with gfortran using mingw32-make.exe.                            #
#                                                                                #
#================================================================================#


   # Name of compiler to use and flags to use.
FC      = gfortran

# 32-bit or 64-bit?
BITS = 32
#BITS = 64


FFLAGS  = -O2 -m$(BITS) -fbacktrace -ffree-line-length-none -x f95-cpp-input
LDFLAGS = -O2 -m$(BITS) -fbacktrace


#.........................................
# preprocessor definitions:
#  -DDOUBLE_PRECISION
#
# some useful gfortran options:
# -Wl,--stack=999999999,--large-address-aware
#
#  -DFPE_TRAP_ENABLED
#  -Wl,--large-address-aware              # not necessary when $(BITS) is 64 (64-bit target)
#  -Wl,--stack=999999999                  # not necessary when $(BITS) is 64
#  -Wconversion-extra -Wconversion
#  -fdefault-real-8                       # note this also promotes doubles to 16-bytes if possible
#  -fcheck=bounds,do,mem,pointer -std=f2003 -O0 -v -Wall
#  -pg                                    # generate profiling and debugging info
# http://gcc.gnu.org/onlinedocs/gfortran/Option-Index.html#Option-Index
#
# makefile options:
#  --warn-undefined-variables
#.........................................


   # Location of source files for FAST, and its modules.

ModName_DIR  = ../source
Driver_DIR   = ../source/drivers
BIN_DIR      = ../bin

NWTC_Lib_DIR = ../subs/NWTC_Library/source
NETLIB_DIR   = $(NWTC_Lib_Dir)/NetLib/lapack


   # Names and locations of the Registry, MAP libraries, and instructions for linking with LAPACK.
   # You will probably need to change these for your system.

ifeq ($(OS),Windows_NT)

   Registry     = $(BIN_DIR)/Registry_win32.exe
   LDFLAGS     := $(LDFLAGS) -llapack -lblas -LC:/LAPACK/win32
   LAPACK_LINK :=

   UNAME       := $(OS)

else

   Registry     = ../bin/registry.exe

   UNAME := $(shell uname -s)

   ifeq ($(UNAME), Darwin)
         # Mac needs some additional flags prepended for linking to get lapack working:
      LDFLAGS     := -framework Accelerate $(LDFLAGS)
      LAPACK_LINK :=
      #DLL_EXT     := dylib
   else
        # ifeq ($(UNAME), Linux)
        # bjj verify this: # I've been told that the lapack linking needs to occur at the end of the line (after all the object files) on linux:
      LAPACK_LINK  = -llapack -lblas
      #DLL_EXT     := so
   endif

endif


   # Destination and RootName for executable

DEST_DIR    = $(BIN_DIR)
OUTPUT_NAME = ModName_Driver


   #==========================================================#
   # You should not need to change anything beyond this point #
   #==========================================================#

   # System-specific settings.

ifeq ($(OS),Windows_NT)

      # Windows
   DEL_CMD   = del
   INTER_DIR = Obj_win$(BITS)
   MD_CMD    = @mkdir
   CP_CMD    = copy
   OBJ_EXT   = .obj
   PATH_SEP  = \\
   SYS_FILE  = SysGnuWin

   ifeq ($(FAST_driver),FAST_Library)
      EXE_EXT   = _gwin$(BITS).dll
   else
      EXE_EXT   = _gwin$(BITS).exe
   endif

else

      # Linux/Mac
   DEL_CMD   = rm -f
   INTER_DIR = Obj_lin$(BITS)
   MD_CMD    = @mkdir -p
   CP_CMD    = cp
   OBJ_EXT   = .o
   PATH_SEP  = /
   SYS_FILE  = SysGnuLinux

   ifeq ($(FAST_driver),FAST_Library)
      EXE_EXT   = _glin$(BITS).so
   else
      EXE_EXT   = _glin$(BITS)
   endif

endif


   # Source files (by module)

LIB_SOURCES =                  \
	SingPrec.f90                \
	NWTC_Base.f90               \
	$(SYS_FILE).f90             \
	NWTC_Library_Types.f90      \
	NWTC_IO.f90                 \
	NWTC_Num.f90                \
	ModMesh_Types.f90           \
	ModMesh.f90                 \
	ModMesh_Mapping.f90         \
	NWTC_Library.f90

NETLIB_SOURCES=                \
	NWTC_LAPACK.f90

MODNAME_SOURCES   =            \
	ModuleName.f90              \
	ModuleName_Types.f90        \
	ModName_Driver.f90



vpath %.f90       $(NWTC_Lib_DIR) $(NETLIB_DIR) $(ModName_DIR) $(Driver_DIR)
vpath %.mod       $(INTER_DIR)
vpath %$(OBJ_EXT) $(INTER_DIR)
vpath %.txt       $(ModName_DIR)


ALL_SOURCES = $(MODNAME_SOURCES) $(LIB_SOURCES) $(NETLIB_SOURCES)
ALL_OBJS    = $(ALL_SOURCES:.f90=$(OBJ_EXT))


   # Rule to do everything.

all:     default
# use this for "all" for debugging: print variables:
#all:     ; $(info $$IceF_SOURCES is [${IceF_SOURCES}]) echo debugging
#
default:  $(INTER_DIR) $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT)


   # General rules for compiling the files.

%$(OBJ_EXT): %.f90
	$(FC) -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)



   # General rules for creating _Types.f90 files from Registry files.

$(ModName_DIR)/ModuleName_Types.f90: $(ModName_DIR)/ModuleName_Registry.txt
	$(Registry) $< -I $(NWTC_Lib_DIR)  -O $(ModName_DIR)



   #  Dependency rules.
#NWTC Library dependency rules:
NWTC_Base$(OBJ_EXT):              SingPrec$(OBJ_EXT)
$(SYS_FILE)$(OBJ_EXT):            NWTC_Base$(OBJ_EXT)
NWTC_Library_Types$(OBJ_EXT):     $(SYS_FILE)$(OBJ_EXT)
NWTC_IO$(OBJ_EXT):                NWTC_Library_Types$(OBJ_EXT)
NWTC_Num$(OBJ_EXT):               NWTC_IO$(OBJ_EXT)
ModMesh_Types$(OBJ_EXT):          NWTC_Num$(OBJ_EXT)
ModMesh$(OBJ_EXT):                ModMesh_Types$(OBJ_EXT)
ModMesh_Mapping$(OBJ_EXT):        ModMesh$(OBJ_EXT) NWTC_LAPACK$(OBJ_EXT)
NWTC_Library$(OBJ_EXT):           ModMesh$(OBJ_EXT)  ModMesh_Mapping$(OBJ_EXT)

NWTC_LAPACK$(OBJ_EXT):            NWTC_Base$(OBJ_EXT)


ModuleName_Types$(OBJ_EXT):       NWTC_Library$(OBJ_EXT) $(ModName_DIR)/ModuleName_Types.f90
ModuleName$(OBJ_EXT):             ModuleName_Types$(OBJ_EXT)
ModName_Driver$(OBJ_EXT):         ModuleName$(OBJ_EXT)   NWTC_Library$(OBJ_EXT)



   # Make sure the destination directory for the intermediate files exist.

$(INTER_DIR):
	$(MD_CMD) $(INTER_DIR)


   # For linking

$(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT): $(ALL_OBJS) | $(INTER_DIR)
	$(FC) $(LDFLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT) \
	$(foreach src, $(ALL_OBJS), $(addprefix $(INTER_DIR)/,$(src))) $(LAPACK_LINK)

   # Cleanup afterwards

clean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*$(OBJ_EXT)

superclean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*$(OBJ_EXT)
	$(DEL_CMD) "$(DEST_DIR)$(PATH_SEP)$(OUTPUT_NAME)$(EXE_EXT)"
	$(DEL_CMD) "$(ModName_DIR)$(PATH_SEP)ModuleName_Types.f90"
	$(DEL_CMD) *_Types.f90



