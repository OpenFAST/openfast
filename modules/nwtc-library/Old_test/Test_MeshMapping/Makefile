#================================================================================#
# This makefile copied from B. Jonkman by J. Michalakes on 29-Jan-2013,          #
# adapted from Crunch (M. Buhl on 25-Jan-2013).                                  #
#                                                                                #
# This makefile has been tested on Windows 7 with gfortran.                      #
# This makefile works with mingw32-make.exe.                                     #
#                                                                                #
# It was designed to be used with:                                               #
#     Module1                 (v1.00.04,      29-Jan-2013)                       #
#     Module2                 (v1.00.04,      29-Jan-2013)                       #
#     NWTC Subroutine Library (CSC version    29-Jan-2013)                       #
#                                                                                #
# Older versions of ModuleName and NWTC Library may not work with this makefile. #
#================================================================================#

   # 32-bit or 64-bit?
BITS = 32
#BITS = 64


   # Location of source files for ModuleName and the NWTC Library.  You will probably need to change these for your system.
#OS = Linux

# these settings are overridden below
ifeq ($(OS),Windows_NT)
   LIB_DIR  = C:/Users/bjonkman/Documents/DATA/DesignCodes/miscellaneous/nwtc_subs/SVNdirectory/trunk/source
   Work_DIR = C:/Users/bjonkman/Documents/DATA/DesignCodes/miscellaneous/nwtc_subs/SVNdirectory/branches/MeshMapping
   REGISTRY = Registry
   LAPACK_DIR   = C:/LAPACK/win32
   NETLIB_DIR   = C:/Users/bjonkman/Documents/DATA/DesignCodes/Simulators/FAST/SVNdirectory/branches/BJonkman/source/dependencies/NetLib
else
#   LIB_DIR = $(HOME)/PC/CAEtools/Miscellaneous/NWTC_Library/trunk/source
   LIB_DIR  = ../NWTC_Library/source
   REGISTRY = ../FAST/Registry/Source/registry.exe
endif

MODNAME_DIR = .



   # Name of compiler to use and flags to use.

FC     = gfortran
FFLAGS = -O3 -m$(BITS) -fbacktrace -ffree-line-length-none -x f95-cpp-input -ffree-line-length-none
FFLAGS = -O3 -fbacktrace -ffree-line-length-none -x f95-cpp-input  -DMESH_DEBUG
LDFLAGS = -O3 -fbacktrace -llapack -lblas -L$(LAPACK_DIR)

   # Precision.

# Use "SingPrec" for single precision and "DoubPrec" for double precision.  You may also need to change an option switch to make constants DP.
PREC = SingPrec
#PREC = DoubPrec

   #==========================================================#
   # You should not need to change anything beyond this point #
   #==========================================================#

   # System-specific settings.

ifeq ($(OS),Windows_NT)
      # Windows
   DEL_CMD   = del
   EXE_EXT   = _gwin$(BITS).exe
   INTER_DIR = Obj_win$(BITS)
   MD_CMD    = @mkdir
   OBJ_EXT   = .obj
   PATH_SEP  = \\
   SYS_FILE  = SysGnuWin
else
      # Linux
   DEL_CMD   = rm -f
   EXE_EXT   = _glin$(BITS)
   INTER_DIR = Obj_lin$(BITS)
   MD_CMD    = @mkdir -p
   OBJ_EXT   = .o
   PATH_SEP  = /
   SYS_FILE  = SysGnuLinux
endif

   # Destination and RootName for executable

OUTPUT_NAME = Test_MeshMapping
DEST_DIR    = .

   # Library files.

LIB_SOURCES =                 \
	$(PREC).f90                \
	NWTC_Base.f90              \
	$(SYS_FILE).f90            \
	NWTC_IO.f90                \
	NWTC_Library_Types.f90     \
	ModMesh_Types.f90          \
	ModMesh.f90                \
	ModMesh_Mapping.f90        \
	NWTC_Library.f90           \
	NWTC_Num.f90

NETLIB_SOURCES =           \
	NWTC_ScaLAPACK.f90      \
	NWTC_LAPACK.f90         \
	DLASRT2.f               \
	SLASRT2.f


MODNAME_SOURCES   =             \
	Test_MeshMapping.f90

vpath %.f90 $(Work_DIR) $(LIB_DIR) $(NETLIB_DIR) $(MODNAME_DIR)
vpath %.mod $(INTER_DIR)
vpath %.obj $(INTER_DIR)


LIB_OBJS    = $(LIB_SOURCES:.f90=.obj)
NETLIB_OBJS   = NWTC_LAPACK.obj
#NWTC_ScaLAPACK.obj DLASRT2.obj SLASRT2.obj
MODNAME_OBJS = $(MODNAME_SOURCES:.f90=.obj)

   # Rule to do everything.

all:     default
default: Test_MeshMapping.f90 $(INTER_DIR) $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT)

   # General rule for making the files.

%.obj: %.f90
	$(FC) -I $(INTER_DIR) $(FFLAGS) -g -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR)

%.obj: %.f
	$(F77) -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)


   #  Dependency rules.

NWTC_LAPACK.obj:         NWTC_Base.obj
NWTC_ScaLAPACK.obj:      NWTC_Base.obj DLASRT2.obj SLASRT2.obj

NWTC_Base.obj:           $(PREC).obj
$(SYS_FILE).obj:         NWTC_Base.obj
NWTC_Library_Types.obj:  $(SYS_FILE).obj
NWTC_IO.obj:             NWTC_Library_Types.obj
NWTC_Num.obj:            NWTC_IO.obj
ModMesh_Types.obj:       NWTC_Num.obj
ModMesh.obj:             ModMesh_Types.obj
ModMesh_Mapping.obj:     ModMesh.obj  NWTC_LAPACK.obj
NWTC_Library.obj:        ModMesh.obj  ModMesh_Mapping.obj


Test_MeshMapping.obj:      ModMesh_Mapping.obj      NWTC_Library.obj


#$(OUTPUT_NAME)$(EXE_EXT): Test_MeshMapping.obj

   # Make sure the destination directory for the intermediate files exist.

$(INTER_DIR):
	$(MD_CMD) $(INTER_DIR)



   # For compiling the driver/glue code.

$(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT): $(LIB_OBJS) $(MODNAME_OBJS) | $(INTER_DIR)
	$(FC) $(LDFLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT) \
	$(foreach src, $(NETLIB_OBJS), $(addprefix $(INTER_DIR)/,$(src))) \
	$(foreach src, $(LIB_OBJS), $(addprefix $(INTER_DIR)/,$(src))) \
	$(foreach src, $(MODNAME_OBJS), $(addprefix $(INTER_DIR)/,$(src)))

   # Cleanup afterwards.

clean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*.obj $(OUTPUT_NAME)$(EXE_EXT)

