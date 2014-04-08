#================================================================================#
# This makefile created by B. Jonkman on 7-Apr-2014,                             #
# adapted from Crunch (M. Buhl on 25-Jan-2013).                                  #
# (c) 2013 National Renewable Energy Laboratory                                  #
#                                                                                #
# This makefile has been tested on Windows 7 with gfortran.                      #
# This makefile works with mingw32-make.exe.                                     #
#                                                                                #
# It was designed to be used with discon.f90 for the Bladed DLL Interface        #
#================================================================================#

   # 32-bit or 64-bit?
BITS = 32
#BITS = 64


   # Location of source files for FAST, AeroDyn, InflowWind, and the NWTC Library.
   # You will probably need to change these for your system.


FAST_DIR     = ../Source
DLL_DIR      = .


ifeq ($(OS),Windows_NT)
#override FAST location on Bonnie's PC:
   FAST_DIR    = C:/Users/bjonkman/Documents/DATA/DesignCodes/simulators/FAST/SVNdirectory/branches/BJonkman/Source
   LAPACK_LINK  = -llapack -lblas -LC:/LAPACK/win32
else
   LAPACK_LINK  = -llapack -lblas
endif

NWTC_Lib_DIR = $(FAST_DIR)/dependencies/NWTC_Library
NETLIB_DIR   = $(FAST_DIR)/dependencies/NetLib
MAP_DIR      = $(FAST_DIR)/dependencies/MAP


   # Name of compiler to use and flags to use.
FC      = gfortran
FFLAGS  = -O2 -m$(BITS) -fbacktrace -ffree-line-length-none -x f95-cpp-input -C
LDFLAGS = -shared -O2 -m$(BITS) -fbacktrace $(LAPACK_LINK)



   # Destination and RootName for executable

OUTPUT_NAME = MAP
DEST_DIR    = .

PREC = SingPrec


   #==========================================================#
   # You should not need to change anything beyond this point #
   #==========================================================#

   # System-specific settings.

ifeq ($(OS),Windows_NT)
      # Windows
   DEL_CMD   = del
   EXE_EXT   = _gwin$(BITS).dll
   INTER_DIR = Obj_win$(BITS)
   MD_CMD    = @mkdir
   OBJ_EXT   = .obj
   PATH_SEP  = \\
   SYS_FILE  = SysGnuWin
else
      # Linux
   DEL_CMD   = rm -f
   EXE_EXT   = _glin$(BITS).so
   INTER_DIR = Obj_lin$(BITS)
   MD_CMD    = @mkdir -p
   OBJ_EXT   = .o
   PATH_SEP  = /
   SYS_FILE  = SysGnuLinux
endif



   # Source files (by module)

DLL_SOURCES =          \
	MAP_dll.f90         \
	MAP_Types_Intf.f90  \
	MAP_C_Types.f90     \
	MAP_Types.f90

NETLIB_SOURCES=            \
	NWTC_ScaLAPACK.f90      \
	NWTC_LAPACK.f90         \
	dlasrt2.f               \
	slasrt2.f

LIB_SOURCES =           \
	$(PREC).f90          \
	NWTC_Base.f90        \
	$(SYS_FILE).f90      \
	NWTC_Library_Types.f90 \
	NWTC_IO.f90          \
	NWTC_Num.f90         \
	ModMesh_Types.f90    \
	ModMesh.f90          \
	ModMesh_Mapping.f90  \
	NWTC_Library.f90



vpath %.f90 $(NWTC_Lib_DIR) $(MAP_DIR) $(NETLIB_DIR) $(DLL_DIR)
vpath %.f   $(NETLIB_DIR)
vpath %.mod $(INTER_DIR)
vpath %.obj $(INTER_DIR)


ALL_SOURCES = $(LIB_SOURCES) $(NETLIB_SOURCES) $(DLL_SOURCES)
tmp_objs1   = $(ALL_SOURCES:.f90=.obj)
ALL_OBJS    = $(tmp_objs1:.f=.obj)


   # Rule to do everything.

all:     default
#all:     ; $(info $$ALL_OBJS is [${ALL_OBJS}]) echo debugging

default:  $(INTER_DIR) $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT)

   # General rules for compliling the files.

%.obj: %.f90
	$(FC) -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)

%.obj: %.f
	$(F77) -I $(INTER_DIR) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR) -B $(INTER_DIR)


   #  Dependency rules.
#NWTC Library dependency rules:
NWTC_Base.obj:           $(PREC).obj
$(SYS_FILE).obj:         NWTC_Base.obj
NWTC_Library_Types.obj:  $(SYS_FILE).obj
NWTC_IO.obj:             NWTC_Library_Types.obj
NWTC_Num.obj:            NWTC_IO.obj
ModMesh_Types.obj:       NWTC_Num.obj
ModMesh.obj:             ModMesh_Types.obj
ModMesh_Mapping.obj:     ModMesh.obj NWTC_LAPACK.obj
NWTC_Library.obj:        ModMesh.obj  ModMesh_Mapping.obj

NWTC_LAPACK.obj:         NWTC_Base.obj
NWTC_ScaLAPACK.obj:      NWTC_Base.obj DLASRT2.obj SLASRT2.obj

# MAP dependency rules:
MAP_C_Types.obj:         NWTC_Library.obj
MAP_Types.obj:           NWTC_Library.obj MAP_C_Types.obj
MAP_Types_Intf.obj:       NWTC_Library.obj MAP_Types.obj
MAP_dll.obj:             MAP_Types.obj

   # Make sure the destination directory for the intermediate files exist.

$(INTER_DIR):
	$(MD_CMD) $(INTER_DIR)


   # For linking DLL.

$(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT): $(ALL_OBJS) | $(INTER_DIR)
	$(FC) $(LDFLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/$(OUTPUT_NAME)$(EXE_EXT) \
	$(foreach src, $(ALL_OBJS), $(addprefix $(INTER_DIR)/,$(src)))

   # Cleanup afterwards.

clean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*.obj

#gfortran -shared -O2 -m32 -fbacktrace -ffree-line-length-none -x f95-cpp-input -C ../CertTest/5MW_Baseline/ServoData/Source/Discon.f90