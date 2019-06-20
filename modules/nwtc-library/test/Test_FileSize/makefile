#================================================================================#
# This makefile was last updated by M. Buhl on 26-Feb-2013.                      #
#                                                                                #
# This makefile has been tested with:                                            #
#    Windows 7                                                                   #
#       GNU Make3.81                                                             #
#       gfortran v4.6.1                                                          #
#    Ubuntu 12.04.0 LTS                                                          #
#       GNU Make3.82                                                             #
#       gfortran v4.6.3                                                          #
#                                                                                #
# It was designed to be used with:                                               #
#     Test_FileSize           (v1.00.00     , 26-Feb-2013)                       #
#     NWTC Subroutine Library (v2.00.00c-bjj, 05-Feb-2013)                       #
#================================================================================#


# WARNING!!!!!  This makefile doesn't actually work and I have no idea why.  It will not compile the Library on Windows.


   # Windows or Linux?  32-bit or 64-bit?

BITS = 32
BITS = 64

   # Location of source files for the NWTC Library.  You will probably need to change these for your system.

ifeq ($(OS),Windows_NT)
   LIB_DIR = C:\Users\mbuhl\UserData\CAEtools\Miscellaneous\NWTC_Library\branches\BJonkman\source
   #LIB_DIR = C:\Users\mbuhl\UserData\CAEtools\Miscellaneous\NWTC_Library\trunk\source
else
   LIB_DIR = $(HOME)/PC/CAEtools/Miscellaneous/NWTC_Library/branches/BJonkman/source
endif

TEST_DIR = .

   # Name of compiler to use and flags to use.

FC     = gfortran
FFLAGS = -O3 -m$(BITS) -fbacktrace -x f95-cpp-input -ffree-line-length-none
LFLAGS = -m$(BITS)

   # Precision.

# Use "SingPrec" for single precision and "DoubPrec" for double precision.  You may also need to change an option switch to make constants DP.
PREC = SingPrec

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

   # Destination for executable.

DEST_DIR = .

   # Library files.

LIB_SOURCES =         \
	$(PREC).f90        \
	$(SYS_FILE).f90    \
	NWTC_IO.f90        \
	NWTC_Num.f90       \
   ModMesh_Types.f90  \
	ModMesh.f90        \
	NWTC_Aero.f90      \
	NWTC_Library.f90

TEST_SOURCES   =      \
	Test_FileSize.f90

vpath %.f90 $(LIB_DIR) $(TEST_DIR)
vpath %.mod $(INTER_DIR)
vpath %.obj $(INTER_DIR)

LIB_OBJS  = $(LIB_SOURCES:.f90=.obj)
TEST_OBJS = $(TEST_SOURCES:.f90=.obj)

   # Rule to do everything.

all:     default
default: $(INTER_DIR) $(DEST_DIR)/Test_FileSize$(EXE_EXT)

   # General rule for making the files.

%.obj: %.f90
	$(FC) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR)

   #  Dependency rules.

$(SYS_FILE).obj:    $(PREC).obj
ModMesh_Types.obj:  $(PREC).obj
ModMesh.obj:        ModMesh_Types.obj
NWTC_IO.obj:        $(SYS_FILE).obj
NWTC_Num.obj:       NWTC_IO.obj
NWTC_Aero.obj:      NWTC_IO.obj      NWTC_Num.obj
NWTC_Library.obj:   NWTC_Aero.obj    ModMesh.obj
Test_FileSize.obj:  NWTC_Library.obj

#Test_FileSize$(EXE_EXT): Test_FileSize.obj

   # Make sure the destination directory for the intermediate files exist.

$(INTER_DIR):
	$(MD_CMD) $(INTER_DIR)

   # For compiling Test_FileSize.

$(DEST_DIR)/Test_FileSize$(EXE_EXT): $(LIB_OBJS) $(TEST_OBJS) | $(INTER_DIR)
	$(FC) $(LFLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/Test_FileSize$(EXE_EXT) \
	$(foreach src, $(LIB_OBJS), $(addprefix $(INTER_DIR)/,$(src))) \
	$(foreach src, $(TEST_OBJS), $(addprefix $(INTER_DIR)/,$(src)))

   # Cleanup afterwards.

clean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*.obj

