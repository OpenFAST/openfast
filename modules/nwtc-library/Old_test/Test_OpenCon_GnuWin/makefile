#==========================================================================================#
# This makefile created by M. Buhl on 26-Dec-2012.                                         #
#                                                                                          #
# This makefile has been tested on Windows 7 with gfortran.                                #
# This makefile works with mingw32-make.exe.                                               #
#                                                                                          #
# It was designed to be used with:                                                         #
#     Test_OpenCon_Gnu        (v1.00.00, 26-Dec-2012)                                      #
#     NWTC Subroutine Library (v1.06.00c-mlb, 12-Dec-2012)                                 #
#                                                                                          #
# Older versions of Test_OpenCon_Gnu and the NWTC Library may not work with this makefile. #
#==========================================================================================#

   # Location of source files for Test_OpenCon_Gnu and the NWTC Library.  You will probably need to change these for your system.

LIB_DIR  = D:/CAEtools/Miscellaneous/NWTC_Library/trunk/source
TEST_DIR = .

   # System-specific settings.

   # Windows
DEL_CMD  = del
MD_CMD   = mkdir
PATH_SEP = \\
EXE_EXT  = .exe
OBJ_EXT  = .obj
SYS_FILE = SysGnuWin

   # Linux
#DEL_CMD  = rm -f
#MD_CMD   = mkdir -p
#PATH_SEP = /
#EXE_EXT  =
#OBJ_EXT  = .o
#SYS_FILE = SysGnuLinux

   # Program name.

PROG = Test_OpenCon_GnuWin

   # Name of compiler to use and flags to use.

FC     = gfortran
FFLAGS = -O3 -fbacktrace

   # Destinations for executable and intermediate files.

DEST_DIR  = .
INTER_DIR = Obj

   # Precision.

# Change to "DoubPrec" for double precision.  You may also need to change an option switch to make constants DP.
PREC = SingPrec

   # Library files.

LIB_SOURCES =         \
	$(PREC).f90        \
	$(SYS_FILE).f90    \
	NWTC_IO.f90        \
	NWTC_Num.f90       \
	ModMesh.f90        \
	NWTC_Aero.f90      \
	NWTC_Library.f90

   #==========================================================#
   # You should not need to change anything beyond this point #
   #==========================================================#

TEST_SOURCES = $(PROG).f90

vpath %.f90 $(LIB_DIR) $(TEST_DIR)
vpath %.mod $(INTER_DIR)
vpath %.obj $(INTER_DIR)

LIB_OBJS  = $(LIB_SOURCES:.f90=.obj)
TEST_OBJS = $(TEST_SOURCES:.f90=.obj)

   # Rule to do everything.

all:     default
default: $(INTER_DIR) $(DEST_DIR)/$(PROG)$(EXE_EXT)

   # General rule for making the files.

%.obj: %.f90
	$(FC) $(FFLAGS) -c $< -o $(INTER_DIR)/$@ -J $(INTER_DIR)

   #  Dependency rules.

$(SYS_FILE).obj:  $(PREC).obj
ModMesh.obj:      $(PREC).obj
NWTC_IO.obj:      $(SYS_FILE).obj
NWTC_Num.obj:     NWTC_IO.obj
NWTC_Aero.obj:    NWTC_IO.obj      NWTC_Num.obj
NWTC_Library.obj: NWTC_Aero.obj    ModMesh.obj
$(PROG).obj:      NWTC_Library.obj
$(PROG).exe:      $(PROG).obj

   # Make sure the destination directory for the intermediate files exist.

$(INTER_DIR):
	$(MD_CMD) $(INTER_DIR)

   # For compiling Test_OpenCon_Gnu.

$(DEST_DIR)/$(PROG)$(EXE_EXT): $(LIB_OBJS) $(TEST_OBJS) | $(INTER_DIR)
	$(FC) $(FLAGS) -I $(INTER_DIR) -o $(DEST_DIR)/$(PROG)$(EXE_EXT) \
	$(foreach src, $(LIB_OBJS), $(addprefix $(INTER_DIR)/,$(src))) \
	$(foreach src, $(TEST_OBJS), $(addprefix $(INTER_DIR)/,$(src)))

   # Cleanup afterwards.

clean:
	$(DEL_CMD) $(INTER_DIR)$(PATH_SEP)*.mod $(INTER_DIR)$(PATH_SEP)*.obj

