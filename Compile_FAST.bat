@ECHO OFF

REM The calling syntax for this script is
REM                  Compile_FAST [dll]
REM
REM Add the "dll" to the command line to compile FAST for the Bladed-style dll.

REM ----------------------------------------------------------------------------
REM -------------------- LOCAL VARIABLES ---------------------------------------
REM ----------------------------------------------------------------------------

SET ROOT_NAME=FAST

SET COMPOPTS=/threads /inline:speed /optimize:5 /traceback /Qzero /Qsave /real_size:32 /assume:byterecl /nostand
rem SET LINKOPTS=/link /stack:64000000
SET LINKOPTS=/link 

REM ----------------------------------------------------------------------------
REM ------------------------- LOCAL PATHS --------------------------------------
REM ----------------------------------------------------------------------------
REM -- USERS WILL NEED TO EDIT THESE PATHS TO POINT TO FOLDERS ON THEIR LOCAL --
REM ------- NOTE: do not use quotation marks around the path names!!!! ---------
REM ----------------------------------------------------------------------------
REM IVF_LOC      is the location of the Intel Visual Fortran files (IFORTVARS.bat and IFORT.exe).
REM AeroDyn_Loc  is the location of the AeroDyn source files
REM Wind_Loc     is the location of the AeroDyn wind inflow source files
REM FAST_LOC     is the location of the FAST source files
REM NWTC_Lib_Loc is the location of the NWTC subroutine library files
REM ----------------------------------------------------------------------------

SET IVF_LOC=C:\Program Files\Intel\Compiler\Fortran\10.1.024\IA32\Bin

SET NWTC_Lib_Loc=D:\DATA\Fortran\Fortran Subs\NWTC_subs\v1.02.00\Source
SET AeroDyn_Loc=D:\DATA\Fortran\IVF Projects\AeroDyn\Release\Source
SET Wind_Loc=D:\DATA\Fortran\IVF Projects\AeroDyn\Release\Source\InflowWind\Source
SET FAST_Loc=D:\DATA\Fortran\IVF Projects\FAST\Release\Source


REM ----------------------------------------------------------------------------
REM -------------------- LIST OF ALL SOURCE FILES ------------------------------
REM ----------------------------------------------------------------------------

SET IVF="%IVF_LOC%\IFORT.EXE"


SET NWTC_Files=
SET NWTC_Files=%NWTC_Files%  "%NWTC_Lib_Loc%\SingPrec.f90"
SET NWTC_Files=%NWTC_Files%  "%NWTC_Lib_Loc%\SysIVF.f90" 
SET NWTC_Files=%NWTC_Files%  "%NWTC_Lib_Loc%\NWTC_IO.f90"
SET NWTC_Files=%NWTC_Files%  "%NWTC_Lib_Loc%\NWTC_Num.f90"
SET NWTC_Files=%NWTC_Files%  "%NWTC_Lib_Loc%\NWTC_Aero.f90"
SET NWTC_Files=%NWTC_Files%  "%NWTC_Lib_Loc%\NWTC_Library.f90"


SET Wind_Files=
SET Wind_Files=%Wind_Files%  "%Wind_Loc%\SharedInflowDefs.f90"
SET Wind_Files=%Wind_Files%  "%Wind_Loc%\HHWind.f90"
SET Wind_Files=%Wind_Files%  "%Wind_Loc%\FFWind.f90"
SET Wind_Files=%Wind_Files%  "%Wind_Loc%\FDWind.f90"
SET Wind_Files=%Wind_Files%  "%Wind_Loc%\CTWind.f90"
SET Wind_Files=%Wind_Files%  "%Wind_Loc%\UserWind.f90"
SET Wind_Files=%Wind_Files%  "%Wind_Loc%\InflowWindMod.f90"


SET AeroDyn_Files=
SET AeroDyn_Files=%AeroDyn_Files%  "%AeroDyn_Loc%\SharedTypes.f90"
SET AeroDyn_Files=%AeroDyn_Files%  "%AeroDyn_Loc%\AeroMods.f90"
SET AeroDyn_Files=%AeroDyn_Files%  "%AeroDyn_Loc%\GenSubs.f90"
SET AeroDyn_Files=%AeroDyn_Files%  "%AeroDyn_Loc%\AeroSubs.f90"
SET AeroDyn_Files=%AeroDyn_Files%  "%AeroDyn_Loc%\AeroDyn.f90"


SET FAST_Files=
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\fftpack.f"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\FFTMod.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\HydroCalc.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\FAST_Mods.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\Noise.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\FAST_IO.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\FAST.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\FAST_Lin.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\FAST2ADAMS.f90"

IF /I "%1"=="DLL" GOTO dllFiles

SET FAST_Files=%FAST_Files%  "%FAST_LOC%\PitchCntrl_ACH.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\UserSubs.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\UserVSCont_KP.f90"

GOTO endFASTfiles

:dllFiles
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\BladedDLLInterface.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\UserSubs_forBladedDLL.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\UserVSCont_KP_forBladedDLL.f90"

REM NOTE: UserSubs_forBladedDLL.f90 is a copy of UserSubs.f90 with SUBROUTINES UserHSSBr() and UserYawCont() commented out
REM       UserVSCont_KP_forBladedDLL.f90 is a copy of UserVSCont_KP.f90 with SUBROUTINE UserVSCont() commented out

SET ROOT_NAME=%ROOT_NAME%_DLL

:endFASTfiles
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\AeroCalc.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\SetVersion.f90"
SET FAST_Files=%FAST_Files%  "%FAST_LOC%\FAST_Prog.f90"



:ivf
REM ----------------------------------------------------------------------------
REM ---------------- COMPILE WITH INTEL VISUAL FORTRAN -------------------------
REM ----------------------------------------------------------------------------

REM                   set compiler internal variables 
IF "%INTEL_SHARED%"=="" CALL "%IVF_LOC%\IFORTVARS.bat"

REM                           compile 

ECHO.
ECHO Compiling FAST, AeroDyn, and NWTC_Library routines to create %ROOT_NAME%.exe:

%IVF% %COMPOPTS% %NWTC_Files% %Wind_Files% %AeroDyn_Files% %FAST_Files% %LINKOPTS% /out:%ROOT_NAME%.exe


:end
REM ----------------------------------------------------------------------------
REM ------------------------- CLEAR MEMORY -------------------------------------
REM ------------- and delete all .mod and .obj files ---------------------------
REM ----------------------------------------------------------------------------
ECHO 

DEL *.mod
DEL *.obj

SET IVF_LOC=
SET IVF=
SET ROOT_NAME=
SET COPTS=

SET NWTC_Files=
SET Wind_Files=
SET AeroDyn_Files=
SET FAST_Files=
SET A2AD_Files=
SET Fixed_Files=

SET NWTC_Lib_Loc=
SET Wind_Loc=
SET AeroDyn_Loc=
SET A2AD_Loc=
SET FAST_Loc=

SET COMPOPTS=
SET LINKOPTS=