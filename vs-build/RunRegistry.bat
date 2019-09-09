@ECHO OFF

set lines=======================================================================
echo %lines%
IF "%1"=="" (
ECHO.
ECHO   The calling syntax for this script is
ECHO             RunRegistry ModuleName [FAST_Root_Loc]
ECHO.
GOTO Done
)


REM ----------------------------------------------------------------------------
REM ------------------------- LOCAL PATHS --------------------------------------
REM ----------------------------------------------------------------------------
REM -- USERS MAY EDIT THESE PATHS TO POINT TO FOLDERS ON THEIR LOCAL MACHINES. -
REM -- NOTE: do not use quotation marks around the path names!!!! --------------
REM ----------------------------------------------------------------------------
REM ----------------------------------------------------------------------------
SET Root_Loc=..\..
IF not "%2"=="" SET Root_Loc=%2

SET Local_Modules_Loc=%Root_Loc%\modules
SET Ext_Modules_Loc=%Root_Loc%\modules
SET Registry=..\..\build\bin\Registry.exe
SET Output_Loc=..\..\build\types-files
SET FAST_Loc=%Local_Modules_Loc%\openfast-library\src
SET ED_Loc=%Local_Modules_Loc%\elastodyn\src
SET AD14_Loc=%Local_Modules_Loc%\aerodyn14\src
SET IfW_Loc=%Local_Modules_Loc%\inflowwind\src
SET HD_Loc=%Local_Modules_Loc%\hydrodyn\src
SET SD_Loc=%Local_Modules_Loc%\subdyn\src
SET MAP_Loc=%Ext_Modules_Loc%\map
SET FEAM_Loc=%Ext_Modules_Loc%\feamooring\src
SET IceF_Loc=%Ext_Modules_Loc%\icefloe\src\interfaces\FAST
SET IceD_Loc=%Ext_Modules_Loc%\icedyn\src
SET MD_Loc=%Ext_Modules_Loc%\moordyn\src
SET OpFM_Loc=%Local_Modules_Loc%\openfoam\src
SET Orca_Loc=%Local_Modules_Loc%\orcaflex-interface\src
SET NWTC_Lib_Loc=%Local_Modules_Loc%\nwtc-library\src
SET ExtPtfm_Loc=%Local_Modules_Loc%\extptfm\src
SET AD_Loc=%Local_Modules_Loc%\aerodyn\src
SET SrvD_Loc=%Local_Modules_Loc%\servodyn\src
SET BD_Loc=%Local_Modules_Loc%\beamdyn\src
SET SC_Loc=%Local_Modules_Loc%\supercontroller\src

SET MAP_Loc_R=%MAP_Loc%\src

SET ALL_FAST_Includes=-I "%FAST_Loc%" -I "%NWTC_Lib_Loc%" -I "%ED_Loc%" -I "%SrvD_Loc%" -I "%AD14_Loc%" -I^
 "%AD_Loc%" -I "%BD_Loc%" -I "%SC_Loc%" -I^
 "%IfW_Loc%" -I "%SD_Loc%" -I "%HD_Loc%" -I "%MAP_Loc_R%" -I "%FEAM_Loc%"  -I^
 "%IceF_Loc%" -I "%IceD_Loc%" -I "%MD_Loc%" -I "%OpFM_Loc%" -I "%Orca_Loc%" -I "%ExtPtfm_Loc%"


SET ModuleName=%1

GOTO %ModuleName%

REM ----------------------------------------------------------------------------
REM ---------------- RUN THE REGISTRY TO AUTO-GENERATE FILES -------------------
REM ----------------------------------------------------------------------------
:MAP
SET CURR_LOC=%MAP_Loc_R%
%REGISTRY% "%CURR_LOC%\%ModuleName%_Registry.txt"  -ccode -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%"
%REGISTRY% "%CURR_LOC%\MAP_Fortran_Registry.txt"  -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%" -noextrap
GOTO checkError

:FAST
SET CURR_LOC=%FAST_Loc%
%REGISTRY% "%CURR_LOC%\FAST_Registry.txt" %ALL_FAST_Includes% -noextrap -O "%Output_Loc%"
GOTO checkError

:BeamDyn
SET CURR_LOC=%BD_Loc%
%REGISTRY% "%CURR_LOC%\Registry_BeamDyn.txt" -I "%NWTC_Lib_Loc%" -O "%Output_Loc%"
GOTO checkError

:SuperController
SET CURR_LOC=%SC_Loc%
%REGISTRY% "%CURR_LOC%\SuperController_Registry.txt" -I "%NWTC_Lib_Loc%" -O "%Output_Loc%" -ccode -noextrap
GOTO checkError

:ElastoDyn
SET CURR_LOC=%ED_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%_Registry.txt" -I "%NWTC_Lib_Loc%" -O "%Output_Loc%"
GOTO checkError

:TMD
:ServoDyn
SET CURR_LOC=%SrvD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%_Registry.txt" -I "%NWTC_Lib_Loc%" -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:Lidar
:InflowWind
SET CURR_LOC=%IfW_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%" -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:IfW_TSFFWind
:IfW_HAWCWind
:IfW_BladedFFWind
:IfW_UserWind
:IfW_4Dext
:IfW_UniformWind
SET CURR_LOC=%IfW_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%" -I "%CURR_LOC%" -noextrap  -O "%Output_Loc%"
GOTO checkError

:OpenFOAM
SET CURR_LOC=%OpFM_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%_Registry.txt" -I "%NWTC_Lib_Loc%" -ccode -O "%Output_Loc%"
GOTO checkError

:AeroDyn
:BEMT
:DBEMT
SET CURR_LOC=%AD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%_Registry.txt" -I "%NWTC_Lib_Loc%" -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:AeroDyn_Driver
SET CURR_LOC=%AD_Loc%
%REGISTRY% "%CURR_LOC%\AeroDyn_Driver_Registry.txt" -I %NWTC_Lib_Loc% -I %CURR_LOC%  -O %Output_Loc% -noextrap
GOTO checkError

:AFI
SET CURR_LOC=%AD_Loc%
%REGISTRY% "%CURR_LOC%\AirfoilInfo_Registry.txt" -I "%NWTC_Lib_Loc%" -I "%CURR_LOC%" -noextrap -O "%Output_Loc%"
GOTO checkError

:UA
SET CURR_LOC=%AD_Loc%
%REGISTRY% "%CURR_LOC%\UnsteadyAero_Registry.txt" -I "%NWTC_Lib_Loc%" -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:AeroDyn14
SET CURR_LOC=%AD14_Loc%
%REGISTRY% "%CURR_LOC%\Registry-AD14.txt" -I "%NWTC_Lib_Loc%" -I "%CURR_LOC%" -I "%IfW_Loc%" -O "%Output_Loc%"
GOTO checkError

:DWM
SET CURR_LOC=%AD14_Loc%
%REGISTRY% "%CURR_LOC%\Registry-DWM.txt" -I "%NWTC_Lib_Loc%" -I "%IfW_Loc%"  -O "%Output_Loc%"
GOTO checkError

:HydroDyn
SET CURR_LOC=%HD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:Current
SET CURR_LOC=%HD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:Waves
SET CURR_LOC=%HD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:Waves2
SET CURR_LOC=%HD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:SS_Excitation
SET CURR_LOC=%HD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:SS_Radiation
SET CURR_LOC=%HD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:Conv_Radiation
SET CURR_LOC=%HD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:WAMIT
SET CURR_LOC=%HD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:WAMIT2
SET CURR_LOC=%HD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:Morison
SET CURR_LOC=%HD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%"  -I "%CURR_LOC%" -O "%Output_Loc%"
GOTO checkError

:SubDyn
SET CURR_LOC=%SD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%_Registry.txt" -I "%NWTC_Lib_Loc%"  -O "%Output_Loc%"
GOTO checkError

:FEAMooring
SET CURR_LOC=%FEAM_Loc%
%REGISTRY% "%CURR_LOC%\FEAM_Registry.txt" -I "%NWTC_Lib_Loc%"  -O "%Output_Loc%"
GOTO checkError

:MoorDyn
SET CURR_LOC=%MD_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%_Registry.txt" -I "%NWTC_Lib_Loc%"  -O "%Output_Loc%"
GOTO checkError

:IceFloe
SET CURR_LOC=%IceF_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%_FASTRegistry.inp" -I "%NWTC_Lib_Loc%"  -O "%Output_Loc%"
GOTO checkError

:IceDyn
SET CURR_LOC=%IceD_Loc%
%REGISTRY% "%CURR_LOC%\Registry_%ModuleName%.txt" -I "%NWTC_Lib_Loc%"  -O "%Output_Loc%"
GOTO checkError

:OrcaFlexInterface
SET CURR_LOC=%Orca_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%.txt" -I "%NWTC_Lib_Loc%" -O "%Output_Loc%"
GOTO checkError

:ExtPtfm_MCKF
SET CURR_LOC=%ExtPtfm_Loc%
%REGISTRY% "%CURR_LOC%\%ModuleName%_Registry.txt" -I "%NWTC_Lib_Loc%" -O "%Output_Loc%"
GOTO checkError

:FASTWrapper
SET CURR_LOC=%Farm_Loc%
%REGISTRY% "%CURR_LOC%\FASTWrapper_Registry.txt" -I %NWTC_Lib_Loc%  %ALL_FAST_INCLUDES% -noextrap -O "%Output_Loc%"
GOTO checkError

:checkError
ECHO.
IF %ERRORLEVEL% NEQ 0 (
ECHO Error running FAST Registry for %ModuleName%.
) ELSE (
ECHO Registry for %ModuleName% completed.
REM COPY /Y "%ModuleName%_Types.f90"   "%CURR_LOC%"
rem IF /I "%ModuleName%"=="MAP" COPY /Y "%ModuleName%_Types.h" "%CURR_LOC%"
)




:end
REM ----------------------------------------------------------------------------
REM ------------------------- CLEAR MEMORY -------------------------------------
REM ----------------------------------------------------------------------------
ECHO. 

SET ModuleName=
SET CURR_LOC=

SET Root_Loc=

SET Subs_Loc=
SET FAST_Loc=
SET Registry=

SET ED_Loc=
SET BD_Loc=
SET AD14_Loc=
SET IfW_Loc=
SET HD_Loc=
SET SD_Loc=
SET MAP_Loc=
SET FEAM_Loc=
SET IceF_Loc=
SET IceD_Loc=
SET MD_Loc=
SET OpFM_Loc=
SET Orca_Loc=
SET NWTC_Lib_Loc=
SET ExtPtfm_Loc=
SET AD_Loc=
SET SrvD_Loc=

SET MAP_Loc_R=
SET ALL_FAST_Includes=

:Done
echo %lines%
set lines=

:PathsOnly