@ECHO OFF


::=======================================================================================================
:: run a script to set all of the path locations needed to copy files

CALL .\Compiling\Set_FAST_Paths.bat

::=======================================================================================================

SET bin_dir=.\bin
SET depend_dir=.\Source\dependencies
SET SW_ModuleOnly=-MODULE

IF /I "%1"=="%SW_ModuleOnly%" goto %2

:BinDir
ECHO Binaries
COPY "%CRUNCH%"                      "%bin_dir%"

IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars

:dependencies

:Registry
ECHO FAST Registry
COPY "%Registry%"        "%bin_dir%\Registry_Win32.exe"

SET src_folder=%REG_Loc%
SET dst_folder=%depend_dir%\Registry
SET list_of_files=%src_folder%\SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars

:NWTC_Library
ECHO NWTC_Library
SET src_folder=%NWTC_Lib_Loc%\..
SET dst_folder=%depend_dir%\NWTC_Library
SET list_of_files=%src_folder%\SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars

:NetLib
ECHO NetLib
SET src_folder=%NETLIB_Loc%\..
SET dst_folder=%depend_dir%\NetLib
SET list_of_files=%src_folder%\SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars

:ElastoDyn
ECHO ElastoDyn
SET src_folder=%ED_Loc%\..
SET dst_folder=%depend_dir%\ElastoDyn
SET list_of_files=%src_folder%\FAST_SourceFiles_ED.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars

:BeamDyn
ECHO BeamDyn
SET src_folder=%BD_Loc%\..
SET dst_folder=%depend_dir%\BeamDyn
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars

:ServoDyn
:TMD
ECHO ServoDyn and TMD
SET src_folder=%SrvD_Loc%\..
SET dst_folder=%depend_dir%\ServoDyn
SET list_of_files=%src_folder%/FAST_SourceFiles_SrvD.txt

@CALL :CopyFileList

SET src_folder=%TMD_Loc%\..
SET list_of_files=%src_folder%/FAST_SourceFiles.txt
for /f %%i in (%list_of_files%) DO copy /Y "%src_folder%\%%i" "%dst_folder%"

IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars


:InflowWind
ECHO InflowWind
SET src_folder=%IfW_Loc%\..
SET dst_folder=%depend_dir%\InflowWind
SET list_of_files=%src_folder%/FAST_SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars


:AeroDyn14
ECHO AeroDyn14
SET src_folder=%AD14_Loc%\..
SET dst_folder=%depend_dir%\AeroDyn14
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars


:AeroDyn
ECHO AeroDyn
SET src_folder=%AD_Loc%\..
SET dst_folder=%depend_dir%\AeroDyn
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars


:HydroDyn
ECHO HydroDyn
SET src_folder=%HD_Loc%\..
SET dst_folder=%depend_dir%\HydroDyn
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars

:SubDyn
ECHO SubDyn
SET src_folder=%SD_Loc%\..
SET dst_folder=%depend_dir%\SubDyn
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars

:IceFloe
ECHO IceFloe
SET src_folder=%IceF_Loc%\..
SET dst_folder=%depend_dir%\IceFloe
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars

:IceDyn
ECHO IceDyn
SET src_folder=%IceD_Loc%\..
SET dst_folder=%depend_dir%\IceDyn
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars

:MAP
ECHO MAP
COPY "%MAP_DLL%"                     "%bin_dir%"
COPY "%MAP_DLL64%"                   "%bin_dir%"

SET src_folder=%MAP_Loc%
SET dst_folder=%depend_dir%\MAP
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

@CALL :CopyFileList
rem Change the case of this source file, if necessary:
MOVE "%dst_folder%\map.f90"   "%dst_folder%\MAP.f90"
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars

:MoorDyn
ECHO MoorDyn
SET src_folder=%MD_Loc%\..
SET dst_folder=%depend_dir%\MoorDyn
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars


:FEAMooring
ECHO FEAMooring
SET src_folder=%FEAM_Loc%\..
SET dst_folder=%depend_dir%\FEAMooring
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars



:OrcaFlex
ECHO OrcaFlex Integration
SET src_folder=%Orca_Loc%\..
SET dst_folder=%depend_dir%\OrcaFlex
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

@CALL :CopyFileList
IF /I "%1"=="%SW_ModuleOnly%" GOTO ClearVars


REM -------------------------------------

goto ClearVars

REM -------------------------------------
:CopyFileList
if exist "%dst_folder%\*" DEL "%dst_folder%\*"
for /f %%i in (%list_of_files%) DO copy /Y "%src_folder%\%%i" "%dst_folder%"

EXIT /B
REM -------------------------------------


:ClearVars
SET bin_dir=
SET depend_dir=
SET src_folder=
SET dst_folder=
SET list_of_files=
SET SW_ModuleOnly=


:Done
