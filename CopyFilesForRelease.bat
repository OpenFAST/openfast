@ECHO OFF


::=======================================================================================================
:: run a script to set all of the path locations needed to copy files

CALL .\Compiling\Set_FAST_Paths.bat

::=======================================================================================================

set bin_dir=.\bin
set depend_dir=.\Source\dependencies


:BinDir
COPY "%REG_Loc%\registry.exe"        "%bin_dir%\Registry_win32.exe"
COPY "%MAP_DLL%"                     "%bin_dir%"
COPY "%CRUNCH%"                      "%bin_dir%"


:dependencies

:NWTC_Library
ECHO NWTC_Library
SET src_folder=%NWTC_Lib_Loc%\..
SET dst_folder=%depend_dir%\NWTC_Library
SET list_of_files=%src_folder%\SourceFiles.txt

if exist "%dst_folder%\*" DEL "%dst_folder%\*"
for /f %%i in (%list_of_files%) DO copy /Y "%src_folder%\%%i" "%dst_folder%"

:NetLib
ECHO NetLib
SET src_folder=%NETLIB_Loc%\..
SET dst_folder=%depend_dir%\NetLib
SET list_of_files=%src_folder%\SourceFiles.txt

if exist "%dst_folder%\*" DEL "%dst_folder%\*"
for /f %%i in (%list_of_files%) DO copy /Y "%src_folder%\%%i" "%dst_folder%"


:ElastoDyn
ECHO ElastoDyn
SET src_folder=%ED_Loc%\..
SET dst_folder=%depend_dir%\ElastoDyn
SET list_of_files=%src_folder%\FAST_SourceFiles_ED.txt

if exist "%dst_folder%\*" DEL "%dst_folder%\*"
for /f %%i in (%list_of_files%) DO copy /Y "%src_folder%\%%i" "%dst_folder%"

:ServoDyn
ECHO ServoDyn
SET src_folder=%SrvD_Loc%\..
SET dst_folder=%depend_dir%\ServoDyn
SET list_of_files=%src_folder%/FAST_SourceFiles_SrvD.txt

if exist "%dst_folder%\*" DEL "%dst_folder%\*"
for /f %%i in (%list_of_files%) DO copy /Y "%src_folder%\%%i" "%dst_folder%"

:InflowWind
ECHO InflowWind
SET src_folder=%IfW_Loc%\..
SET dst_folder=%depend_dir%\InflowWind
SET list_of_files=%src_folder%/FAST_SourceFiles.txt

if exist "%dst_folder%\*" DEL "%dst_folder%\*"
for /f %%i in (%list_of_files%) DO copy /Y "%src_folder%\%%i" "%dst_folder%"

:AeroDyn
ECHO AeroDyn
SET src_folder=%AD_Loc%\..
SET dst_folder=%depend_dir%\AeroDyn
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

if exist "%dst_folder%\*" DEL "%dst_folder%\*"
for /f %%i in (%list_of_files%) DO copy  "%src_folder%\%%i" "%dst_folder%"

:HydroDyn
ECHO HydroDyn
SET src_folder=%HD_Loc%\..
SET dst_folder=%depend_dir%\HydroDyn
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

if exist "%dst_folder%\*" DEL "%dst_folder%\*"
for /f %%i in (%list_of_files%) DO copy /Y "%src_folder%\%%i" "%dst_folder%"

:SubDyn
ECHO SubDyn
SET src_folder=%SD_Loc%\..
SET dst_folder=%depend_dir%\SubDyn
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

if exist "%dst_folder%\*" DEL "%dst_folder%\*"
for /f %%i in (%list_of_files%) DO copy /Y "%src_folder%\%%i" "%dst_folder%"

:IceFloe
ECHO IceFloe
SET src_folder=%IceF_Loc%\..
SET dst_folder=%depend_dir%\IceFloe
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

if exist "%dst_folder%\*" DEL "%dst_folder%\*"
for /f %%i in (%list_of_files%) DO copy /Y "%src_folder%\%%i" "%dst_folder%"

:MAP
ECHO MAP
SET src_folder=%MAP_Loc%\..
SET dst_folder=%depend_dir%\MAP
SET list_of_files=%src_folder%\FAST_SourceFiles.txt

if exist "%dst_folder%\*" DEL "%dst_folder%\*"
for /f %%i in (%list_of_files%) DO copy /Y "%src_folder%\%%i" "%dst_folder%"

:FEAMooring
ECHO Skipping FEAMooring


:ClearVars
SET bin_dir=
SET depend_dir=
SET src_folder=
SET dst_folder=
SET list_of_files=

:Done
