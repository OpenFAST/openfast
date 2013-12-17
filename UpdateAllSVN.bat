@ECHO OFF

::=======================================================================================================
IF "%COMPUTERNAME%"=="APLATT-21846S" GOTO APLATT-21846S
IF "%COMPUTERNAME%"=="BJONKMAN-23080S" GOTO BJONKMAN-23080S
IF "%COMPUTERNAME%"=="MBUHL-20665S" GOTO MBUHL-20665S
IF "%COMPUTERNAME%"=="WIND-WAS13" GOTO WIND-WAS13
GOTO NotFound


:APLATT-21846S
@SET AeroDyn=
@SET ElastoDyn=
@SET FAST=
@SET HydroDyn=
@SET InflowWind=
@SET MAP=
@SET NWTC_Library=
@SET ServoDyn=
@SET SubDyn=
GOTO UpdateSVN

:BJONKMAN-23080S
CALL ./Compiling/Set_FAST_Paths.bat
@SET AeroDyn=%AD_Loc%\..\
@SET ElastoDyn=%ED_Loc%\..\
@SET FAST=%FAST_Loc%\..\
@SET HydroDyn=%HD_Loc%\..\
@SET InflowWind=%IfW_Loc%\..\
@SET MAP=%MAP_Loc%\..\..\
@SET NWTC_Library=%NWTC_Lib_Loc%\..\
@SET ServoDyn=%SrvD_Loc%\..\
@SET SubDyn=%SD_Loc%\..\
GOTO UpdateSVN

:MBUHL-20665S
@SET AeroDyn="M:\CAEtools\Simulators\AeroDyn\branches\Framework"
@SET ElastoDyn="M:\CAEtools\Simulators\FAST\branches\BJonkman"
@SET FAST="M:\CAEtools\Simulators\FAST\branches\BJonkman"
@SET HydroDyn="M:\CAEtools\Simulators\HydroDyn\branches\HydroDyn_Modularization"
@SET InflowWind="M:\CAEtools\Simulators\InflowWind\branches\modularization"
@SET MAP="M:\CAEtools\Simulators\MAP\trunk"
@SET NWTC_Library="M:\CAEtools\Miscellaneous\NWTC_Library\trunk"
@SET ServoDyn="M:\CAEtools\Simulators\FAST\branches\BJonkman"
@SET SubDyn="M:\CAEtools\Simulators\SubDyn\branches\v0.4"
GOTO UpdateSVN

:WIND-WAS13
IF "%USERNAME%"=="mbuhlx" GOTO Mbuhl_windwas13
GOTO NotFound

:Mbuhl_windwas13
@SET AeroDyn="M:\CAEtools\Simulators\AeroDyn\branches\Framework"
@SET ElastoDyn="M:\CAEtools\Simulators\FAST\branches\BJonkman"
@SET FAST="M:\CAEtools\Simulators\FAST\branches\BJonkman"
@SET HydroDyn="M:\CAEtools\Simulators\HydroDyn\branches\HydroDyn_Modularization"
@SET InflowWind="M:\CAEtools\Simulators\InflowWind\branches\modularization"
@SET MAP="M:\CAEtools\Simulators\MAP\trunk"
@SET NWTC_Library="M:\CAEtools\Miscellaneous\NWTC_Library\trunk"
@SET ServoDyn="M:\CAEtools\Simulators\FAST\branches\BJonkman"
@SET SubDyn="M:\CAEtools\Simulators\SubDyn\branches\v0.4"
GOTO UpdateSVN

:: Update all the SVN repositories used for this project.

:UpdateSVN
@ECHO.
@ECHO Using folder definitions for %HOSTNAME% to update all modules in this project.
@ECHO.

SVN update %AeroDyn%
SVN update %ElastoDyn%
SVN update %FAST%
SVN update %HydroDyn%
SVN update %InflowWind%
SVN update %MAP%
SVN update %NWTC_Library%
SVN update %ServoDyn%
SVN update %SubDyn%
GOTO End

:: Delete temporary environment variables.

:NotFOund
@ECHO.
@ECHO No path settings for %USERNAME% on %HOSTNAME%.  Edit this script.
@ECHO.

:End
@SET AeroDyn=
@SET ElastoDyn=
@SET FAST=
@SET HydroDyn=
@SET InflowWind=
@SET MAP=
@SET NWTC_Library=
@SET ServoDyn=
@SET SubDyn=

PAUSE
