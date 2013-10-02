@ECHO OFF

::=======================================================================================================
IF "%COMPUTERNAME%"=="APLATT-21846S" GOTO APLATT-21846S
IF "%COMPUTERNAME%"=="BJONKMAN-23080S" GOTO BJONKMAN-23080S
IF "%COMPUTERNAME%"=="MBUHL-20665S" GOTO MBUHL-20665S
IF "%COMPUTERNAME%"=="wind-was13" GOTO MBUHL-20665S

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

:: Delete temporary environment variables.

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
