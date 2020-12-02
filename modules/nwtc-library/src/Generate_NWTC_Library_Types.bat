@ECHO OFF

SET REG_Loc=..\..\..\build\bin
::SET Registry=%REG_Loc%\Registry_Win32.exe
SET Registry=%REG_Loc%\Registry.exe

SET ModuleName=%1

GOTO %ModuleName%

REM ----------------------------------------------------------------------------
REM ---------------- RUN THE REGISTRY TO AUTO-GENERATE FILES -------------------
REM ----------------------------------------------------------------------------
ECHO on
:mesh
%REGISTRY% Registry_NWTC_Library_typedef_mesh.txt  -noextrap
goto end

:nomesh
%REGISTRY% Registry_NWTC_Library_typedef_nomesh.txt  -noextrap


:end
REM ----------------------------------------------------------------------------
REM ------------------------- CLEAR MEMORY -------------------------------------
REM ----------------------------------------------------------------------------
ECHO. 


SET REGISTRY=
SET REG_Loc=

SET ModuleName=
:Done
