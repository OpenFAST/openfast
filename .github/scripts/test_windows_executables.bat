@call "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" %VS_VER%

for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

echo off
setlocal enabledelayedexpansion

:: Initialize a variable to store failed solutions
set "FailedExecs="
set "OverallErrorLevel=0"

echo "Directory listing of executables"
dir build\bin
echo.

:: test OpenFAST
echo on
build\bin\OpenFAST.exe -h
echo off
if %ERRORLEVEL% NEQ 0 (
    set "FailedExecs=!FailedExecs!OpenFAST  "
    set "OverallErrorLevel=1"
    echo OpenFAST failed to run!
)
echo on


:: test TurbSim
echo on
build\bin\TurbSim.exe -h
echo off
if %ERRORLEVEL% NEQ 0 (
    set "FailedExecs=!FailedExecs!TurbSim  "
    set "OverallErrorLevel=1"
    echo TurbSim failed to run!
)
echo on


:: test FAST.Farm
echo on
build\bin\FAST.Farm.exe -h
echo off
if %ERRORLEVEL% NEQ 0 (
    set "FailedExecs=!FailedExecs!FAST.Farm  "
    set "OverallErrorLevel=1"
    echo FAST.Farm failed to run!
)
echo on



echo.
echo Test Summary:
echo off
if defined FailedExecs (
    echo The following executables failed to run:
    echo %FailedExecs%
) else (
    echo All executables ran successfully.
)
echo on

:: Set the final error level based on the overall build status
exit /b %OverallErrorLevel%
