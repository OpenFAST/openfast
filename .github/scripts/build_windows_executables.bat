@call "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" %VS_VER%

for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

@REM Make the script that generates the git version description ignore dirty
@REM since building the Visual Studio projects modifies files
powershell -command "(Get-Content -Path '.\vs-build\CreateGitVersion.bat') -replace '--dirty', '' | Set-Content -Path '.\vs-build\CreateGitVersion.bat'"

setlocal enabledelayedexpansion

:: Initialize a variable to store failed solutions
set "FailedSolutions="
set "OverallErrorLevel=0"


echo "Build all projects (Release|64)"
devenv vs-build/OpenFAST.sln /Build "Release|x64"
if %ERRORLEVEL% NEQ 0 (
    set "FailedSolutions=!FailedSolutions!Release  "
    set "OverallErrorLevel=1"
    echo Build of OpenFAST.sln Release failed!
)


echo "Build all OpenMP projects (OpenMP_Release|64)"
devenv vs-build/OpenFAST.sln /Build "OpenMP_Release|x64"
if %ERRORLEVEL% NEQ 0 (
    set "FailedSolutions=!FailedSolutions!OpenMP_Release  "
    set "OverallErrorLevel=1"
    echo Build of OpenFAST.sln OpenMP_Release failed!
)


echo "Build OpenFAST-Simulink shared library (Matlab_Release|x64)"
devenv vs-build/OpenFAST.sln /Build "Matlab_Release|x64"
if %ERRORLEVEL% NEQ 0 (
    set "FailedSolutions=!FailedSolutions!Matlab_Release  "
    set "OverallErrorLevel=1"
    echo Build of OpenFAST.sln Matlab_Release failed!
)


echo "Build Summary:"
if defined FailedSolutions (
    echo The following solutions failed to build:
    echo %FailedSolutions%
) else (
    echo All solutions built successfully.
)

@echo off
setlocal enabledelayedexpansion

cd /d build\bin || exit /b 1

for %%F in (*_Release*) do (
    set "name=%%~nxF"
    set "newname=!name:_Release=!"
    if not "!name!"=="!newname!" (
        ren "%%F" "!newname!"
    )
)
for %%F in (*_Matlab*) do (
    set "name=%%~nxF"
    set "newname=!name:_Matlab=!"
    if not "!name!"=="!newname!" (
        ren "%%F" "!newname!"
    )
)

endlocal

echo "List executables in build\bin"
dir build\bin

:: Set the final error level based on the overall build status
exit /b %OverallErrorLevel%
