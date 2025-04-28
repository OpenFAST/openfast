@call "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" %VS_VER%

for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

@REM Make the script that generates the git version description ignore dirty
@REM since building the Visual Studio projects modifies files
powershell -command "(Get-Content -Path '.\vs-build\CreateGitVersion.bat') -replace '--dirty', '' | Set-Content -Path '.\vs-build\CreateGitVersion.bat'"

echo on

@REM Build all solutions
devenv vs-build/OpenFAST.sln /Build "Release|x64"
devenv vs-build/OpenFAST.sln /Build "Release_OpenMP|x64"

@REM Build MATLAB solution last
devenv vs-build/OpenFAST.sln /Build "Release_Matlab|x64"

@REM Copy controllers to bin directory
@REM xcopy .\reg_tests\r-test\glue-codes\openfast\5MW_Baseline\ServoData\*.dll .\build\bin\ /y

exit /b %ERRORLEVEL%
