@call "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" %VS_VER%

for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

@REM Make the script that generates the git version description ignore dirty
@REM since building the Visual Studio projects modifies files
powershell -command "(Get-Content -Path '.\vs-build\CreateGitVersion.bat') -replace '--dirty', '' | Set-Content -Path '.\vs-build\CreateGitVersion.bat'"

echo on

@REM Build all solutions
devenv vs-build/AeroDisk/AeroDisk_Driver.sln /Build "Release|x64"
devenv vs-build/AeroDyn/AeroDyn_Driver.sln /Build "Release|x64"
devenv vs-build/AeroDyn/AeroDyn_Driver.sln /Build "Release_OpenMP|x64"
devenv vs-build/AeroDyn_Inflow_c_binding/AeroDyn_Inflow_c_binding.sln /Build "Release|x64"
devenv vs-build/AeroDyn_Inflow_c_binding/AeroDyn_Inflow_c_binding.sln /Build "Release_OpenMP|x64"
devenv vs-build/BeamDyn/BeamDyn-w-registry.sln /Build "Release|x64"
devenv vs-build/Discon/Discon.sln /Build "Release|x64" 
devenv vs-build/FAST-farm/FAST-Farm.sln /Build "Release|x64"
devenv vs-build/FAST-farm/FAST-Farm.sln /Build "Release_OpenMP|x64"
devenv vs-build/HydroDyn/HydroDynDriver.sln /Build "Release|x64"
devenv vs-build/HydroDyn_c_binding/HydroDyn_c_binding.sln /Build "Release|x64"
devenv vs-build/InflowWind_c_binding/InflowWind_c_binding.sln /Build "Release|x64"
devenv vs-build/InflowWind/InflowWind_driver.sln /Build "Release|x64"
devenv vs-build/InflowWind/InflowWind_driver.sln /Build "Release_OpenMP|x64"
devenv vs-build/MoorDyn/MoorDynDriver.sln /Build "Release|x64"
devenv vs-build/MoorDyn_c_binding/MoorDyn_c_binding.sln /Build "Release|x64"
devenv vs-build/FAST/FAST.sln /Build "Release|x64"
devenv vs-build/SeaState/SeaStateDriver.sln /Build "Release|x64"
devenv vs-build/SeaState_c_binding/SeaState_c_binding.sln /Build "Release|x64"
devenv vs-build/SimpleElastoDyn/SimpleElastoDyn_Driver.sln /Build "Release|x64"
devenv vs-build/SubDyn/SubDyn.sln /Build "Release|x64"
devenv vs-build/TurbSim/TurbSim.vfproj /Build "Release|x64"
devenv vs-build/UnsteadyAero/UnsteadyAero.sln /Build "Release|x64"

@REM Build MATLAB solution last
devenv vs-build/FAST/FAST.sln /Build "Release_Matlab|x64"

@REM Copy controllers to bin directory
xcopy .\reg_tests\r-test\glue-codes\openfast\5MW_Baseline\ServoData\*.dll .\build\bin\ /y

exit /b %ERRORLEVEL%
