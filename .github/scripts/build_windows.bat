@call "C:\Program Files (x86)\Intel\oneAPI\setvars-vcvarsall.bat" %VS_VER%

for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

devenv vs-build/FAST/FAST.sln /Build "Release|x64"
devenv vs-build/AeroDyn/AeroDyn_Driver.sln /Build "Release|x64"
devenv vs-build/AeroDyn/AeroDyn_Driver.sln /Build "Release_OpenMP|x64"
devenv vs-build/AeroDyn_Inflow_c_binding/AeroDyn_Inflow_c_binding.sln /Build "Release|x64"
devenv vs-build/AeroDyn_Inflow_c_binding/AeroDyn_Inflow_c_binding.sln /Build "Release_OpenMP|x64"
devenv vs-build/BeamDyn/BeamDyn-w-registry.sln /Build "Release|x64"
devenv vs-build/Discon/Discon.sln /Build "Release|x64" 
devenv vs-build/TurbSim/TurbSim.vfproj /Build "Release|x64"

exit /b %ERRORLEVEL%
