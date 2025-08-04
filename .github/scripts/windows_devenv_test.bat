@REM Check if devenv actually works 
for /f "usebackq tokens=1* delims=: " %%i in (`"C:\Program Files (x86)\Microsoft Visual Studio\Installer\vswhere.exe" -latest`) do (
  if /i "%%i"=="productPath" set devenv=%%j
)

@REM above command finds devenv.exe, but that opens the VS instance.  We need the devenv.com version
set devenv=%devenv:devenv.exe=devenv.com%

echo Using Visual Studio: %devenv%

"%devenv%" /?

exit /b %ERRORLEVEL%
