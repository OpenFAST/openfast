@ECHO.
@ECHO -------------------------------------------------------------------------
@ECHO Updating bin/Registry_win32.exe.
@ECHO -------------------------------------------------------------------------
@ECHO.

cd bin
svn export --force --depth files "https://windsvn.nrel.gov/FAST/branches/FAST_Registry" "./temp"
cd temp
move registry.exe ..
del /q *.*
cd ..
rmdir temp
@IF EXIST Registry_win32.exe DEL Registry_win32.exe
rename registry.exe Registry_win32.exe

@ECHO.
@ECHO -------------------------------------------------------------------------
@ECHO Updating source dependencies.
@ECHO -------------------------------------------------------------------------
@ECHO.

cd ..
cd Source\dependencies

@ECHO -------------------------------------------------------------------------
@ECHO Updating NetLib source files.
@ECHO -------------------------------------------------------------------------
@ECHO.

svn export --force "https://windsvn.nrel.gov/FAST/branches/BJonkman/Source/dependencies/NetLib" "./NetLib"

@ECHO -------------------------------------------------------------------------
@ECHO Updating NWTC_Library source files.
@ECHO -------------------------------------------------------------------------
@ECHO.

svn export --force --depth files "https://windsvn.nrel.gov/NWTC_Library/trunk/source" "./NWTC_Library"
