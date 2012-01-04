@ECHO OFF

@SET ARCHROOT=NWTC_Lib
@SET PROGNAME=NWTC_Library

@SET WINZIP="C:\Program Files\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files\WinZip Self-Extractor\wzipse32.exe"

@IF NOT "%1"==""  GOTO DeleteOld

@ECHO 
@ECHO  The syntax for creating an archive is "Archive <version>"
@ECHO.
@ECHO  Example:  "archive 1.01.00"

@GOTO Done


:DeleteOld
@IF EXIST ARCHTMP.zip DEL ARCHTMP.zip
@IF EXIST %ARCHNAME%.exe DEL %ARCHNAME%.exe


:DoIt
@ECHO.
@ECHO ------------------------------------------
@ECHO Archiving %PROGNAME% for general distribution.
@ECHO ------------------------------------------
@ECHO.
@%WINZIP% -a -o -P ARCHTMP @ArcFiles.txt
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt
@COPY ARCHTMP.exe Archive\%ARCHROOT%_v%1.exe
@DEL ARCHTMP.zip, ARCHTMP.exe

@ECHO.
@ECHO ---------------------------------
@ECHO Archiving %PROGNAME% for maintenance.
@ECHO ---------------------------------
@ECHO.
@%WINZIP% -a -o -P ARCHTMP @ArcFiles.txt @ArcMaint.txt
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt
@COPY ARCHTMP.exe Archive\%ARCHROOT%_v%1_Maint.exe
@DEL ARCHTMP.zip, ARCHTMP.exe


:Done
@SET ARCHROOT=
@SET PROGNAME=
@SET WINZIP=
@SET WINZIPSE=