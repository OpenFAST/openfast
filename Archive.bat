@ECHO OFF

@SET ARCHPATH=Archive
@SET PROGNAME=FAST
@SET ARCHNAME=FAST_v%1

@SET WINZIP="C:\Program Files\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files\WinZip Self-Extractor\wzipse32.exe"


IF NOT "%1"==""  GOTO DeleteOld

@ECHO 
@ECHO  The syntax for creating an archive is "Archive <version>"
@ECHO.
@ECHO  Example:  "archive 7.01.00"

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

@COPY ARCHTMP.exe %ARCHPATH%\%ARCHNAME%.exe
@DEL ARCHTMP.zip, ARCHTMP.exe

@ECHO.
@ECHO ------------------------------------------
@ECHO Archiving %PROGNAME% for internal use by including certification tests.
@ECHO ------------------------------------------
@ECHO.

@%WINZIP% -a -o -P ARCHTMP @ArcAll.txt
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt

@COPY ARCHTMP.exe %ARCHPATH%\%ARCHNAME%_all.exe
@DEL ARCHTMP.zip, ARCHTMP.exe


:Done
@SET ARCHNAME=
@SET ARCHPATH=
@SET PROGNAME=
@SET WINZIP=
@SET WINZIPSE=