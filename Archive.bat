@ECHO OFF

IF NOT "%1"==""  GOTO SetVars

ECHO 
ECHO  The syntax for creating an archive is "Archive <version>"
ECHO.
ECHO  Example:  "archive 602 [all]"

GOTO Done

:SetVars
SET ARCHPATH=Archive
SET PROGNAME=HydroDyn
SET WINZIP="C:\Program Files\WinZip\WZZip"
SET WINZIPSE="C:\Program Files\WinZip Self-Extractor\wzipse32.exe"


IF "%2"=="all"  GOTO ArchAll

ECHO Archiving %PROGNAME% for distribution.
SET ARCHNAME=%PROGNAME%_v%1
SET FILELIST=ArcFiles.txt

GOTO DeleteOld

:ArchAll
ECHO Archiving %PROGNAME% for internal use by including certification tests.
SET ARCHNAME=%PROGNAME%_v%1_all
SET FILELIST=ArcAll.txt


:DeleteOld
IF EXIST ARCHTMP.zip DEL ARCHTMP.zip
IF EXIST %ARCHNAME%.exe DEL %ARCHNAME%.exe


:DoIt
%WINZIP% -a -o -P ARCHTMP @%FILELIST%
ECHO Creating self-extracting archive...
%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt 

ECHO The self-extracting archive has been created.
COPY ARCHTMP.exe %ARCHPATH%\%ARCHNAME%.exe
DEL ARCHTMP.zip, ARCHTMP.exe


:Done
SET ARCHNAME=
SET ARCHPATH=
SET PROGNAME=
SET FILELIST=
SET WINZIP=
SET WINZIPSE=