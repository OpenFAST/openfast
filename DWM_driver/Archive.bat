@ECHO OFF

@SET ARCHPATH=Archive
@SET PROGNAME=DWM
@SET ARCHNAME=DWM_v%1


IF "%COMPUTERNAME%"=="APLATT-21846S" GOTO APLATT-21846S
IF "%COMPUTERNAME%"=="BJONKMAN-23080S" GOTO BJONKMAN-23080S
IF "%COMPUTERNAME%"=="GHAYMAN-17919S" GOTO GHAYMAN-17919S
IF "%COMPUTERNAME%"=="GHAYMAN-26326S" GOTO GHAYMAN-17919S


:APLATT-21846S
@SET WINZIP="C:\Program Files (x86)\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files (x86)\WinZip Self-Extractor\wzipse32.exe"
@SET SEVENZIP="C:\Program Files\7-Zip\7z.exe"
GOTO CheckSyntax

:BJONKMAN-23080S
@SET WINZIP="C:\Program Files (x86)\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files (x86)\WinZip Self-Extractor\WZIPSE22\wzipse32.exe"
@SET SEVENZIP="C:\Program Files\7-Zip\7z.exe"
GOTO CheckSyntax

:GHAYMAN-17919S
@SET WINZIP="C:\Program Files\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files (x86)\WinZip Self-Extractor\wzipse32.exe"
@SET SEVENZIP="C:\Program Files\7-Zip\7z.exe"
GOTO CheckSyntax

::=======================================================================================================

:CheckSyntax
@IF NOT "%1"==""  GOTO DeleteOld

@ECHO 
@ECHO  The syntax for creating an archive is "Archive <version>"
@ECHO.
@ECHO  Example:  "archive 1.01.00"

@GOTO Done

:DeleteOld
@IF EXIST ARCHTMP.zip DEL ARCHTMP.zip
@IF EXIST ARCHTMP.exe DEL ARCHTMP.exe
@IF EXIST %PROGRAM%.tar DEL %PROGRAM%.tar
@IF EXIST %PROGRAM%.tar.gz DEL %PROGRAM%.tar.gz


:DoIt
@ECHO.
@ECHO -------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for general Windows distribution.
@ECHO -------------------------------------------------------------------------
@ECHO.

@%WINZIP% -a -o -P ARCHTMP @ArcFiles.txt
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%"

@COPY ARCHTMP.exe %ARCHPATH%\%ARCHNAME%.exe
@DEL ARCHTMP.zip, ARCHTMP.exe

:: DWM uses Windows-specifc system calls, so we're not going to have a non-Windows archive
:: @ECHO.
:: @ECHO -------------------------------------------------------------------------
:: @ECHO Archiving %PROGNAME% for general distribution (tar.gz).
:: @ECHO -------------------------------------------------------------------------
:: @ECHO.
:: @rem first create a tar file, then compress it (gzip allows only one file)
:: @%SEVENZIP% a -ttar %PROGNAME% @ArcFiles.txt
:: @%SEVENZIP% a -tgzip %PROGNAME%.tar.gz %PROGNAME%.tar
:: @COPY %PROGNAME%.tar.gz %ARCHPATH%\%ARCHNAME%.tar.gz
:: @DEL %PROGNAME%.tar, %PROGNAME%.tar.gz



:Done
@SET ARCHPATH=
@SET ARCHNAME=
@SET PROGNAME=
@SET WINZIP=
@SET WINZIPSE=
@SET SEVENZIP=