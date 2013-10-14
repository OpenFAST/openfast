@ECHO OFF

@SET ARCHROOT=AD
@SET PROGNAME=AeroDyn

::=======================================================================================================
IF "%COMPUTERNAME%"=="APLATT-21846S" GOTO APLATT-21846S
IF "%COMPUTERNAME%"=="BJONKMAN-23080S" GOTO BJONKMAN-23080S
IF "%COMPUTERNAME%"=="MBUHL-20665S" GOTO MBUHL-20665S

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

:MBUHL-20665S
@SET WINZIP="C:\Program Files (x86)\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files (x86)\WinZip Self-Extractor\wzipse32.exe"
@SET SEVENZIP="C:\Program Files (x86)\7-Zip\7z.exe"
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
@IF EXIST ARCHTMP.tar DEL ARCHTMP.tar
@IF EXIST ARCHTMP.tar.gz DEL ARCHTMP.tar.gz


:DoIt
@ECHO.
@ECHO -------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for general Windows distribution.
@ECHO -------------------------------------------------------------------------
@ECHO.
@%WINZIP% -a -o -P ARCHTMP @ArcFiles.txt @SourceFiles.txt
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt
@COPY ARCHTMP.exe Archive\%ARCHROOT%_v%1.exe
@DEL ARCHTMP.zip, ARCHTMP.exe

:: @ECHO.
:: @ECHO -------------------------------------------------------------------------
:: @ECHO Archiving %PROGNAME% for maintenance.
:: @ECHO -------------------------------------------------------------------------
:: @ECHO.
:: @%WINZIP% -a -o -P ARCHTMP @ArcFiles.txt @ArcMaint.txt @SourceFiles.txt
:: @%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt
:: @COPY ARCHTMP.exe Archive\%ARCHROOT%_v%1_Maint.exe
:: @DEL ARCHTMP.zip, ARCHTMP.exe

@ECHO.
@ECHO -------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for general distribution (tar.gz).
@ECHO -------------------------------------------------------------------------
@ECHO.
@rem first create a tar file, then compress it (gzip allows only one file)
@%SEVENZIP% a -ttar ARCHTMP @ArcFiles.txt @SourceFiles.txt
@%SEVENZIP% a -tgzip ARCHTMP.tar.gz ARCHTMP.tar
@COPY ARCHTMP.tar.gz Archive\%ARCHROOT%_v%1.tar.gz
@DEL ARCHTMP.tar, ARCHTMP.tar.gz



:Done
@SET ARCHROOT=
@SET PROGNAME=
@SET WINZIP=
@SET WINZIPSE=
@SET SEVENZIP=