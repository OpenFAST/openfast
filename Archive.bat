@ECHO OFF

@SET PROGNAME=AeroDyn_Driver
@SET ARCHNAME=AD_Driver_v%1
@SET ARCHPATH=Archive

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
@IF EXIST %PROGNAME%.tar DEL %PROGNAME%.tar
@IF EXIST %PROGNAME%.tar.gz DEL %PROGNAME%.tar.gz


:DoIt
@ECHO.
@ECHO -------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for general Windows distribution.
@ECHO -------------------------------------------------------------------------
@ECHO.
@%WINZIP% -a -o -P ARCHTMP @ArcFiles.txt @ArcWin.txt
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt
@COPY ARCHTMP.exe %ARCHPATH%\%ARCHNAME%.exe
@DEL ARCHTMP.zip, ARCHTMP.exe


@ECHO --------------------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for general distribution (tar.gz).
@ECHO --------------------------------------------------------------------------------------
@ECHO.
@rem first create a tar file, then compress it (gzip allows only one file)
@%SEVENZIP% a -ttar %PROGNAME% @ArcFiles.txt
@%SEVENZIP% a -tgzip %PROGNAME%.tar.gz %PROGNAME%.tar
@COPY %PROGNAME%.tar.gz %ARCHPATH%\%ARCHNAME%.tar.gz
@DEL %PROGNAME%.tar, %PROGNAME%.tar.gz


:Done
@SET ARCHNAME=
@SET ARCHPATH=
@SET PROGNAME=
@SET WINZIP=
@SET WINZIPSE=
@SET SEVENZIP=