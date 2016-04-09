@ECHO OFF

@SET ARCHPATH=Archive
@SET PROGNAME=OrcaFlexInterface_Samples
@SET ARCHNAME=OrcaFlexInterface_Samples-FAST_v%1

::=======================================================================================================
IF "%COMPUTERNAME%"=="APLATT-21846S" GOTO APLATT-21846S
IF "%COMPUTERNAME%"=="BJONKMAN-23080S" GOTO BJONKMAN-23080S

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

::=======================================================================================================

:CheckSyntax
IF NOT "%1"==""  GOTO DeleteOld

@ECHO 
@ECHO  The syntax for creating an archive is "Archive <version>"
@ECHO.
@ECHO  Example:  "archive 7.01.00"

@GOTO Done

:DeleteOld
@IF EXIST ARCHTMP.zip DEL ARCHTMP.zip
@IF EXIST ARCHTMP.exe DEL ARCHTMP.exe

:DoIt
@ECHO.
@ECHO --------------------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for general distribution on Windows.
@ECHO --------------------------------------------------------------------------------------
@ECHO.

@%WINZIP% -a -o -P ARCHTMP @ArcFiles.txt
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt

@COPY ARCHTMP.exe %ARCHPATH%\%ARCHNAME%.exe
@COPY ARCHTMP.zip %ARCHPATH%\%ARCHNAME%.zip
@DEL ARCHTMP.zip, ARCHTMP.exe


:Done
@SET ARCHNAME=
@SET ARCHPATH=
@SET PROGNAME=
@SET WINZIP=
@SET WINZIPSE=
@SET SEVENZIP=