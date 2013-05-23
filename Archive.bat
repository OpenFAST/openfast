@ECHO OFF

@SET ARCHROOT=NWTC_Lib
@SET PROGNAME=NWTC_Library

::=======================================================================================================
:: PLEASE NOTE
:: -----------
:: If you need to change the file locations to work for your system, please copy an existing set and put
:: your modified version at the end of this section.  You do not need to comment out other versions.
:: Don't forget to name your set.

:: Someone's setup:
@SET WINZIP="C:\Program Files (x86)\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files (x86)\WinZip Self-Extractor\WZIPSE22\wzipse32.exe"
@SET TARZIP="C:\Program Files\7-Zip\7z.exe"

:: Marshall's setup:
@SET WINZIP="C:\Program Files (x86)\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files (x86)\WinZip Self-Extractor\wzipse32.exe"
@SET TARZIP="C:\Program Files\7-Zip\7z.exe"
::=======================================================================================================


@IF NOT "%1"==""  GOTO DeleteOld

@ECHO 
@ECHO  The syntax for creating an archive is "Archive <version>"
@ECHO.
@ECHO  Example:  "archive 1.01.00"

@GOTO Done


:DeleteOld
@IF EXIST ARCHTMP.zip DEL ARCHTMP.zip
@IF EXIST %ARCHNAME%.exe DEL %ARCHNAME%.exe
@IF EXIST ARCHTMP.tar DEL ARCHTMP.tar
@IF EXIST %ARCHNAME%.tar.gz DEL %ARCHNAME%.tar.gz


:DoIt
@ECHO.
@ECHO -------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for general Windows distribution.
@ECHO -------------------------------------------------------------------------
@ECHO.
@%WINZIP% -a -o -P ARCHTMP @ArcFiles.txt
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt
@COPY ARCHTMP.exe Archive\%ARCHROOT%_v%1.exe
@DEL ARCHTMP.zip, ARCHTMP.exe

@ECHO.
@ECHO -------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for maintenance.
@ECHO -------------------------------------------------------------------------
@ECHO.
@%WINZIP% -a -o -P ARCHTMP @ArcFiles.txt @ArcMaint.txt
@%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt
@COPY ARCHTMP.exe Archive\%ARCHROOT%_v%1_Maint.exe
@DEL ARCHTMP.zip, ARCHTMP.exe

@ECHO.
@ECHO -------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for general distribution (tar.gz).
@ECHO -------------------------------------------------------------------------
@ECHO.
@rem first create a tar file, then compress it (gzip allows only one file)
@%TARZIP% a -ttar ARCHTMP @ArcFiles.txt
@%TARZIP% a -tgzip ARCHTMP.tar.gz ARCHTMP.tar
@COPY ARCHTMP.tar.gz Archive\%ARCHROOT%_v%1.tar.gz
@DEL ARCHTMP.tar, ARCHTMP.tar.gz



:Done
@SET ARCHROOT=
@SET PROGNAME=
@SET WINZIP=
@SET WINZIPSE=
@SET TARZIP=