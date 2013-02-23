@ECHO OFF

@IF NOT "%1"==""  GOTO SetVars

@ECHO 
@ECHO  The syntax for creating an archive is "Archive <version>"
@ECHO.
@ECHO  Example:  "archive 13.01.00"

@GOTO Done

:SetVars
@SET ARCHNAME=AD_v%1
@SET PROGNAME=AeroDyn
@SET FILELIST=ArcFiles.txt
@SET ARCHPATH=Archive

@SET WINZIP="C:\Program Files (x86)\WinZip\WZZip"
@SET WINZIPSE="C:\Program Files (x86)\WinZip Self-Extractor\WZIPSE22\wzipse32.exe"
@SET SEVENZIP="C:\Program Files\7-Zip\7z.exe"


@IF EXIST ARCHTMP.zip DEL ARCHTMP.zip
@IF EXIST ARCHTMP.exe DEL ARCHTMP.exe
@IF EXIST ARCHTMP.tar DEL ARCHTMP.tar
@IF EXIST ARCHTMP.tar.gz DEL ARCHTMP.tar.gz


:DoIt

@ECHO.
@ECHO -------------------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for Windows distribution.
@ECHO --------------------------------------------------------------------------------------
@ECHO.

%WINZIP% -a -o -P ARCHTMP @%FILELIST%
%WINZIPSE% ARCHTMP.zip -d. -y -win32 -le -overwrite -st"Unzipping %PROGNAME%" -m Disclaimer.txt

@COPY ARCHTMP.exe %ARCHPATH%\%ARCHNAME%.exe
@DEL ARCHTMP.zip, ARCHTMP.exe


@ECHO.
@ECHO -------------------------------------------------------------------------------------
@ECHO Archiving %PROGNAME% for general (Linux) distribution. (tar.gz)
@ECHO --------------------------------------------------------------------------------------
@ECHO.

rem NOTE that adding ./ in front of the file names in ArcFiles.txt will remove the relative paths in the archive
@rem first create a tar file, then compress it (gzip allows only one file)
%SEVENZIP% a -ttar ARCHTMP @%FILELIST%
%SEVENZIP% a -tgzip ARCHTMP.tar.gz ARCHTMP.tar

@COPY ARCHTMP.tar.gz %ARCHPATH%\%ARCHNAME%.tar.gz
@DEL ARCHTMP.tar, ARCHTMP.tar.gz

:Done
@SET ARCHNAME=
@SET ARCHPATH=
@SET PROGNAME=
@SET FILELIST=
@SET WINZIP=
@SET WINZIPSE=
@SET SEVENZIP=




