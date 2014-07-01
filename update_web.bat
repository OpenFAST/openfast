@ECHO OFF

set LOC=Y:\Wind\WindWeb\designcodes\simulators\aerodyn\alpha


@IF /I "%1"=="-NOFILES" GOTO UpdateWeb
@IF NOT "%1"=="" GOTO CopyFiles


@ECHO 
@ECHO  The syntax for updating files on the AeroDyn page is "update_web <version>"
@ECHO.
@ECHO  Example:  "update_web 13.01.00"
@ECHO.
@ECHO  The syntax for updating the web page without copying the other files is "update_web -nofiles"
@ECHO.

@GOTO Done


:CopyFiles
copy ChangeLog.txt             %LOC%\ChangeLog.txt
copy Archive\AD_v%1.exe        %LOC%\AD_v%1.exe
copy Archive\AD_v%1.tar.gz     %LOC%\AD_v%1.tar.gz

rem copy UserGuideAddendum_AeroDyn13Interface.pdf %LOC%\UserGuideAddendum_AeroDyn13Interface.pdf

:UpdateWeb
copy CreatePage.pl             %LOC%\CreatePage.pl

perl %LOC%\CreatePage.pl


:Done
set LOC=
