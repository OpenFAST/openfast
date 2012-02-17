@ECHO OFF

@IF NOT "%1"==""  GOTO Copy

@ECHO 
@ECHO  The syntax for updating files on the AeroDyn page is "update_web <version>"
@ECHO.
@ECHO  Example:  "update_web 410"

@GOTO Done

:Copy

rem set LOC=\\netapp01\websites\wind\designcodes\simulators\aerodyn
set LOC=Y:\Wind\WindWeb\designcodes\simulators\aerodyn

copy CreatePage.pl             %LOC%\CreatePage.pl
copy ChangeLog.txt             %LOC%\ChangeLog.txt
copy Archive\AD_v%1.exe        %LOC%\AD_v%1.exe

perl %LOC%\CreatePage.pl

set LOC=

:Done
