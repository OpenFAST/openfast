@ECHO OFF

@IF NOT "%1"==""  GOTO Copy

@ECHO 
@ECHO  The syntax for updating files on the FAST page is "update_web <version>"
@ECHO.
@ECHO  Example:  "update_web 7.01.00"

@GOTO Done

:Copy

set LOC=Y:\Wind\WindWeb\designcodes\simulators\fast\alpha

copy CreatePage.pl             %LOC%\CreatePage.pl
copy ChangeLog.txt             %LOC%\ChangeLog.txt
copy Archive\FAST_v%1.exe      %LOC%\FAST_v%1.exe
copy Archive\FAST_v%1.tar.gz   %LOC%\FAST_v%1.tar.gz
copy Archive\FAST_v%1_all.exe  %LOC%\FAST_v%1_all.exe
copy README_FAST8.pdf          %LOC%\README_FAST8.pdf

perl %LOC%\CreatePage.pl

set LOC=

:Done
