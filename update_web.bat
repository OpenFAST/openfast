@ECHO OFF

@IF NOT "%1"==""  GOTO Copy

@ECHO 
@ECHO  The syntax for updating files on the HydroDyn page is "update_web <version>"
@ECHO.
@ECHO  Example:  "update_web 410"

@GOTO Done

:Copy

set LOC=y:\wind\WindWeb\designcodes\simulators\hydrodyn

copy CreatePage.pl                 %LOC%\CreatePage.pl
copy ChangeLog.txt                 %LOC%\ChangeLog.txt
copy Archive\HydroDyn_v%1.exe      %LOC%\HydroDyn_v%1.exe
copy HydroDyn.pdf                  %LOC%\HydroDyn.pdf

perl %LOC%\CreatePage.pl

set LOC=

:Done
