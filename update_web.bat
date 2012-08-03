@ECHO OFF

@IF NOT "%1"==""  GOTO Copy

@ECHO 
@ECHO  The syntax for updating files on the InflowWind page is "update_web <version>"
@ECHO.
@ECHO  Example:  "update_web 1.03.01"

@GOTO Done

:Copy

set LOC=y:\wind\WindWeb\designcodes\miscellaneous\inflowwind

copy CreatePage.pl                 %LOC%\CreatePage.pl
copy ChangeLog.txt                 %LOC%\ChangeLog.txt
copy Archive\InflowWind_v%1.exe    %LOC%\InflowWind_v%1.exe
copy InflowWind.pdf                %LOC%\InflowWind.pdf

perl %LOC%\CreatePage.pl

set LOC=

:Done
