@ECHO OFF

set LOC=y:\wind\WindWeb\nwtc\docs


:CopyFiles
copy ChangeLog.txt                  %LOC%\InflowWind_ChangeLog.txt
copy Docs\InflowWind_Manual.pdf     %LOC%\InflowWind_Manual.pdf

:Done
set LOC=


