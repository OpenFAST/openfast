@ECHO OFF

set LOC=Y:\Wind\WindWeb\nwtc\docs

:CopyFiles
copy ChangeLog.txt             %LOC%\AeroDyn14_ChangeLog.txt

:Done
set LOC=
