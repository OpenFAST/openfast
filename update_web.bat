@ECHO OFF

set LOC=Y:\Wind\WindWeb\nwtc\docs

:CopyFiles
copy ChangeLog.txt             %LOC%\AeroDyn-ChangeLog.txt

:Done
set LOC=
