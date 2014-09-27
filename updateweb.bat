@ECHO OFF

set LOC=Y:\Wind\WindWeb\nwtc\docs

:CopyFiles
copy ChangeLog.txt                    %LOC%\FAST_ChangeLog.txt
copy README_FAST8.pdf                 %LOC%\README_FAST8.pdf

:Done
set LOC=
