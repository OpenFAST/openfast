@ECHO OFF

set LOC=Y:\Wind\WindWeb\nwtc\docs


:CopyFiles
copy ChangeLog.txt                            %LOC%\NWTC_Library_ChangeLog.txt
copy Documentation\NWTC_Library_Routines.txt  %LOC%\NWTC_Library_Routines.txt


:Done
set LOC=

