@ECHO OFF

@IF NOT "%1"==""  GOTO Copy

@ECHO 
@ECHO  The syntax for updating files on the WT_Perf page is "updateweb <version>"
@ECHO.
@ECHO  Example:  "updateweb 1.01.08"

@GOTO Done

:Copy

set LOC=Y:\Wind\WindWeb\designcodes\miscellaneous\nwtc_subs

copy ChangeLog.txt                        %LOC%\ChangeLog.txt
copy CreatePage.pl                        %LOC%\CreatePage.pl
copy Archive\NWTC_Lib_v%1.exe             %LOC%\NWTC_Lib_v%1.exe
copy Archive\NWTC_Lib_v%1.tar.gz          %LOC%\NWTC_Lib_v%1.tar.gz
copy Archive\NWTC_Lib_v%1_Maint.exe       %LOC%\NWTC_Lib_v%1_Maint.exe

perl %LOC%\CreatePage.pl

set LOC=

:Done
