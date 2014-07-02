@ECHO OFF

set LOC=Y:\Wind\WindWeb\designcodes\simulators\FAST8

@IF /I "%1"=="-NOFILES" GOTO UpdateWeb
@IF NOT "%1"=="" GOTO CopyFiles

@ECHO 
@ECHO  The syntax for updating files on the FAST page is "updateweb <version>"
@ECHO.
@ECHO  Example:  "updateweb 7.01.00"
@ECHO.
@ECHO  The syntax for updating the web page without copying the other files is "updateweb -nofiles"
@ECHO.

@GOTO Done

:CopyFiles
copy ChangeLog.txt             %LOC%\ChangeLog.txt
copy Archive\FAST_v%1.exe      %LOC%\FAST_v%1.exe
copy Archive\FAST_v%1.tar.gz   %LOC%\FAST_v%1.tar.gz
copy Archive\FAST_v%1_all.exe  %LOC%\FAST_v%1_all.exe
copy README_FAST8.pdf          %LOC%\README_FAST8.pdf

:UpdateWeb
copy CreatePage.pl             %LOC%\CreatePage.pl
perl %LOC%\CreatePage.pl

:Done
set LOC=
