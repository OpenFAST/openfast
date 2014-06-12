@ECHO OFF

set LOC=Y:\Wind\WindWeb\designcodes\miscellaneous\NWTC_Library

@IF /I "%1"=="-NOFILES" GOTO UpdateWeb
@IF NOT "%1"=="" GOTO CopyFiles

@ECHO 
@ECHO  The syntax for updating files on the NWTC Subroutine Library page is "updateweb <version>"
@ECHO.
@ECHO  Example:  "updateweb 1.01.08"
@ECHO.
@ECHO  The syntax for updating the web page without copying the other files is "updateweb -nofiles"
@ECHO.

@GOTO Done

:CopyFiles
copy ChangeLog.txt                            %LOC%\ChangeLog.txt
copy Documentation\NWTC_Library_Routines.txt  %LOC%\NWTC_Library_Routines.txt
copy Archive\NWTC_Lib_v%1.exe                 %LOC%\NWTC_Lib_v%1.exe
copy Archive\NWTC_Lib_v%1.tar.gz              %LOC%\NWTC_Lib_v%1.tar.gz

:UpdateWeb
copy CreatePage.pl                            %LOC%\CreatePage.pl
perl %LOC%\CreatePage.pl

:Done
set LOC=

