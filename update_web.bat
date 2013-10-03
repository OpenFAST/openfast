@ECHO OFF

set LOC=y:\wind\WindWeb\designcodes\simulators\inflowwind\alpha

@IF /I "%1"=="-NOFILES" GOTO UpdateWeb
@IF NOT "%1"=="" GOTO CopyFiles

@ECHO 
@ECHO  The syntax for updating files on the InflowWind page is "update_web <version>"
@ECHO.
@ECHO  Example:  "update_web 1.03.01"
@ECHO.
@ECHO  The syntax for updating the web page without copying the other files is "update_web -nofiles"
@ECHO.

@GOTO Done

:CopyFiles
copy ChangeLog.txt                 %LOC%\ChangeLog.txt
copy Archive\InflowWind_v%1.exe    %LOC%\InflowWind_v%1.exe
copy Archive\InflowWind_v%1.tar.gz %LOC%\InflowWind_v%1.tar.gz
copy InflowWind.pdf                %LOC%\InflowWind.pdf


:UpdateWeb
copy CreatePage.pl                 %LOC%\CreatePage.pl
perl %LOC%\CreatePage.pl

:Done
set LOC=


