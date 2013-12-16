:: Some path settings relative to here
@SET LocationDocs="..\Docs"
@SET LocationArchives="Archive"

@ECHO OFF

@IF NOT "%1"==""  GOTO Copy

@ECHO 
@ECHO  The syntax for updating files on the SubDyn page is "update_web 0.04.00a-rrd"
@ECHO.
@ECHO  Example:  "update_web 0.04.00a-rrd"

@GOTO Done

:Copy

set LOC=Y:\Wind\WindWeb\designcodes\simulators\SubDyn
set ProgName=SubDyn

copy CreatePage.pl                                                   %LOC%\CreatePage.pl

copy %LocationArchives%\%ProgName%_v%1.exe                           %LOC%\%ProgName%_v%1.exe
copy %LocationArchives%\%ProgName%_v%1.tar.gz                        %LOC%\%ProgName%_v%1.tar.gz
copy %LocationDocs%\SUBDYN_README_v0.04.00a-rrd.pdf                  %LOC%\SUBDYN_README_v0.04.00a-rrd.pdf
copy ChangeLog.txt                                                   %LOC%\ChangeLog.txt

perl %LOC%\CreatePage.pl

set LOC=

:Done
