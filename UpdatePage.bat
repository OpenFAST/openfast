@ECHO OFF

set LOC=Y:\Wind\WindWeb\designcodes\miscellaneous\nwtc_subs

copy CreatePage.pl %LOC%\CreatePage.pl

perl %LOC%\CreatePage.pl

set LOC=

:Done
