@ECHO OFF

set LOC=\\netapp01\websites\wind\designcodes\simulators\fast\verification

copy CreatePageVerif.pl    %LOC%\CreatePageVerif.pl
copy CertTest\Test01.pdf   %LOC%\Test01.pdf
copy CertTest\Test02.pdf   %LOC%\Test02.pdf
copy CertTest\Test03.pdf   %LOC%\Test03.pdf
copy CertTest\Test04.pdf   %LOC%\Test04.pdf
copy CertTest\Test05.pdf   %LOC%\Test05.pdf
copy CertTest\Test06.pdf   %LOC%\Test06.pdf
copy CertTest\Test07.pdf   %LOC%\Test07.pdf
copy CertTest\Test08.pdf   %LOC%\Test08.pdf
copy CertTest\Test09.pdf   %LOC%\Test09.pdf
copy CertTest\Test10.pdf   %LOC%\Test10.pdf
copy CertTest\Test11.pdf   %LOC%\Test11.pdf
copy CertTest\Test12.pdf   %LOC%\Test12.pdf
copy CertTest\Test13.pdf   %LOC%\Test13.pdf
copy CertTest\Test14.pdf   %LOC%\Test14.pdf

perl %LOC%\CreatePageVerif.pl

set LOC=

:Done
