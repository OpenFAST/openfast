@ECHO OFF

set LOC=Y:\Wind\WindWeb\nwtc\docs

:CopyFiles
copy docs\Manual\BeamDyn_Manual.pdf             %LOC%\BeamDyn_Manual.pdf

:Done
set LOC=
