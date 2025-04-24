@ECHO off
SET IncludeFile=..\gitVersionInfo.h

<NUL SET /p IncludeTxt=#define GIT_VERSION_INFO '> %IncludeFile%
FOR /f %%a IN ('git describe --abbrev^=8 --always --tags --dirty') DO <NUL SET /p IncludeTxt=%%a>> %IncludeFile%
git describe --abbrev^=8 --always --tags --dirty > NUL
IF %ERRORLEVEL%==0 ( ECHO '>> %IncludeFile% ) else ( ECHO Unversioned from $Format:%H$ '>> %IncludeFile% )

EXIT /B 0