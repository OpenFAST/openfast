@ECHO off
SETLOCAL ENABLEDELAYEDEXPANSION
SET IncludeFile=..\gitVersionInfo.h

<NUL SET /p IncludeTxt=#define GIT_VERSION_INFO '> %IncludeFile%
FOR /f %%a IN ('git describe --abbrev^=8 --always --tags --dirty 2^>NUL') DO <NUL SET /p IncludeTxt=%%a>> %IncludeFile%
git describe --abbrev^=8 --always --tags --dirty > NUL 2>&1
IF !ERRORLEVEL!==0 (
    ECHO '>> %IncludeFile%
) ELSE (
    REM git describe failed (e.g. shallow clone). Fall back to short commit hash.
    SET GitHash=
    FOR /f %%a IN ('git rev-parse --short^=8 HEAD 2^>NUL') DO (
        SET GitHash=%%a
    )
    IF DEFINED GitHash (
        git diff-index --quiet HEAD -- > NUL 2>&1
        IF !ERRORLEVEL!==0 (
            <NUL SET /p IncludeTxt=!GitHash!>> %IncludeFile%
        ) ELSE (
            <NUL SET /p IncludeTxt=!GitHash!-dirty>> %IncludeFile%
        )
        ECHO '>> %IncludeFile%
    ) ELSE (
        ECHO Unversioned from $Format:%H$ '>> %IncludeFile%
    )
)

ENDLOCAL
EXIT /B 0