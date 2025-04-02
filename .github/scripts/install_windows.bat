REM SPDX-FileCopyrightText: 2022 Intel Corporation
REM
REM SPDX-License-Identifier: MIT

set URL=%1
set COMPONENTS=%2

curl.exe --output %TEMP%\webimage.exe --url %URL% --retry 5 --retry-delay 5
start /b /wait %TEMP%\webimage.exe -s -x -f webimage_extracted --log extract.log
del %TEMP%\webimage.exe
webimage_extracted\bootstrapper.exe -s --action install --components=%COMPONENTS% --eula=accept -p=NEED_VS2017_INTEGRATION=0 -p=NEED_VS2019_INTEGRATION=0 -p=NEED_VS2022_INTEGRATION=1 --log-dir=.

set installer_exit_code=%ERRORLEVEL%
rd /s/q "webimage_extracted"
exit /b %installer_exit_code%
