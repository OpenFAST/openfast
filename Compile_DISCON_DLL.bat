@ECHO OFF

REM ----------------------------------------------------------------------------
REM                   set compiler internal variables 
REM ----------------------------------------------------------------------------
REM    You can run this bat file from the IVF compiler's command prompt (and not 
REM    do anything in this section). If you choose not to run from the IVF command
REM    prompt, you must call the compiler's script to set internal variables.
REM    TIP: Right click on the IVF Compiler's Command Prompt shortcut, click
REM    properties, and copy the target (without cmd.exe and/or its switches) here:

ECHO.
 
IF "%INCLUDE%"=="" ( 

IF /I "%1"=="64" ( ECHO //// Using intel64 compiler \\\\
CALL "C:\Program Files (x86)\Intel\Composer XE 2011 SP1\bin\ipsxe-comp-vars.bat" intel64 vs2010
) ELSE ( ECHO  //// Using ia32 compiler \\\\
CALL "C:\Program Files (x86)\Intel\Composer XE 2011 SP1\bin\ipsxe-comp-vars.bat" ia32 vs2010
)

) ELSE ( ECHO //// Using existing compiler settings \\\\
)

REM ----------------------------------------------------------------------------
REM                   set compiler and linker options
REM ----------------------------------------------------------------------------

REM Use /libs:static if the DLL will be used on other PCs. (will result in a 
REM   larger dll)

rem SET CompOpts= /nologo /O2 /inline /traceback /libs:static /threads 
SET CompOpts= /nologo /O2 /inline /traceback /libs:dll /threads 
rem SET LinkOpts=/link

SET FileName=DISCON

REM ----------------------------------------------------------------------------
REM                   compile DLL
REM ----------------------------------------------------------------------------

ECHO.
ECHO Creating %FileName%.dll:

CALL IFORT /dll %CompOpts% %FileName%.f90 /exe:%FileName%.dll

ECHO.


REM ----------------------------------------------------------------------------
REM                   clear variables
REM ----------------------------------------------------------------------------

DEL %FileName%.obj

SET CompOpts=
SET LinkOpts=
SET FileName=



