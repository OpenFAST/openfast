%% INSTRUCTIONS
% Before running this script, you must have compiled OpenFAST for Simulink to create a DLL (i.e., a shared library like .so, .dylib, .lib, etc.).
% - If cmake was used, make sure the install directory is specified properly in the `installDir` variable below,
%   and set `built_with_visualStudio` to false (necessary on Windows only).
% - If the Visual Studio Solution file contained in the vs-build directory was used to create the DLL on Windows,
%   make sure `built_with_visualStudio` is set to true.
% - The name of the library that was generated must match the `libname` variable below
%   and should be located in the directory specified by `libDir`.
% - The `includeDir` variable must specify the directory that contains the following header files:
%   "FAST_Library.h", "OpenFOAM_Types.h", and "SuperController_Types.h"
% - The `outDir` variable indicates where the resulting mex file will reside.
%
% Run `mex -setup` in Matlab to configure a C compiler if you have not already done so.


mexname = 'FAST_SFunc'; % base name of the resulting mex file

built_with_visualStudio = true; %if the libraries were built with cmake, set to false


if (ispc && built_with_visualStudio)   
%% defaults for visual studio builds:

    libDir = '../../../build/bin';
    includeDir = '../../../modules/openfast-library/src';  % needed for visual studio builds to find "FAST_Library.h"
    outDir = libDir;
        
    switch computer('arch')
        case 'win64'
            % this is set up for files generated using the x64 configuration of vs-build
            libName = 'OpenFAST-Simulink_x64';

        case 'win32' 
            % this is set up for files generated using the x86
            % configuration of vs-build (win32 will work only on older versions of Matlab)
            libName = 'OpenFAST-Simulink_Win32';
    end
    
else    
%% defaults for cmake builds:

    if ( ispc ) % Windows PC
        installDir = '../../../install';
        outDir = fullfile(installDir, 'lib');
        % If there are shared libraries does it work for outDir to be the local directory?
    else
        installDir = '/usr/local';
        outDir = '.';
    end

    libDir = fullfile(installDir, 'lib');
    includeDir = fullfile(installDir, 'include');
    libName = 'openfastlib';
end

%% BUILD COMMAND
fprintf( '\n----------------------------\n' );
fprintf( 'Creating %s\n\n', [outDir filesep mexname '.' mexext] );

mex('-largeArrayDims', ...
    ['-L' libDir], ...
    ['-l' libName], ...
    ['-I' includeDir], ...
    '-I../../../modules/supercontroller/src', ... % needed for visual studio builds to find "SuperController_Types.h"
    '-I../../../modules/openfoam/src',        ... % needed for visual studio builds to find "OpenFOAM_Types.h"
    '-outdir', outDir, ...
    'COMPFLAGS=$COMPFLAGS -MT -D', ...
    ['S_FUNCTION_NAME=' mexname], ...
    '-output', mexname, ...
    'FAST_SFunc.c');
