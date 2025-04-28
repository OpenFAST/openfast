%% INSTRUCTIONS
% This script is used to manually build a Simulink mex file on Windows with Visual Studio. It uses the openfastlib shared
% library (.dll).
%
% If you are using Windows and building with CMake, do not use this script.  Instead use cmake to build the FAST_SFunc.mexXXXX directly.
%
% If you are not using Windows, do not use this script.  Instead use cmake to build the FAST_SFunc.mexXXXX directly.
%
% Alternative building with CMAKE:
% specify -DBUILD_OPENFAST_SIMULINK_API=ON when running cmake
%  - "make FAST_SFunc" will place the resulting mex file at <build-dir>/glue-codes/simulink/FAST_SFunc.mexXXXX
%  - "make install" will place the resulting mex file at install/bin/FAST_SFunc.mexXXXX
%
%
% Before running this script, you must have compiled OpenFAST for Simulink to create a DLL (i.e., a shared library .lib).
% - If the Visual Studio Solution file contained in the vs-build directory was used to create the DLL on Windows,
%   make sure `built_with_visualStudio` is set to true.
% - The name of the library that was generated must match the `libname` variable below
%   and should be located in the directory specified by `libDir`.
% - The `includeDir` variable must specify the directory that contains the following header files:
%   "FAST_Library.h", "OpenFOAM_Types.h", and "ExtLoadsDX_Types.h" 
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

   fprintf( '\n----------------------------\n' );
   fprintf( 'Do not use this script with Mac/Linux.  Follow the CMake instructions at the top of the script instead.' );
   fprintf( '\n----------------------------\n' );

end



if ispc () % Windows PC
   %% BUILD COMMAND
   fprintf( '\n----------------------------\n' );
   fprintf( 'Creating %s\n\n', [outDir filesep mexname '.' mexext] );


    mex('-largeArrayDims', ...
        ... % '-v', ... %add this line for "verbose" output (for debugging)
        ['-L' libDir], ...
        ['-l' libName], ...
        ['-I' includeDir], ...
        '-I../../../modules/externalinflow/src',  ... % needed for visual studio builds to find "ExternalInflow_Types.h"
        '-I../../../modules/extloads/src', ... % needed for visual studio builds to find "ExtLoadsDX_Types.h"
        '-outdir', outDir, ...
        ['COMPFLAGS=$COMPFLAGS -MT -DS_FUNCTION_NAME=' mexname], ...
        '-output', mexname, ...
        'FAST_SFunc.c');

end
