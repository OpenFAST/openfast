% Before running this script, you must have compiled OpenFAST for Simulink to create a DLL (i.e., a shared library like .so, .dylib, .lib, etc.).
% The name of the library that was generated must match the `libname` variable below, and should be located in the directory specified by `binDir`.

mexname = 'FAST_SFunc';

switch computer('arch')
    case 'win64'
        % this is set up for files generated using the x64 configuration of vs-build
        binDir = '../../../build/bin';
        libname = 'OpenFAST-Simulink_x64';
        outDir = binDir;
        
    case 'win32'
        % this is set up for files generated using the x86 configuration of vs-build
        binDir = '../../../build/bin';
        libname = 'OpenFAST-Simulink_Win32';
        outDir = binDir;
        
    case 'maci64'
        binDir = '/usr/local/lib';
        libname = 'openfastlib';
        outDir = '.';
        
%         mex -L/usr/local/lib -lopenfastlib ...    
%             -I/usr/local/include -outdir . COMPFLAGS='$COMPFLAGS -MT' FAST_SFunc.c 
            % This builds the mex file in the current directory, ./openfast/glue-codes/simulink/src/
    otherwise
        error('Unexpected computer architecture type.')
end

    mex( '-largeArrayDims', ['-L' binDir], ['-l' libname] ,...
        '-I../../../modules-local/openfast-library/src', '-I../../../build/types-files', ['-I' binDir], '-outdir', outDir, ...
        ['COMPFLAGS=$COMPFLAGS /MT /D S_FUNCTION_NAME=' mexname], '-output', mexname, ...
        'FAST_SFunc.c');


%%
 %addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\bin');
 %addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\Simulink\Samples');
% use mex -setup to configure a C compiler