% Before running this script, you must have compiled OpenFAST for Simulink to create a DLL (i.e., a shared library like .so, .dylib, .lib, etc.).
% The name of the library that was generated must match the `libname` variable below, and should be located in the directory specified by `binDir`.

mexname = 'FAST_SFunc';

switch computer('arch')
    case 'win64'
        % this is set up for files generated using the x64 configuration of vs-build
        libDir = '../../../build/bin';
        includeDir = '/usr/local/include';
        libName = 'OpenFAST-Simulink_Win32';
        outDir = libDir;
        
    case 'win32'
        % this is set up for files generated using the x86 configuration of vs-build
        libDir = '../../../build/bin';
        includeDir = '/usr/local/include';
        libName = 'OpenFAST-Simulink_Win32';
        outDir = libDir;
        
    case 'maci64'
        libDir = '/usr/local/lib';
        includeDir = '/usr/local/include';
        libName = 'openfastlib';
        outDir = '.';

    otherwise
        error('Unexpected computer architecture type.')
end

mex('-largeArrayDims', ...
    ['-L' libDir], ...
    ['-l' libName], ...
    ['-I' includeDir], ...
    '-outdir', outDir, ...
    'COMPFLAGS=$COMPFLAGS -MT -D', ...
    ['S_FUNCTION_NAME=' mexname], ...
    '-output', mexname, ...
    'FAST_SFunc.c');
