
mexname = 'FAST_SFunc';

switch computer('arch')
    case 'win64'
        binDir = '../../../build/bin';
        libname = 'OpenFAST-Simulink_x64';
        outDir = binDir;
        
    case 'win32'
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