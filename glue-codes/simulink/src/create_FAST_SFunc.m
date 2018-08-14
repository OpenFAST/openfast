
switch computer('arch')
    case 'win64'
        mex -v -L../../bin -lFAST_Library_x64 ...
            -I../../Source -I../../Source/dependencies/OpenFOAM -outdir ../../bin COMPFLAGS='$COMPFLAGS /MT' FAST_SFunc.c        
    case 'win32'
        mex -L../../bin -lFAST_Library_Win32 ...
            LINKFLAGS='$LINKFLAGS /LARGEADDRESSAWARE' ...
            -I../../Source -I../../Source/dependencies/OpenFOAM -outdir ../../bin COMPFLAGS='$COMPFLAGS /MT' FAST_SFunc.c
    case 'maci64'
        mex -L/usr/local/lib -lopenfastlib ...    
            -I/usr/local/include -outdir . COMPFLAGS='$COMPFLAGS -MT' FAST_SFunc.c 
            % This builds the mex file in the current directory, ./openfast/glue-codes/simulink/src/
    otherwise
        error('Unexpected computer architecture type.')
end

%%
 %addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\bin');
 %addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\Simulink\Samples');
