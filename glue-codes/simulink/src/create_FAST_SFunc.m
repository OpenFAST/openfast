
switch computer('arch')
    case 'win64'
        mex -v -L../../bin -lFAST_Library_x64 ...
            -I../../Source -I../../Source/dependencies/OpenFOAM -outdir ../../bin COMPFLAGS='$COMPFLAGS /MT' FAST_SFunc.c        
    case 'win32'
        mex -L../../bin -lFAST_Library_Win32 ...
            LINKFLAGS='$LINKFLAGS /LARGEADDRESSAWARE' ...
            -I../../Source -I../../Source/dependencies/OpenFOAM -outdir ../../bin COMPFLAGS='$COMPFLAGS /MT' FAST_SFunc.c
    otherwise
        error('Unexpected computer architecture type.')
end

%%
 %addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\bin');
 %addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\Simulink\Samples');
