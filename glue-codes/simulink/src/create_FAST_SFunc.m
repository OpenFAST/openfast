%% 32-bit Matlab:
% mex -L../../bin -lFAST_Library_Win32 ...
%     LINKFLAGS='$LINKFLAGS /STACK:999999999 /LARGEADDRESSAWARE' ...
%     -I../../Source -I../../Source/dependencies/OpenFOAM -outdir ../../bin COMPFLAGS='$COMPFLAGS /MT' FAST_SFunc.c

%% 64-bit Matlab:
mex -v -L../../bin -lFAST_Library_x64 ...
    -I../../Source -I../../Source/dependencies/OpenFOAM -outdir ../../bin COMPFLAGS='$COMPFLAGS /MT' FAST_SFunc.c

%%
 %addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\bin');
 %addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\Simulink\Samples');
