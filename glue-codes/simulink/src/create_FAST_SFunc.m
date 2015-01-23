mex -L../../bin -lFAST_Library_x64.lib -compatibleArrayDims -outdir ../../bin FAST_SFunc.c 
% mex FAST_gateway.c -l..\..\bin\FAST_Library_x64.lib
%%
FAST_InputFileName = '..\..\CertTest\Test03.fst';
TMax = 20;

 addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\bin');