mex -L../../bin -lFAST_Library_Win32 -outdir ../../bin FAST_SFunc.c 
%mex -L../../bin -lFAST_Library_x64 -outdir ../../bin FAST_SFunc.c 

%%
 %addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\bin');
 %addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\Simulink\Samples');
  
% *************************************************************************** 
%   Warning: MEX-files generated using Microsoft Visual C++ 2010 require 
%            that Microsoft Visual Studio 2010 run-time libraries be  
%            available on the computer they are run on. 
%            If you plan to redistribute your MEX-files to other MATLAB 
%            users, be sure that they have the run-time libraries. 
% ***************************************************************************  