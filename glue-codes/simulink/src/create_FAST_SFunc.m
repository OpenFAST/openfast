mex -L../../bin -lFAST_Library_x64.lib -compatibleArrayDims -outdir ../../bin FAST_SFunc.c 
% mex FAST_gateway.c -l..\..\bin\FAST_Library_x64.lib
%%
% FAST_InputFileName = '..\..\CertTest\Test03.fst';
% TMax = 20;

 addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\bin');
 addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\Simulink\Samples');
 
 %%
CertTest_Dir = '..\..\CertTest';

 for iTest = 4 % [1 3:13 15:17]
    
        %------------------------------------------------------------------       
        % Set up and run the Simulink OpenLoop model
        %------------------------------------------------------------------       
%     clear FAST_Sfunc;       %perhaps this requirement could be removed in the future, after we have cleaned up all the SAVEd variables and INITIALIZATION of variables in modules, and deal with AeroDyn's internal states.
    
    FileRoot   = sprintf( 'Test%02.0f', iTest );
    
    disp('***********************************************');
    disp( ['FAST_SFunc certification test for ' FileRoot] );
    disp('***********************************************');
    
    FAST_InputFileName = [CertTest_Dir filesep FileRoot '.fst'];
    TMax       = 20;
    
    sim('OpenLoop.mdl',[0,TMax]);
   
    
end
