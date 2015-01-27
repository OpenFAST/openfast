% make sure the FASTv8\bin directory is in the MATLAB path
%    (relative path names are not recommended in addpath()):
% addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\bin');


CertTest_Dir = '..\..\CertTest';

 for iTest = 4 % [1 3:13 15:17]
    
        %------------------------------------------------------------------       
        % Set up and run the Simulink OpenLoop model
        %------------------------------------------------------------------       
    
    FileRoot   = sprintf( 'Test%02.0f', iTest );
    
    disp('***********************************************');
    disp( ['FAST_SFunc certification test for ' FileRoot] );
    disp('***********************************************');
    
    FAST_InputFileName = [CertTest_Dir filesep FileRoot '.fst'];
    TMax               = 20;
    
    sim('OpenLoop.mdl',[0,TMax]);
       
 end