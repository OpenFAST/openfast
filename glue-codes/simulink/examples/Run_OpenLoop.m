% make sure the FASTv8\bin directory is in the MATLAB path
%    (relative path names are not recommended in addpath()):
% addpath('C:\Users\bjonkman\Documents\CAETools\FASTv8\bin');


CertTest_Dir = '..\..\CertTest';

CertTest_TMax=[20, 20, 20, 70, 30, ...
               35, 70, 20, 40, 25, ...
               20, 20, 40,  0, 20, ...
               20, 70, 60, 60, 60, ...
               60, 60, 60, 60, 60, ...
               20                  ];

 for iTest = [1:13 15:26]  
     
        %------------------------------------------------------------------       
        % Set up and run the Simulink OpenLoop model
        %------------------------------------------------------------------       
    
    FileRoot   = sprintf( 'Test%02.0f', iTest );
    
    disp('***********************************************');
    disp( ['FAST_SFunc certification test for ' FileRoot] );
    disp('***********************************************');
    
    FAST_InputFileName = [CertTest_Dir filesep FileRoot '.fst'];
    TMax               = CertTest_TMax(iTest);
    
    sim('OpenLoop.mdl',[0,TMax]);
       
 end