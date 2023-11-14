% make sure the OpenFAST directory where the FAST_SFunc.mex* file is located
% is in the MATLAB path (also make sure any other OpenFAST library files that
% are needed are on the MATLAB path)
%    (relative path names are not recommended in addpath()):
% addpath('../../../build/bin'); % install location for Windows Visual Studio builds
% addpath(genpath('../../../install')); % cmake default install location


% Simple Induction Generator Example ======================================
% To model a simple induction generator in Simulink use model Test01_SIG.mdl.  
% The following parameters duplicate those used in Certification Test #01.  

% Change AWT_YFix_WSt.fst (formerly Test01.fst) as follows:
% set VSContrl = 4 in ../../../reg_tests/r-test/glue-codes/openfast/AWT_YFix_WSt/AWT_YFix_WSt_ServoDyn.dat


GenEff   =  100.0;          % - Generator efficiency [ignored by the Thevenin and user-defined generator models] (%)
GBRatio  =   22.5;          % - Gearbox ratio (-)
SIG_SlPc =    1.5125;       % - Rated generator slip percentage [>0] (%)              Now HSS side!
SIG_SySp = 1200.0;          % - Synchronous (zero-torque) generator speed [>0] (rpm)  Now HSS side!
SIG_RtTq = 1367.9;          % - Rated torque [>0] (N-m)                               Now HSS side!
SIG_PORt =    2.0;          % - Pull-out ratio (Tpullout/Trated) [>1] (-)% 
SIG_SySp = SIG_SySp*pi/30;  % convert to rad/s
SIG_RtSp = SIG_SySp*(1.0+0.01*SIG_SlPc);
SIG_POS1=SIG_PORt*(SIG_RtSp-SIG_SySp);
SIG_POTq=SIG_RtTq*SIG_PORt;
SIG_Slop=SIG_RtTq/(SIG_RtSp - SIG_SySp);


% parameters required for the S-Function block:
OpenFASTRoot = '../../../reg_tests/r-test/glue-codes/openfast/AWT_YFix_WSt/';
FAST_InputFileName = [ OpenFASTRoot 'AWT_YFix_WSt.fst' ];
TMax = 20;

OutList = {'Time','Wind1VelX','Wind1VelY','Wind1VelZ','LSSGagVxa','LSSGagPxa','TeetDefl',...
           'TipDxb2','TipDyb2','TipALxb2','TipALyb2','Spn2ALxb1','Spn2ALyb1','YawBrRDxt',...
           'YawBrRDyt','YawBrRVxp','YawBrRVyp','YawBrRAxp','YawBrRAyp','RootMyc1','RootMxc1',...
           'RootFxc2','RootFyc2','Spn3MLxb1','Spn3MLyb1','RotTorq','YawBrMzn','TwrBsMzt'};

% run the model
sim('Test01_SIG.mdl',[0,TMax]);

% look at results:
% PlotFASToutput({[ OpenFASTRoot 'AWT_YFix_WSt.SFunc.out'],[ OpenFASTRoot 'windows-intel/AWT_YFix_WSt.out']},{'SFunc','exe'});