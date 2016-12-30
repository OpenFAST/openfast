% make sure the FASTv8\bin directory is in the MATLAB path
%    (relative path names are not recommended in addpath()):
% cd ..\..\
% FASTv8_root_directory = pwd;
% cd Simulink\Samples
% addpath([ FASTv8_root_directory '\bin']);


% Simple Induction Generator Example ======================================
% To model a simple induction generator in Simulink use model Test01_SIG.mdl.  
% The following parameters duplicate those used in Certification Test #01.  

% Change Test01.fst as follows:
% set VSContrl = 4 in ..\..\CertTest\AWT27\Test01_ServoDyn.dat


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
FAST_InputFileName = '..\..\CertTest\Test01.fst';
TMax = 20;

% run the model
sim('Test01_SIG.mdl',[0,TMax]);

% look at results:
% PlotFASToutput({'../../CertTest/Test01.SFunc.out','../../CertTest/Test01.out'},{'SFunc','exe'});