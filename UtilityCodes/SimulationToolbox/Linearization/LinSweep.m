%Linearization Sweep Using Simulation Toolbox
%Paul Fleming, July 2011
%Adapted from code written by Jan-Willem van Wingerden of TU DELFT in summer 2010

clc; clear;

%% Parameter specification
fstFile = 'NREL_G97.fst'; %FST file to use, ensure it is already set up for linearizations, the script will be adjusting only the reference rotor speeds
iptFile = 'dlc12g_1_AD.ipt'; %The AD file, also ensure this is configured for linearization, just the con_wind speeds will be adjusted
linFile2 = 'lin_file2'; %The linearization input file for region 2,
linFile3 = 'lin_file3'; %The linearization input file for region 3,
windSweep = [0 5:5:25]; %An array of wind speeds to sweep across, increments can't be smaller than 0.1 m/s
ratedWindSpeed = 9.5; %Wind speed where operating switches from region2 to region3
bladeFinePitch = 0.74; %Fine pitch angle of blades
baseFileName = 'CART3Sweep';



%% First generate the wind files for use in this sweep
disp('Generating Constant Wind Speed Files')

if ~isdir('Wind')
    system('mkdir Wind');% make the wind subdirectory if it doesn't already exist
end
cd Wind %enter wind file directory
for i = 1:length(windSweep)
    windList{i} = GenConstWindFile('ConWind',windSweep(i));
end

disp('Generating Step Wind Speed File')
[junk startReg3] = find(windSweep>=ratedWindSpeed,1);
numSteps = startReg3 - 1; %numSteps in region 2
if numSteps >= 1
    stepWindFile = 'StepWind.wnd';
    GenStepWindFile(stepWindFile,windSweep(1:numSteps));
else
    disp('No step test needed (all region 3)')
end


cd .. %exit wind file directory

%% Now run the step change wind file to get the RPM set points for sweep
if numSteps >=1

%Make an FST file for the wind sweep
stepFSTFile = 'step.fst';
stepIPTFile = 'step.ipt';

FP = Fast2Matlab(fstFile,4);

%Modify parameters for this sweep
FP = SetFastPar(FP,'AnalMode',1);%Set to run time marching simulation
FP = SetFastPar(FP,'BlPitch(1)',bladeFinePitch);%Set blades to fine pitch
FP = SetFastPar(FP,'BlPitch(2)',bladeFinePitch);%Set blades to fine pitch
FP = SetFastPar(FP,'BlPitch(3)',bladeFinePitch);%Set blades to fine pitch
FP = SetFastPar(FP,'ADFile',['"' stepIPTFile '"']);%Set IPT file for this step test
FP = SetFastPar(FP,'TMax',100 * numSteps);%Set TMax file for this step test
FP = SetFastPar(FP,'DT',.005);%1 Hz seems reasonable here
FP = SetFastPar(FP,'RotSpeed',10);%1 Hz seems reasonable here

%Setup the outlist to have certain outputs in certain columns
FP.OutList{1} = '"HSShftTq"';
FP.OutList{2} = '"LSSTipVxa"';
FP.OutList{3} = '"WindVxi"';
FP.OutList{4} = '"TSR"';

%Generate new FST file
Matlab2FAST(FP,fstFile,stepFSTFile);
%Generate a new IPT file
ADP = AeroDyn2Matlab(iptFile,1);
%Modify parameters for this sweep
ADP = setFastPar(ADP,'WindFile',['"Wind\' stepWindFile '"']);%Select step wind file
Matlab2AD(ADP,iptFile,stepIPTFile);

%now can run the time marching simulation
system(['fast ' stepFSTFile]);

%now read in the output taking advantage of the fact that the outputs were
%setup to be in specific locations
importfile([stepFSTFile(1:end-3) 'out']);
t=data(:,1);
HSS_torque=data(:,2)*1000; %get into NM
LSS_RPM=data(:,3);
V=data(:,4);
TSR = data(:,5);

figure
title('Results of step testing')
a(1)=subplot(4,1,1);
plot(t,HSS_torque)
legend('HSS Torque (Nm)')
a(2)=subplot(4,1,2);
plot(t,LSS_RPM)
legend('LSS RPM (RPM)')
a(3)=subplot(4,1,3);
plot(t,V)
legend('Wind Speed (m/s)')
a(4)=subplot(4,1,4);
plot(t,TSR)
legend('TSR (-)')
linkaxes(a,'x')

%% Find operationg points for region2

l = length(LSS_RPM);
stepLength = floor(l/numSteps);
for i = 1:numSteps
    testRange = ((i-1) * stepLength + floor(stepLength *3/4)):(i *stepLength);
    RPMSetPoint(i) = mean(LSS_RPM(testRange));
    if windSweep(i) == 0
        RPMSetPoint(i) = 0 %exception for 0 wind case
    end
end

end
%% Derive the trim list, which says whether a given linearization is region
% 2 or region 3
trimCase = ones(size(windSweep));
for i = 1:length(trimCase)
    if windSweep(i) >= ratedWindSpeed
        trimCase(i) = 3;
    else
        trimCase(i) = 2;
    end
end

%% Now sweep across the wind speeds

sweepFSTFile = 'sweep.fst';
sweepIPTFile = 'sweep.ipt';
linResultFolder = 'LinResult';


linFile = [sweepFSTFile(1:end-3) 'lin'];
if ~isdir(linResultFolder)
    system(['mkdir ' linResultFolder]);% make the wind subdirectory if it doesn't already exist
end


for i=1:length(windSweep)

    disp(['Running wind case: ' windList{i}])
    disp('---------------------------------------------------')
    
    %Generate a new FST file
    FP = FAST2Matlab(fstFile,4);

    %Modify parameters for this sweep
    FP = setFastPar(FP,'AnalMode',2);%Set to find linearization

    %Set rotor speed and trim case
    if trimCase(i) == 2
        FP = setFastPar(FP,'RotSpeed',RPMSetPoint(i));%Set to region 2
        FP = setFastPar(FP,'LinFile',linFile2);
    else
        FP = setFastPar(FP,'RotSpeed',37.0679);%Set to region 3
        FP = setFastPar(FP,'LinFile',linFile3);
    end
    
    %Set to match IPTfile
    FP = setFastPar(FP,'ADFile',sweepIPTFile);
    
    %Write new FST file
    Matlab2FAST(FP,fstFile,sweepFSTFile);
    
    %Generate new IPT file
    ADP = AeroDyn2Matlab(iptFile,1);
    %Modify parameters for this sweep
    ADP = setFastPar(ADP,'WindFile',['"Wind\' windList{i} '"']);%Select step wind file
    Matlab2AD(ADP,iptFile,sweepIPTFile);
    
    %run Fast
    disp('Running linearization')
    system(['fast ' sweepFSTFile]);
    
    
    %following the sweep, move the lin file into a result direct
    linList{i} = [linResultFolder '\' linFile(1:end-4) '_' windList{i} '.lin'];
    system(['copy ' linFile ' ' linList{i} ]);
end



