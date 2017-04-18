function PlotResults(azimuth)

WTP_Results_File = 'C:\Dev\NREL_SVN\WT_Perf\branches\v4.x\CertTest\ccBlade_UAE.WTP.out';
% Setup data which both files have in common
%r = [0.5680,0.8801,1.2321,1.5087,1.7099,1.9278,2.1457,2.3469,2.5480,2.7660,2.9840,3.1850,3.3862,3.6041,3.8220,4.0232,4.2244,4.4004,4.5764,4.7776,4.9536];
r             = [0.660, 0.883, 1.008, 1.067, 1.133, 1.257, 1.343, 1.510, 1.648,1.952, 2.257, 2.343, 2.562, 2.867, 3.172, 3.185, 3.476, 3.781, 4.023, 4.086, 4.391,4.696, 4.780, 5.000];
r_FAST        = [0.6060,0.8801,1.2321,1.5087,1.7099,1.9278,2.1457,2.3469,2.5480,2.7660,2.9840,3.1850,3.3862,3.6041,3.8220,4.0232,4.2244,4.4004,4.5764,4.7776,4.9536];

R             = 5.029;
r_over_R      = r ./ R;
r_over_R_FAST = r_FAST ./ R;
yawAngles     = [0.0, 20.0, 30.0, 60.0, 75.0, 90.0];
windSpeeds    = [5.0, 7.0, 10.0];
steps_per_rev = 16;
numYaw        = 6;
numWindSpeeds = length(windSpeeds);
numNodes   = length(r);
nFASTNodes = length(r_FAST);
psi = [0:steps_per_rev-1]*360/steps_per_rev;


% find the first time step associated with the last
% revolution prior to the ramping wind speed


%rev4_plus_one = steps_per_rev*4 + 1;

[ FAST_AOA, FAST_Fx, FAST_Fy ] = Read_FAST_Results(azimuth);

[ WTP_AOA, WTP_Fx, WTP_Fy ] = Read_WTP_Results(WTP_Results_File);

[ r_over_R_EXP, EXP_Fx, EXP_Fy ] = Read_Experimental_Results(windSpeeds, yawAngles, azimuth);

dataCCblade = Read_CCblade_Results(numNodes,numYaw,numWindSpeeds, azimuth);


%Chans = strmatch('AOA',ChanName);

switch azimuth
    case 0.0
        offset = 1;
    case 90.0
        offset =   steps_per_rev / 4 + 1;
    case 180.0
        offset =   steps_per_rev / 2 + 1;
    case 270.0
        offset = 3*steps_per_rev / 4 + 1;
    otherwise
        return
end

%Plot data for 7 m/s, 90 deg azimuth, and various yaw angles
% step 5 of revolution

for j=1:numWindSpeeds
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 1 scrsz(3) scrsz(4)])
    
   % text(.5,.5, );

    for i=1:numYaw
        sp1 = subplot(2,3,i);
       % title(sp1,['Azimuth = 90.0 deg, Uinf= ' num2str(windSpeeds(j)) ] );
        plot(r_over_R_FAST, WTP_Fx(steps_per_rev*(i-1) + (numYaw*steps_per_rev*(j-1)) + offset, 1:nFASTNodes),'r','LineWidth',1.5);
        hold on;
        plot(r_over_R_FAST,-WTP_Fy(steps_per_rev*(i-1) + (numYaw*steps_per_rev*(j-1)) + offset, 1:nFASTNodes),'k','LineWidth',1.5);

        plot(r_over_R,dataCCblade(:,1,numYaw*(j-1)+i)','r--','LineWidth',1.5); 
        plot(r_over_R,dataCCblade(:,2,numYaw*(j-1)+i)','k--','LineWidth',1.5); 

        plot(r_over_R_FAST,FAST_Fx(j+(i-1)*numWindSpeeds,:),'g-.','LineWidth',1.5);
        
        plot(r_over_R_FAST,FAST_Fy(j+(i-1)*numWindSpeeds,:),'b-.','LineWidth',1.5);

        plot(r_over_R_EXP(j+(i-1)*numWindSpeeds,:),EXP_Fx(j+(i-1)*numWindSpeeds,:),'ro','MarkerFaceColor','r','MarkerSize',5);
        plot(r_over_R_EXP(j+(i-1)*numWindSpeeds,:),EXP_Fy(j+(i-1)*numWindSpeeds,:),'ko','MarkerFaceColor','k','MarkerSize',5);

        xlabel('radius fraction');
        ylabel('force per unit length (N/m)');

        if i==3
            legend('BEMT-flapwise','BEMT-lead-lag','CCblade-flapwise','CCblade-lead-lag','FAST-flapwise','FAST-lead-lag','UAE-flapwise','UAE-lead-lag','location','NorthEastOutside');
        end;
        ylim([-100, 300]);
        title(['Yaw = ' num2str(yawAngles(i))]);

    end
    [ax1,h1]=suplabel('Pitt Peters correction: ON');
    [ax4,h3]=suplabel(['Azimuth = ' num2str(azimuth) ' deg, Uinf= ' num2str(windSpeeds(j)) ]  ,'t');
    set(h3,'FontSize',30)
    
end

%Plot data for 7 m/s, 270 deg azimuth, and various yaw angles
% step 13 of revolution




%plot(r_over_R(3:end),WTP_AOA(1,3:numNodes));
% figure;
% plot(r_over_R(),WTP_Fx(1,1:numNodes));
% hold on;
% plot(r_over_R(),-WTP_Fy(1,1:numNodes));


% Chans = strmatch('Fx',ChanName);
% Fx = Channels(1:steps_per_rev,Chans);
% 
% Chans = strmatch('Fy',ChanName);
% Fy = Channels(1:steps_per_rev,Chans);


%plot(r_over_R(3:end),alphaData(3:end))
%plot(r_over_R(1:2),alphaData(1:2))

% Now read the BEM module results and plot with the FAST8 results


%hold on;
%plot(r_over_R(3:21),alphaDataBEM(3:21))
%plot(r_over_R(1:2),alphaDataBEM(1:2))

%axis tight;
% xlabel('nameX');
% ylabel('nameY');
% legend('name1','name2')
% title('asdfasdfasd')
% gtext('asdfasdf')   

end

function [ r_over_R_exp, EXP_Fx, EXP_Fy ] = Read_Experimental_Results(wind,yaw, azimuth)

    basename = 'C:\Dev\aerodyn_skew\results\';
    nYaw = length(yaw);
    nWind = length(wind);
    for j=1:nWind
        for i=1:nYaw
            file = [num2str(wind(j)) '_Y' num2str(yaw(i),'%3.1f') '_A' num2str(azimuth) '_'];
            name = [basename file 'r_exp.mat'];
            load(name);
            r_over_R_exp(j+(i-1)*nWind,:) = r_exp';
            name = [basename file 'tq.mat'];
            load(name);
            EXP_Fy(j+(i-1)*nWind,:) = tq_mean';
            name = [basename file 'th.mat'];
            load(name);
            EXP_Fx(j+(i-1)*nWind,:) = th_mean';
        end
    end


end


function [ FAST_AOA, FAST_Fx, FAST_Fy ] = Read_FAST_Results(azimuth)

    steps_per_rev  = 160;
    
    HeaderRows     = 3;
    NameLine       = 2;
    UnitsLine      = 3;
    delim          = '';
    switch azimuth
        case 0.0
            offset = 0;
        case 90.0
            offset =   steps_per_rev / 4 - 1;
        case 180.0
            offset =   steps_per_rev / 2 - 1;
        case 270.0
            offset = 3*steps_per_rev / 4 - 1;
        otherwise
            return
    end
    five_mps   =  4*steps_per_rev + offset;
    seven_mps  = 14*steps_per_rev + offset;
    ten_mps    = 24*steps_per_rev + offset;
    keyRows = [five_mps, seven_mps, ten_mps];
    dr=[0.1961,0.3520,0.3520,0.2012,0.2012,0.2347,0.2012,0.2012,0.2012,0.2349,0.2010,0.2012,0.2012,0.2347,0.2012,0.2012,0.2012,0.1509,0.2012,0.2012,0.1509];
    drmat = repmat(dr,18,1);
    
    % These are the Yaw = 0.0 results for the 3 windspeeds
    file = 'C:\Dev\Fastv8.08 download\BEM Verification\UAE_VI\Test01\UAE_Test1.AD.out';
    [FASTChannels, FASTChanName, ~, ~] = ReadFASTtext(file, delim, HeaderRows, NameLine, UnitsLine );
    
    % Find all the AOA, Fx, Fy channels
    
    chans = strmatch('Alpha',FASTChanName);
    FAST_AOA = FASTChannels(keyRows,chans);
    chans = strmatch('ForcN',FASTChanName);
    FAST_Fx = FASTChannels(keyRows,chans);
    chans = strmatch('ForcT',FASTChanName);
    FAST_Fy = FASTChannels(keyRows,chans);
    
    % These are the Yaw = 20.0 results for the 3 windspeeds
    file = 'C:\Dev\Fastv8.08 download\BEM Verification\UAE_VI\Test02\UAE_Test2.AD.out';
    [FASTChannels, FASTChanName, ~, ~] = ReadFASTtext(file, delim, HeaderRows, NameLine, UnitsLine );
    
    chans = strmatch('Alpha',FASTChanName);
    FAST_AOA(4:6,:) = FASTChannels(keyRows,chans);
    chans = strmatch('ForcN',FASTChanName);
    FAST_Fx(4:6,:) = FASTChannels(keyRows,chans);
    chans = strmatch('ForcT',FASTChanName);
    FAST_Fy(4:6,:) = FASTChannels(keyRows,chans);
    
    
    % These are the Yaw = 30.0 results for the 3 windspeeds
    file = 'C:\Dev\Fastv8.08 download\BEM Verification\UAE_VI\Test03\UAE_Test3.AD.out';
    [FASTChannels, FASTChanName, ~, ~] = ReadFASTtext(file, delim, HeaderRows, NameLine, UnitsLine );
    
    chans = strmatch('Alpha',FASTChanName);
    FAST_AOA(7:9,:) = FASTChannels(keyRows,chans);
    chans = strmatch('ForcN',FASTChanName);
    FAST_Fx(7:9,:) = FASTChannels(keyRows,chans);
    chans = strmatch('ForcT',FASTChanName);
    FAST_Fy(7:9,:) = FASTChannels(keyRows,chans);
    
    
    % These are the Yaw = 60.0 results for the 3 windspeeds
    file = 'C:\Dev\Fastv8.08 download\BEM Verification\UAE_VI\Test04\UAE_Test4.AD.out';
    [FASTChannels, FASTChanName, ~, ~] = ReadFASTtext(file, delim, HeaderRows, NameLine, UnitsLine );
    
    chans = strmatch('Alpha',FASTChanName);
    FAST_AOA(10:12,:) = FASTChannels(keyRows,chans);
    chans = strmatch('ForcN',FASTChanName);
    FAST_Fx(10:12,:) = FASTChannels(keyRows,chans);
    chans = strmatch('ForcT',FASTChanName);
    FAST_Fy(10:12,:) = FASTChannels(keyRows,chans);
    
    
    % These are the Yaw = 75.0 results for the 3 windspeeds
    file = 'C:\Dev\Fastv8.08 download\BEM Verification\UAE_VI\Test05\UAE_Test5.AD.out';
    [FASTChannels, FASTChanName, ~, ~] = ReadFASTtext(file, delim, HeaderRows, NameLine, UnitsLine );
    
    chans = strmatch('Alpha',FASTChanName);
    FAST_AOA(13:15,:) = FASTChannels(keyRows,chans);
    chans = strmatch('ForcN',FASTChanName);
    FAST_Fx(13:15,:) = FASTChannels(keyRows,chans);
    chans = strmatch('ForcT',FASTChanName);
    FAST_Fy(13:15,:) = FASTChannels(keyRows,chans);
    
    % These are the Yaw = 90.0 results for the 3 windspeeds
    file = 'C:\Dev\Fastv8.08 download\BEM Verification\UAE_VI\Test06\UAE_Test6.AD.out';
    [FASTChannels, FASTChanName, ~, ~] = ReadFASTtext(file, delim, HeaderRows, NameLine, UnitsLine );
    
    chans = strmatch('Alpha',FASTChanName);
    FAST_AOA(16:18,:) = FASTChannels(keyRows,chans);
    chans = strmatch('ForcN',FASTChanName);
    FAST_Fx(16:18,:) = FASTChannels(keyRows,chans);
    chans = strmatch('ForcT',FASTChanName);
    FAST_Fy(16:18,:) = FASTChannels(keyRows,chans);
    
    
    FAST_Fx = FAST_Fx./drmat;
    FAST_Fy = FAST_Fy./drmat;
end

function [ WTP_AOA, WTP_Fx, WTP_Fy ] = Read_WTP_Results(WTP_Results_File)


% Read the WT_Perf data first
HeaderRows    = 5;
NameLine      = 4;
UnitsLine     = 5;
delim         = '';

[WTPChannels, WTPChanName, ~, ~] = ReadFASTtext(WTP_Results_File, delim, HeaderRows, NameLine, UnitsLine );

WTP_AOA = ExtractData(WTPChanName, WTPChannels, 'AOA');
WTP_Fx  = ExtractData(WTPChanName, WTPChannels, 'Fx');
WTP_Fy  = ExtractData(WTPChanName, WTPChannels, 'Fy');


end

function data = Get_ccblade_data(file, NumIC)


    fid = fopen( file, 'rt' );

      if ( fid < 0 )
         beep
         error( '    Could not open "%s" for reading.', file );
      end



         % Read the numeric data and store it in the time-series matrix.

      temp = textscan(fid,repmat('%f ',1,NumIC),'CollectOutput',1);
      data = temp{1}';
      

         % Close the data file.

      fclose ( fid );
end

function dataCCblade = Read_CCblade_Results(numNodes,numYaw,numWindSpeeds,azimuth)

dataCCblade = zeros(numNodes,2,numYaw*numWindSpeeds);

file = ['C:\Dev\aerodyn_skew\results\5_Y0.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,1) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\7_Y0.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw+1) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\10_Y0.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw*2+1) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\5_Y20.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,2) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\7_Y20.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw+2) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\10_Y20.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw*2+2) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\5_Y30.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,3) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\7_Y30.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw+3) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\10_Y30.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw*2+3) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\5_Y60.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,4) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\7_Y60.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw+4) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\10_Y60.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw*2+4) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\5_Y75.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,5) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\7_Y75.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw+5) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\10_Y75.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw*2+5) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\5_Y90.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,6) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\7_Y90.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw+6) = Get_ccblade_data(file, numNodes);

file = ['C:\Dev\aerodyn_skew\results\10_Y90.0_A' num2str(azimuth,'%3.0f') '.txt'];
dataCCblade(:,:,numYaw*2+6) = Get_ccblade_data(file, numNodes);

end

function data = ExtractData(InputNames,InputData, str)

    k = strfind(InputNames,str);
    count = 1;
    
    for i=1:length(InputNames)
        if ~isempty(k{i})
            found(count) = i;
            count = count+1;
        end
    end

    data = InputData(1:end,found);
    
end
