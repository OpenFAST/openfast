function [YawManRat, PitManRat] = CalculateYawAndPitchRates(inputFile, outputFile)
%function [YawManRat, PitManRat] = CalculateYawAndPitchRates(inputFile, outputFile)
%
% This routine takes FAST pre-v8 input/output files and converts the 
% pitch and yaw maneuver end-times to equivalent rates to be specified in
% the FAST v8.0 input files.
%
%Inputs:
%   inputFile  - the primary FAST input file (pre-v8.*)
%   outputFile - the FAST output file containing the yaw position and/or
%                blade pitch outputs necessary to perform the rate calculations.
%Outputs:
% YawManRat - the equivalent yaw maneuver rate
% PitManRat - the equivalent pitch maneuver rates for all 3 blades (0 if
%             not defined)
%..........................................................................

% FAST pre-v8.* input file
FP = Fast2Matlab(inputfile,4); %FP are Fast Parameters, specify 4 lines of header

% output file
if isempty( strfind( lower(outputfile), '.outb') )
    [ChannelData, ChannelNames] = ReadFASTtext(outputfile);
else
    [ChannelData, ChannelNames] = ReadFASTbinary(FileName);
end

%% ........................................................................
% Calculate yaw maneuver rate:
%..........................................................................

% yaw position output names:
YawPosNames = {'YawPzn', 'YawPzp', 'NacYawP', 'NacYaw', 'YawPos'};

%find the column with the yaw maneuver output:
YawColumn = -1;
for i=1:length(YawPosNames)
    indx = strcmpi(ChannelNames,YawPosNames{i});
    if any( indx )
        YawColumn = find(indx,1,'first');
        break;
    end
end

if YawColumn <= 0
    disp( ['Error: could not find yaw position column in output file, ' outputfile] );
    YawManRat = 0;
else
    NacYawF  = GetFastPar(FP,'NacYawF' );
    TYawManE = GetFastPar(FP,'TYawManE');
    TYawManS = GetFastPar(FP,'TYawManS');
    
        % find the YawPos at time TYawManS
    timeIndx = find( ChannelData(:,1) >= TYawManS, 1, 'first' );
    YawPos   = ChannelData(timeIndx,YawColumn);
    
    YawManRat = abs( (NacYawF - YawPos) / (TYawManE - TYawManS) );
    fprintf( 'YawManRat = %f\n', YawManRat );
end

    
%% ........................................................................
% Calculate pitch maneuver rate:
%..........................................................................
NumBl  = GetFastPar(FP,'NumBl' );
PitManRat = zeros(NumBl,1);

% pitch position output names:
BlPitchNames = {'PtchPMzc', 'PtchPMzb', 'BldPitch', 'BlPitch'};

for k=1:NumBl
    blStr = num2str(k);
    
    %find the column with the yaw maneuver output:
    PitchColumn = -1;
    for i=1:length(BlPitchNames)
        indx = strcmpi(ChannelNames,[BlPitchNames{i} blStr]);
        if any( indx )
            PitchColumn = find(indx,1,'first');
            break;
        end
    end

    if PitchColumn <= 0
        disp( ['Error: could not find blade pitch column for blade ' blStr ' in output file, ' outputfile] )
    else
        BlPitchF = GetFastPar(FP,['BlPitchF(' blStr ')'] );
        TPitManE = GetFastPar(FP,['TPitManE(' blStr ')'] );
        TPitManS = GetFastPar(FP,['TPitManS(' blStr ')'] );

            % find the BlPitch at time TPitManS
        timeIndx = find( ChannelData(:,1) >= TPitManS, 1, 'first' );
        BlPitch  = ChannelData(timeIndx,PitchColumn);

        PitManRat(k) = abs( (BlPitchF - BlPitch) / (TPitManE - TPitManS) );
        fprintf( 'PitManRat(%1.0f) = %f\n', k, PitManRat(k) );
    end
end

return


