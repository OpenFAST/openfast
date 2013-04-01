%% Sweep Plot
%Run this script directly following LinSweep

clear freqMat freqMatNoZero %get rid of sizing
%empty MBC3 variables
clear AC ADP AK AMat AvgAMat AvgBMat AvgBdMat  AvgCMat  AvgDMat  AvgDdMat  Avgxdop  Avgxop  Avgyop  AzimB1Up  Azimuth  BMat  BdMat  CMat  DMat  DdMat  DescCntrlInpt  DescDisturbnc  DescOutput  DescStates  FP  MBC_A  MBC_AvgA  MBC_B  MBC_Bd  MBC_C  MBC_D  MBC_DampedFrequency  MBC_DampedFrequencyHz  MBC_DampingRatio  MBC_Dd  MBC_DecrementRate  MBC_EigenVects  MBC_Evals  MBC_ModeShapeMagnitude  MBC_ModeShapePhaseDeg  MBC_NaturalFrequency  MBC_NaturalFrequencyHz  MBC_eigenVals  MBC_eigenVects  MdlOrder  N  NActvDOF NAzSteps  NAzimStep  NDisturbs  NInputs  NRotTripletCntrlInpt  NRotTripletOutput  NRotTripletStates  NumOuts  Omega  OmegaDot  OmgDot  Period  ProgName  RPMSetPoint RotTripletIndicesCntrlInpt  RotTripletIndicesOutput  RotTripletIndicesStates  T1  T1c  T1o  T1ov  T1v  T2  T3  az1  az2  az3  azm a



%% First process initial windspeed

%First, grab the first data set and store the eigenvectors
RootName = linList{1}(1:end-4); %define the first RootName

%Load linearization and analyze
GetMats;
MBC3;

%Grab frequencies
freq =  MBC_DampedFrequencyHz;

% %Grab eigenvectors
% initialEig = MBC_ModeShapeMagnitude;
% %Sort frequencies and eigenvectors
% [junk idx] = sort(freq);
% freq = freq(idx);
% initalEig = initialEig(:,idx);

%Now place in first column of matrix
if windSweep(1) == 0 %catch missing modes in stopped case
    freqMat(:,1) = [0;0;freq]; %pad in missing rigid body modes
else
    freqMat(:,1) = freq;
end


%% Loop through remaining wind speeds
for i = 2:length(linList)
    
    %clear MBC vars
    clear AC ADP AK AMat AvgAMat AvgBMat AvgBdMat  AvgCMat  AvgDMat  AvgDdMat  Avgxdop  Avgxop  Avgyop  AzimB1Up  Azimuth  BMat  BdMat  CMat  DMat  DdMat  DescCntrlInpt  DescDisturbnc  DescOutput  DescStates  FP  MBC_A  MBC_AvgA  MBC_B  MBC_Bd  MBC_C  MBC_D  MBC_DampedFrequency  MBC_DampedFrequencyHz  MBC_DampingRatio  MBC_Dd  MBC_DecrementRate  MBC_EigenVects  MBC_Evals  MBC_ModeShapeMagnitude  MBC_ModeShapePhaseDeg  MBC_NaturalFrequency  MBC_NaturalFrequencyHz  MBC_eigenVals  MBC_eigenVects  MdlOrder  N  NActvDOF NAzSteps  NAzimStep  NDisturbs  NInputs  NRotTripletCntrlInpt  NRotTripletOutput  NRotTripletStates  NumOuts  Omega  OmegaDot  OmgDot  Period  ProgName  RPMSetPoint RotTripletIndicesCntrlInpt  RotTripletIndicesOutput  RotTripletIndicesStates  T1  T1c  T1o  T1ov  T1v  T2  T3  az1  az2  az3  azm a
    
    RootName = linList{i}(1:end-4); %define the current RootName
    
    %Load linearization and analyze
    GetMats;
    MBC3;
    
    
    %grab frequencies for this iteration
    freq = MBC_DampedFrequencyHz;
    
    %grab eigenvectors for this iteration
    curEig = MBC_ModeShapeMagnitude;
    
%     %sort by MAC comparison
%     %loop through current eigenvector matrix, find best match
%     for j = 1:size(curEig,2)
%         testEig = curEig(:,j);
%         for k = 1:size(initialEig,2)
%             testInitEig = initalEig(:,k);
%             testResult(k) = ((testEig' * testInitEig)^2)/( (testEig' * testEig) * (testInitEig' * testInitEig) );
%         end
%         
%         [junk maxTestIdx] = max(testResult);
%         idx(j) = maxTestIdx;
%     end
%     freq = freq(idx);
    
    %Assign into matrix
    freqMat(:,i) = freq;

end

%% Remove all 0 frequencies, these aren't too helpful
for i = 1:size(freqMat,2)
    freqCol = freqMat(:,i);
    idx = find(freqCol > 0.5)
    freqMatNoZero(:,i) = freqCol(idx);
end


%% Plot
figure
hold all
title('Modal Frequencies vs Windspeed')
ylabel('Frequency (Hz)')
xlabel('Windspeed (m/s)')
set(gca,'LineStyleOrder', '-d|-*|-v|-+|-o|-.|-x|-s|-^|-h')
styleVector = {'b-d','k-*','r-v','b-+','k-o','r.-','b-x','k-s','r-^','b-h','m-d','g-*','m-v','g-+','m-o','g.-','m-x','g-s','m-^','g-h'};
%set(gca,'ColorOrder','k|b')
for i = 1:size(freqMatNoZero,1)
    plot(windSweep,squeeze(freqMatNoZero(i,:)),styleVector{i})
    
end

hold off