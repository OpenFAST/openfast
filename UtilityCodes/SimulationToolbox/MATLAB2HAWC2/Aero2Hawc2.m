% Fucntion for writing HAWC2 AE and PC files from Aerodyne data extracted
% to matlab
%
% In:   OutPath         -       Path to folder holding hawc2 files
%       OutFileNames    -       Names of output files .ae and .pc
%       AeroDyn_AE      -       Aerodyne data extracted using
%       AeroDyn2Matlab.m
%       PC_Data         -       Profile coefficient extracted using
%       Fast_PC2Matlab.m
%       FastPar         -       Extracted using Fast2Matlab.m
%
% Knud A. Kragh

function Aero2Hawc2(OutPath,OutFileNames,Hawc2_Ae,AeroDyn_AE,PC_Data,FastPar)

mkdir(OutPath,'Data');
% Writing ae file
fid = fopen([OutPath 'Data\' OutFileNames '.ae' ], 'w+');
fprintf(fid,'1        Chord[m]  T/C  Set no. OBS: T/C number are not real!! \n');
fprintf(fid,['1 ' num2str(size(Hawc2_Ae,1)) ' Blade \n']);

for i=1:size(Hawc2_Ae,1)
    fprintf(fid,[num2str(Hawc2_Ae(i,1)) '\t' num2str(Hawc2_Ae(i,2)) '\t' num2str(Hawc2_Ae(i,3)) '\t' '1' '\n']);       
end
fclose(fid);

% Writing PC file
PC_active=unique(AeroDyn_AE.BldNodes(:,end));

fid = fopen([OutPath 'Data\' OutFileNames '.pc' ], 'w+');
fprintf(fid,'1 Blade aerodynamic coefficients for converted from FAST OBS: Average values   \n');
fprintf(fid,[ num2str(length(PC_active))  '\n']);

for i=1:length(PC_active)
    fprintf(fid,[num2str(i) '  ' num2str(size(PC_Data(PC_active(i)).DataAv,1)) ' ' num2str(round(PC_active(i))) '  '  '-------------------------------------------------------------\n']);
    
    for j=1:size(PC_Data(PC_active(i)).DataAv,1)
        fprintf(fid,[num2str(PC_Data(PC_active(i)).DataAv(j,1)) '\t' num2str(PC_Data(PC_active(i)).DataAv(j,2)) '\t' num2str(PC_Data(PC_active(i)).DataAv(j,3)) '\t' num2str(0.0) '\n']);
    end
end



fclose(fid);
