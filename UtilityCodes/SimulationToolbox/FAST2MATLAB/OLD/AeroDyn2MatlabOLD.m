% File for reading AeroDyn Input file
% In    -   Aerodyn input file name
% Out   -   Data structure with pc-data in Aerodyn format
%
% Knud A. Kragh

function DataOut=AeroDyn2Matlab(AED_file)

%----------------------Read FST main file----------------------------------
fid = fopen([AED_file],'r');
if fid == -1
    Flag = 0;
    disp('  ')
    disp('==============================================================')
    disp(['AED file could not be found'])
    disp('--------------------------------------------------------------')
    return
end

for i=1:17
    tline = fgets(fid); % Line n
end

DataOut.nFiles=fscanf(fid,'%f',1);
tline = fgets(fid); % Line n

for i=1:DataOut.nFiles
    DataOut.Files{i}=fscanf(fid,'%s',1);
    DataOut.Files{i}=DataOut.Files{i}(2:end-1);
    tline = fgets(fid); % Line n
end

DataOut.nAESec=fscanf(fid,'%f',1);

tline = fgets(fid); % Line n
tline = fgets(fid); % Line n
DataOut.Labels=(tline);

for i=1:DataOut.nAESec
    DataOut.AEData(i,:)=fscanf(fid,'%f',5)';
    tline = fgets(fid); % Line n
end

end
