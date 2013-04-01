% Function for getting tower data from FAST tower file
%
% In:   -   Name of FAST tower file
% Out:  -   Data structure containing all structural blade data
%
% Note: Not a very general function, take advantages of blade file format,
% check format of blade file in example. Data section must contain numbers
% in all columns, also those anly needed for ADAMS (just put zeros if you
% dont have the additional data)
%
% Knud Abildgaard Kragh


function DataOut=Tower2Matlab(TowerFile)

fid = fopen([TowerFile],'r');
if fid == -1
    Flag = 0;
    disp('  ')
    disp('==============================================================')
    disp(['Tower file could not be found'])
    disp('--------------------------------------------------------------')
    return
end
for i=1:4
    tline = fgets(fid); % Line n
end

DataOut.nSec=fscanf(fid,'%f',1);

for i=1:11
    tline = fgets(fid); % Line n
end

DataOut.MassTune=fscanf(fid,'%f',1);
tline = fgets(fid); % Line n
DataOut.FAStiffTune=fscanf(fid,'%f',1);
tline = fgets(fid); % Line n
DataOut.SSStiffTune=fscanf(fid,'%f',1);
tline = fgets(fid); % Line n
tline = fgets(fid); % Line n

for i=1:10
    DataOut.Labels{i}=fscanf(fid,'%s',1);
end

tline = fgets(fid); % Line n
tline = fgets(fid); % Line n

for i=1:DataOut.nSec
    for ii=1:10
        DataOut.Data(i,ii)=fscanf(fid,'%f',1);
    end
    tline = fgets(fid); % Line n
end

fclose(fid);


end