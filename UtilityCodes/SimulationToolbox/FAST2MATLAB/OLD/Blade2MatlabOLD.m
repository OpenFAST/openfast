% Function for getting blade data from FAST blade file
%
% In:   -   Name of FAST blade file
% Out:  -   Data structure containing all structural blade data
%
% Note: Not a very general function, take advantages of blade file format,
% check format of blade file in example. Data section must contain numbers
% in all columns, also those anly needed for ADAMS (just put zeros if you
% dont have the additional data)
%
% Knud Abildgaard Kragh


function DataOut=Blade2Matlab(BladeFile)

fid = fopen([BladeFile],'r');
if fid == -1
    Flag = 0;
    disp('  ')
    disp('==============================================================')
    disp(['Blade file could not be found'])
    disp('--------------------------------------------------------------')
    return
end
for i=1:4
    tline = fgets(fid); % Line n
end

DataOut.nSec=fscanf(fid,'%f',1);

for i=1:8
    tline = fgets(fid); % Line n
end

DataOut.MassTune=fscanf(fid,'%f',1);
tline = fgets(fid); % Line n
DataOut.FlapStiffTune=fscanf(fid,'%f',1);
tline = fgets(fid); % Line n
DataOut.EdgeStiffTune=fscanf(fid,'%f',1);
tline = fgets(fid); % Line n
tline = fgets(fid); % Line n

for i=1:17
    DataOut.Labels{i}=fscanf(fid,'%s',1);
end
 
tline = fgets(fid); % Line n
tline = fgets(fid); % Line n

for i=1:DataOut.nSec
    for ii=1:17
        DataOut.Data(i,ii)=fscanf(fid,'%f',1);
    end
    tline = fgets(fid); % Line n
end

fclose(fid);


end