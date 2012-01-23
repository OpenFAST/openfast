function [t,cts,dat,lab] = readCTSfile(fileName)
%%
fid = fopen(fileName);

if fid < 1
    error(['readCTSfile::Error reading file "' fileName '"']);
end

nh = 9;

typ = textscan(fid,'%s%*s%s',1);

d=textscan(fid,'%f%*s%s',nh-1);
% dat = d{1};
% lab = d{2}; 

if strcmpi(typ{1}{1},'les')
    dat = [1; d{1}];
else
    dat = [0; d{1}];
end
lab = {typ{2}{1}, d{2}{:}};

[tmp] = textscan(fid,'%f%f',d{1}(end));

t   = tmp{1};
cts = tmp{2};

fclose(fid);


end

