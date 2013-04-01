%ImportFastOut
%Paul Fleming
%2/6/2013
%Quick tool for importing a FAST output file

function ImportFastOut(filename)

if nargin < 1
    listing = dir('*.out');
    filename = listing.name;
end
    
fprintf('Importing %s\n',filename);    
    

A = importdata(filename, '\t',  8);
channels = strtrim(regexp(A.textdata{7,1}, '\t', 'split'));
units = A.colheaders;
data = A.data;
assignin('base','channels',channels)
assignin('base','units',units)
assignin('base','data',data)