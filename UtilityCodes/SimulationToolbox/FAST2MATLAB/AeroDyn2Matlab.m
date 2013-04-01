% AeroDyn2Matlab
% Function for reading AD input files in to a MATLAB struct.
%
%
%This function returns a structure DataOut, which contains 2 cell arrays,
%.Val: An array of values
%.Label: An array of matching labels
%.FoilNm: A Cell array of foil names
%.BldNodes A matrix of blade nodes with columns RNodes, AeroTwst DRNodes
%Chord and Nfoil
%.PrnElm An array determining whether or not to print a given element
%These arrays are extracted from the FAST input file
%
% In:   AD_file    -   Name of AD input file
%       hdrLines    -   Number of lines to skip at the top (optional)
%
% Knud A. Kragh, May 2011, NREL, Boulder
%
% Modfied by Paul Fleming, JUNE 2011

function DataOut = AeroDyn2Matlab(AD_file,hdrLines)

if nargin < 2
    hdrLines = 0;
end

%----------------------Read FST main file----------------------------------
fid = fopen([AD_file],'r');
if fid == -1
    error('AD file could not be found')
end

%skip hdr
for hi = 1:hdrLines
    DataOut.HdrLines{hi,1} = fgetl(fid); %bjj: added field for storing header lines
end

%PF: Commenting this out, not sure it's necessary
%DataOut.Sections=0;

%Loop through the file line by line, looking for value-label pairs
%Stop once we've reached the FoilNm, meaning we've reached the list of foil
%names
count=1;


while true %loop until discovering FoilNm, than break
    
    line = fgetl(fid);
    
    [DataOut.Val{count,1}, DataOut.Label{count,1}] = ParseFASTInputLine( line );
    
    %check if the last label read is the NumFoil
    if strcmpi(DataOut.Label{count},'NumFoil')
        numFoil = DataOut.Val{count};
        break; %if it does end the loop
    end
    
    %if not NumFoil iterate
    count=count+1;
    
    
    % end %endif %PF no skipping in AD
end %end while


%% Read in FoilNm list
% numFoil = GetFastPar(DataOut,'NumFoil');
for i = 1:numFoil
    line = fgetl(fid);
    DataOut.FoilNm{i,1} = ParseFASTInputLine( line );
end

%% Read in BldNodes and PrnElm

%Get the number of blade nodes
line = fgetl(fid);
count = count + 1;
[DataOut.Val{count}, DataOut.Label{count}] = ParseFASTInputLine( line );
BldNodes = DataOut.Val{count};

%skip the header row of table
fgets(fid);

%Now loop through and get all the data
for i = 1:BldNodes
    fgets(fid); %Advance to the next line
    for col = 1:5
        DataOut.BldNodes(i,col) = fscanf(fid,'%f',1);
    end
    DataOut.PrnElm{i} = fscanf(fid,'%s',1);
end

fclose(fid); %close file

end %end function