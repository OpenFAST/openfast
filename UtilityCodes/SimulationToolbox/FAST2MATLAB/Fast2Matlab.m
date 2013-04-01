function DataOut = Fast2Matlab(FST_file,hdrLines,DataOut)
%% Fast2Matlab
% Function for reading FAST input files in to a MATLAB struct.
%
%
%This function returns a structure DataOut, which contains the following 
% cell arrays:
%.Val              An array of values
%.Label            An array of matching labels
%.HdrLines         An array of the header lines (size specified at input)
%
% The following cell arrays may or may not be part of the DataOut structure
% (depending on the type of file being read):
%.OutList          An array of variables to output
%.OutListComments  An array of descriptions of the .OutList values
%.TowProp          A matrix of tower properties with columns .TowPropHdr 
%.TowPropHdr       A cell array of headers corresponding to the TowProp table
%.BldProp          A matrix of blade properties with columns .BldPropHdr
%.BldPropHdr       A cell array of headers corresponding to the BldProp table
%
%These arrays are extracted from the FAST input file
%
% In:   FST_file    -   Name of FAST input file
%       hdrLines    -   Number of lines to skip at the top (optional)
%
% Knud A. Kragh, May 2011, NREL, Boulder
%
% Modified by Paul Fleming, JUNE 2011
% Modified by Bonnie Jonkman, February 2013 (to allow us to read the 
% platform file, too)
%%
if nargin < 2
    hdrLines = 0;
end

%----------------------Read FST main file----------------------------------
fid = fopen(FST_file,'r');

if fid == -1
    error(['FST file, ' FST_file ', could not be found'])
end

%skip hdr
for hi = 1:hdrLines
    if nargin == 3
        fgetl(fid);
    else
        DataOut.HdrLines{hi,1} = fgetl(fid); %bjj: added field for storing header lines
    end
end

%PF: Commenting this out, not sure it's necessary
%DataOut.Sections=0;

%Loop through the file line by line, looking for value-label pairs
%Stop once we've reached the OutList which this function is the last
%occuring thing before the EOF
if nargin == 3
    count = max(length(DataOut.Label),length(DataOut.Var))+1;
else
    count = 1;
end


while true %loop until discovering Outlist or end of file, than break
    
    line = fgetl(fid);
    
    if isnumeric(line) % we reached the end of the file
        break
    end
    
        % Check to see if the value is Outlist
    if ~isempty(strfind(upper(line),upper('OutList'))) 
        [DataOut.OutList DataOut.OutListComments] = ParseFASTOutList(fid);
        break; %bjj: we could continue now if we wanted to assume OutList wasn't the end of the file...
    end         
        
    [value, label, isComment, descr, fieldType] = ParseFASTInputLine( line );    
    

    if ~isComment
        
        if strcmpi(value,'"HtFract"') %we've reached the distributed tower properties table (and we think it's a string value so it's in quotes)
            NTwInpSt = GetFastPar(DataOut,'NTwInpSt');        
            [DataOut.TowProp, DataOut.TowPropHdr] = ParseFASTTable(line, fid, NTwInpSt);
            continue; %let's continue reading the file
        elseif strcmpi(value,'"BlFract"') %we've reached the distributed blade properties table (and we think it's a string value so it's in quotes)
            NBlInpSt = GetFastPar(DataOut,'NBlInpSt');        
            [DataOut.BldProp, DataOut.BldPropHdr] = ParseFASTTable(line, fid, NBlInpSt);
            continue; %let's continue reading the file
        else                
            DataOut.Label{count,1} = label;
            DataOut.Val{count,1}   = value;
            count = count + 1;
        end
    end
    
end %end while

fclose(fid); %close file

return
end %end function
%%
function [OutList OutListComments] = ParseFASTOutList( fid )

    %Now loop and read in the OutList
    
    outCount = 0;
    while true
        line = fgetl(fid);
        [outVarC, position] = textscan(line,'%q',1); %we need to get the entire quoted line
        outVar  = outVarC{1}{1};    % this will not have quotes around it anymore...

        if isnumeric(line) %loop until we reach the word END or hit the end of the file
            break;
        else
            indx = strfind(upper(outVar),'END');
            if (~isempty(indx) && indx == 1) %we found "END" so that's the end of the file
                break;
            else
                outCount = outCount + 1;
                OutList{outCount,1} = ['"' outVar '"'];
                if position < length(line)
                  OutListComments{outCount,1} = line((position+1):end);
                else
                  OutListComments{outCount,1} = ' ';
                end
            end
        end
    end %end while   

    if outCount == 0
        disp( 'WARNING: no outputs found in OutList' );
        OutList = [];
        OutListComments = '';
    end
    
end %end function
%%
function [Table, Headers] = ParseFASTTable( line, fid, InpSt  )

    % we've read the line of the table that includes the header 
    % let's parse it now, getting the number of columns as well:
    TmpHdr  = textscan(line,'%s');
    Headers = TmpHdr{1};
    nc = length(Headers);

    % read the units line:
    fgetl(fid); 
        
    % now initialize Table and read its values from the file:
    Table = zeros(InpSt, nc);   %this is the size table we'll read
    i = 0;                      % this the line of the table we're reading    
    while i < InpSt
        
        line = fgetl(fid);
        if isnumeric(line)      % we reached the end prematurely
            break
        end        

        i = i + 1;
        Table(i,:) = sscanf(line,'%f',nc);       

    end
    
end %end function


