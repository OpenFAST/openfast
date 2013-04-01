% Matlab2FAST
% Function for creating a new FAST file given:
% 1) An input template file
% 2) A FAST parameter structure
%
% In:   FastPar         -   A FAST parameter list
%       TemplateFile    -   A .fst file to use as a template
%       OutputFilename  -   Desired filename of output .fst file
%
% Paul Fleming, JUNE 2011
% Using code copied from functions written by Knud Kragh
% Modified by Bonnie Jonkman, Feb 2013
% Note: the template must not contain the string "OutList" anywhere
%  except for the actual "OutList" output parameter lists.

function Matlab2FAST(FastPar,TemplateFile,OutputFilename,hdrLines)

if nargin < 4
    hdrLines = 0;
end

% we're going to get the appropriate newline character(s) from the template file

HaveNewLineChar = false;
% newline = '\r'; %mac
% newline = '\n'; %linux
% newline = '\r\n'; %windows
ContainsOutList = false;


%Declare file IDs from the template to the resulting file
fidIN = fopen(TemplateFile,'r');
if fidIN == -1
    error([' Template file, ' TemplateFile ', not found.'])
end
fidOUT = fopen(OutputFilename,'w');
if fidOUT == -1
    error(['Can''t create file ' OutputFilename])
end

for hi = 1:hdrLines
    line    = fgets(fidIN); 
    
        % get the line feed characters(s) here: CHAR(13)=CR(\r) CHAR(10)=LF(\n)
        % so our file has consistant line endings
    if ~HaveNewLineChar
        HaveNewLineChar = true;
        indx = min( strfind(line,char(13)), strfind(line,char(10)) );
        if ~isempty(indx)
            newline = line(indx:end);
        end
    end

        % print appropriate header lines:
    if isfield(FastPar, 'HdrLines') && hi ~= 1
        fprintf(fidOUT,'%s',FastPar.HdrLines{hi,1}); %does not contain the line ending
        fprintf(fidOUT,newline);                     %so print it here instead
    else
        fprintf(fidOUT,'%s',line);
    end
end

printTable = false; %assume we'll get the tables from the FastPar data structure;

%loop through the template up until OUTLIST or end of file
while true
    
    line = fgets(fidIN); %get the next line from the template
    
    if isnumeric(line) %we reached the end of the file
        break;
    end
    
        % get the line feed characters(s) here: CHAR(13)=CR(\r) CHAR(10)=LF(\n)
        % so our file has consistant line endings
    if ~HaveNewLineChar
        HaveNewLineChar = true;
        indx = min( strfind(line,char(13)), strfind(line,char(10)) );
        if ~isempty(indx)
            newline = line(indx:end);
        end
    end
    
    if ~isempty(strfind(upper(line),upper('OutList'))) 
        ContainsOutList = true;
        fprintf(fidOUT,'%s',line); %if we've found OutList, write the line and break 
        break;
    end  
    
    [value, label, isComment, ~, ~] = ParseFASTInputLine(line);
            
    if ~isComment && length(label) > 0        
        
        if strcmpi(value,'"HtFract"') %we've reached the distributed tower properties table (and we think it's a string value so it's in quotes)            
            if ~isfield(FastPar,'TowProp')
                disp(  ['WARNING: tower properties table not found in the FAST data structure.'] );
                printTable = true;
            else
                WriteFASTTable(line, fidIN, fidOUT, FastPar.TowProp, FastPar.TowPropHdr, newline);
                continue; %let's continue reading the template file            
            end

        elseif strcmpi(value,'"BlFract"') %we've reached the distributed blade properties table (and we think it's a string value so it's in quotes)
            if ~isfield(FastPar,'BldProp')
                disp( 'WARNING: blade properties table not found in the FAST data structure.' );
                printTable = true;
            else
                WriteFASTTable(line, fidIN, fidOUT, FastPar.BldProp, FastPar.BldPropHdr, newline);
                continue; %let's continue reading the template file            
            end
        else

            indx = strcmpi( FastPar.Label, label );
            if any( indx )

                if sum(indx) > 1 % we found more than one....
                    disp( ['WARNING: multiple occurrences of ' label ' in the FAST data structure.'] );
                end

                % The template label matches a label in FastPar
                %  so let's use the old value.
                indx2 = find(indx,1,'first');       
                val2Write = FastPar.Val{indx2}; 

                    % print the old value at the start of the line,
                    % using an appropriate format
                if isnumeric(val2Write)
                    writeVal = sprintf('%11G',val2Write(1));
                    if ~isscalar(val2Write) %Check for the special case of an array
                        writeVal = [writeVal sprintf(',%11G',val2Write(2:end)) ' '];
                    end
                else
                    writeVal = [val2Write repmat(' ',1,max(1,11-length(val2Write)))];
                end


                idx = strfind( line, label ); %let's just take the line starting where the label is first listed            
                line = [writeVal '   ' line(idx(1):end)];            

            else
                disp( ['WARNING: ' label ' not found in the FAST data structure. Default value listed below (from template file, ' TemplateFile ') will be used instead:'] )
                disp( value );
                disp( '' );            
            end

        end
    else % isComment || length(label) == 0 
        if isComment
            printTable = false;     % we aren't reading a table (if we were, we reached the end) 
        else
            if ~printTable
                continue;           % don't print this line without a label
            end
        end
    end       
    
    %write this line into the output file
    fprintf(fidOUT,'%s',line);
end

if ContainsOutList
    if isfield(FastPar,'OutList')
        OutListChar = char(FastPar.OutList);  %this will line up the comments nicer
        spaces      = repmat(' ',1,max(1,28-size(OutListChar,2)));
        %Now add the Outlist
        for io = 1:length(FastPar.OutList)
            fprintf(fidOUT,[OutListChar(io,:) spaces FastPar.OutListComments{io} newline]);
        end
    else
        disp( 'WARNING: OutList was not found in the FAST data structure. The OutList field will be empty.' );        
    end

        %Now add the close of file
    fprintf(fidOUT,'END of input file (the word "END" must appear in the first 3 columns of this last OutList line)');
    fprintf(fidOUT,newline);
    fprintf(fidOUT,'---------------------------------------------------------------------------------------');
    fprintf(fidOUT,newline);
end

fclose(fidIN);
fclose(fidOUT);
end %end function

function WriteFASTTable( HdrLine, fidIN, fidOUT, Table, Headers, newline )

    % we've read the line of the template table that includes the header 
    % let's parse it now:
    TmpHdr = textscan(HdrLine,'%s');
    TemplateHeaders = TmpHdr{1};
    nc = length(TemplateHeaders);

    fprintf(fidOUT,'%s',HdrLine);           % print the new headers
    fprintf(fidOUT,'%s',fgets(fidIN));      % print the new units (we're assuming they are the same)
    
    % let's figure out which columns in the old Table match the headers
    % in the new table:
    ColIndx = ones(1,nc);
    

    for i=1:nc
        indx = strcmpi(TemplateHeaders{i}, Headers);
        if sum(indx) > 0
            ColIndx(i) = find(indx,1,'first');
            if sum(indx) ~= 1
                disp( ['WARNING: Multiple instances of ' TemplateHeaders{i} ' column found in FAST table.'] );
            end
        else
           error( [ TemplateHeaders{i} ' column not found in FAST table. Cannot write the table.'] );
        end                
    end
    
    
    % now we'll write the table:
    for i=1:size(Table,1) 
        fprintf(fidOUT, '%11.7E  ', Table(i,ColIndx) );  %write all of the columns
        fprintf(fidOUT, newline);
    end
              
end

