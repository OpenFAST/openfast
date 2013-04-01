% Matlab2FASTLin
% Function for creating a new FAST Linear file given:
% 1) An input template file
% 2) A FAST Linear parameter structure
%
%In:    LinPar          -   A FAST parameter list
%       TemplateFile    -   A .fst file to use as a template
%       OutputFilename: Desired filename of output .fst file
%
% Paul Fleming, JUNE 2011
% Using code copied from functions written by Knud Kragh


function Matlab2FASTLin(LinPar,TemplateFile,OutputFilename)


%Declare file IDs from the template to the resulting file
fidIN = fopen(TemplateFile,'r');
if fidIN == -1
    error([TemplateFile ' not found'])
end
fidOUT = fopen(OutputFilename,'w');
if fidOUT == -1
    error(['Cant create ' OutputFilename])
end


%loop through the template up until OUTLIST
while true
    line = fgets(fidIN); %get the next line from the template
    
    if ~ischar(line) %if found the end of file
        break;  %EOF
    end
    
    %now see if there is a match to a known parameter
    idx = 9999; %9999 code is no match
    for i = 1:length(LinPar.Label)
        idxTemp = strfind(line,LinPar.Label{i});
        if ~isempty(idxTemp) %if there is a match, and it occurs before other matches
            idx = min(idx,idxTemp(1));
            val2Write = LinPar.Val{i};
        end
    end
    %if we found a match
    if (idx ~= 9999)
        if(line(1)~='-') %not a comment
            
            %build a string version of the value
            writeVal = num2str(val2Write);
            
            %             %Check for the special case of an array
            %             if isnumeric(val2Write)
            %                 if ~isscalar(val2Write)
            %                     writeVal = [];
            %                     for ii = 1:length(val2Write)
            %                         writeVal = [writeVal ' ' num2str(val2Write(ii))];
            %                     end
            %                     writeVal(1) = []; %get rid of leading comma
            %                 end
            %             end
            
            %add space padding to make it look nice
            if length(writeVal) < 11
                for i = 1:12 - length(writeVal)
                    writeVal = [writeVal ' '];
                end
            else
                writeVal = [writeVal '     '];
            end
            
            line = [writeVal line(idx(1):end)];
        end %endif
    end %endif
    
    
    %write this line into the output file
    fprintf(fidOUT,'%s',line);
end

fclose(fidIN);
fclose(fidOUT);
end %end function
