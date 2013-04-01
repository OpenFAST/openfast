% MATLAB2AD
% Function for creating a new AD file given:
% 1) An input template file
% 2) A FAST parameter structure
%
% In:   ADPar         -   An AD parameter list
%       TemplateFile    -   A .ipt file to use as a template
%       OutputFilename  -   Desired filename of output .ipt file
%
% Paul Fleming, JULY 2011
% Using code copied from functions written by Knud Kragh, 2011


function Matlab2AD(ADPar,TemplateFile,OutputFilename)


%Declare file IDs from the template to the resulting file
fidIN = fopen(TemplateFile,'r');
if fidIN == -1
    error([TemplateFile ' not found'])
end
fidOUT = fopen(OutputFilename,'w');
if fidOUT == -1
    error(['Cant create ' OutputFilename])
end


%loop through the template up until FoilNm
while true
    line = fgets(fidIN); %get the next line from the template
    
    if ~isempty(strfind(line,'FoilNm'))
        break;
    end
    
    %now see if there is a match to a known parameter
    idx = 9999; %9999 code is no match
    label2Write = []; %keep track of the label to make sure we find the longest matching label
    for i = 1:length(ADPar.Label)
        idxTemp = strfind(line,ADPar.Label{i});
        if ~isempty(idxTemp) %if there is a match,
            if idxTemp(1) < idx %if this match occurs before other matches use it
                
                label2Write = ADPar.Label{i}; %store this label
                idx =idxTemp(1);
                val2Write = ADPar.Val{i};
            
            elseif idxTemp(1) == idx %if this match equals other matches check if it is the longer match
                if (length(ADPar.Label{i}) > length(label2Write)) %and this label is longer than others (avoid coincident substring matches)
                    label2Write = ADPar.Label{i}; %store this label
                    idx =idxTemp(1);
                    val2Write = ADPar.Val{i};
                end
            end
        end
    end
    %if we found a match
    if (idx ~= 9999)
        if(line(1)~='-') %not a comment
            
            %build a string version of the value
            writeVal = num2str(val2Write);
            
            %PF: Not necessary in AD
%             %Check for the special case of an array
%             if isnumeric(val2Write)
%                 if ~isscalar(val2Write)
%                     writeVal = [];
%                     for ii = 1:length(val2Write)
%                         writeVal = [writeVal ',' num2str(val2Write(ii))];
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

%Now write out the Foilnm list

%The first line is a special case, with FoilNm appended
fprintf(fidOUT,'%s',ADPar.FoilNm{1});
fprintf(fidOUT,'\tFoilNm   - Names of the airfoil files [NumFoil lines] (quoted strings)\n');

%now loop through the rest
for i = 2:length(ADPar.FoilNm)
    fprintf(fidOUT,'%s',ADPar.FoilNm{i});
    fprintf(fidOUT,'\n');
end

%Now write the BldNodes line
fprintf(fidOUT,[num2str(GetFastPar(ADPar,'BldNodes')) '\tBldNodes - Number of blade nodes used for analysis (-)\n']);

%Write the BNodes HDR
fprintf(fidOUT,'RNodes\tAeroTwst\tDRNodes\tChord\tNFoil\tPrnElm\n');

%now loop through and pring out blade nodes matrix


for i = 1:length(ADPar.PrnElm)
    fprintf(fidOUT,[num2str(ADPar.BldNodes(i,1)),'\t',num2str(ADPar.BldNodes(i,2)),'\t',num2str(ADPar.BldNodes(i,3)),'\t',num2str(ADPar.BldNodes(i,4)),'\t',num2str(ADPar.BldNodes(i,5)),'\t',ADPar.PrnElm{i},'\n']);
end

%now print the word ReNum
fprintf(fidOUT,'ReNum');

fclose(fidIN);
fclose(fidOUT);
end %end function
