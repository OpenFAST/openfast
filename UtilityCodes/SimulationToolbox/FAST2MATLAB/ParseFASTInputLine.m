function [value, label, isComment, descr, fieldType] = ParseFASTInputLine( line )
% This routine parses a line from a FAST type input file. 
% Comment lines are assumed to start with any of the following individual characters:
%   #!= 
% or this combination of characters:
%  --
% If the line is not a comment, it is assumed to be of the form:
%  value [,Array values] <old values> label descr
% If the value is an array, it must be separated by commas.
% If there are multiple values separated by white space, it is assumed to
%   contain old values (instead of an array)
%--------------------------------------------------------------------------
% Inputs:
%   line        - a line of text
% Outputs:
%   value       - the value of the parameter
%   label       - the name of the parameter/variable/field
%   isComment   - logical value that says if the line is a comment line
%   descr       - the description of the line
%   fieldType   - text saying the field is either a
%                   "Comment", "Logical", "Character", "Numeric" variable
%--------------------------------------------------------------------------

        % first check that this isn't a blank line...
    if isempty(line) || length(strtrim(line)) < 1
        value = '';
        label = '';
        isComment = true;
        descr = '';
        fieldType = 'Comment';
        return
    end

    trueFalseValues = {'true','false','t','f'};

    % determine if this is a comment line:
    first2Chars = sscanf(strtrim(line),'%c',2); %read the first two not-starting-with-whitespace characters
    firstChar   = first2Chars(1);
    
    if ~isempty( strfind( '#!=', firstChar ) ) || strcmp( first2Chars, '--' ) %comments start with any of these characters: # ! - =
        value = strtrim(line);
        label = '';
        isComment = true;
        descr = value;
        fieldType = 'Comment';
    else
        isComment = false;

        % Get the Value, number or string
        [value, cnt, ~, nextindex] = sscanf(line,'%f', 1);  %First check if line begins with a number

        if cnt == 0 % we didn't find a number so...
            [testVal, position] = textscan(line,'%q',1);  %look for a string instead
            if any( strcmpi(testVal{1}{1},trueFalseValues) )
                value = testVal{1}{1};                %this is a logical input
                fieldType = 'Logical';            
            else
                value = ['"' testVal{1}{1} '"']; %add quotes back around the string
                fieldType = 'Character';            
            end
            nextindex = position + 1;
        else
            fieldType = 'Numeric';            
        end

        % Now get the label     

            % Some looping is necessary because often times,
            % old values for FAST parameters and kept next to new
            % ones seperated by a space and need to be ignored


        IsLabel = false;
        label   = '';  %initialize in case the line doesn't have a label
        descr   = '';
        
        while ~IsLabel && length(line) >= nextindex
            line = line(nextindex:end);

            [~, cnt, ~, nextindex] =sscanf(line,'%f',1);
            if cnt == 0 %if we've reached something besides a number - or we're at the end of the line

                [testVal, cnt, ~, nextindex] = sscanf(line,'%s',1);
                if cnt == 1
                    if any( strcmpi(testVal,trueFalseValues) )
                        %this is a logical input
                    elseif strcmpi(testVal(1),',') 
                        % commas are an indication that this parameter is a list
                        if strcmpi(fieldType, 'Numeric') 
                            [testVal, cnt, ~, nextindex] = sscanf(line(2:end),'%f',1);
                            if cnt == 1
                                value = [value; testVal];
                            end
                        end
                    else
                        IsLabel = true;
                        label = testVal;
                        descr = strtrim(line(nextindex:end));
                    end
                end
            end

        end %while 
            
    end %not a comment
           
    
return