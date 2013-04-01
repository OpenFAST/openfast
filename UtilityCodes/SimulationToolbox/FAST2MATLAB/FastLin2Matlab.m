% FastLin2Matlab
% Function for reading FAST input files in to a MATLAB struct.
%
% Knud A. Kragh, May 2011, NREL, Boulder
%
% Modfied by Paul Fleming, JUNE 2011

%This function returns a structure LinParOut, which contains 2 cell arrays,
%.Val: An array of values
%.Label: An array of matching labels
%These arrays are extracted from the FAST input file

function LinParOut = FastLin2Matlab(Lin_file,hdrLines)

if nargin < 2
    hdrLines = 0;
end

%----------------------Read FST main file----------------------------------
fid = fopen(Lin_file,'r');
if fid == -1
    Flag = 0;
    error('FST file could not be found')
end

%skip hdr
for hi = 1:hdrLines
    fgets(fid);
end


%Loop through the file line by line, looking for value-label pairs
count=1;
while ~feof(fid)
    fgets(fid); %Advanced to the next line
    skipLine = false; %reset skipline
    %Label=[]; %Re-initialize label  PF: Temp disabling this
    
    
    % Get the Value, number or string
    testVal=fscanf(fid,'%f',1);  %First check if line begins with a number
    if isempty(testVal)
        testVal=fscanf(fid,'%s',1);  %If not look for a string instead
        
        if isempty(testVal) %if there's nothing we're at EOF
            break;
        end
        
        %now check to see if the string in test val makes sense as a value
        if strcmpi(testVal,'false') || strcmpi(testVal,'true') || testVal(1)=='"'
            %disp(testVal) %this is a parameter
        else
            skipLine = true;
            if testVal(1)~='-' %test if this non parameter not a comment
                LinParOut.Label{count}=testVal;  %if not a comment, make the value the label
                LinParOut.Val{count}=' ';
                count=count+1;
            end
            
        end
    end
    
    
    
    if ~skipLine %if this is actually a parameter line add it
        LinParOut.Val{count}=testVal; %assign Val
        
        
        % Now get the label or possibly the remaining numbers in this list
        test=0;
        while test==0
            testVal=fscanf(fid,'%f',1);
            if isempty(testVal) %if we've reached something besides a number this is the label
                testVal=fscanf(fid,'%s',1);
                test=1;
            else %if we've hit another number append it LinParOut.Val
                LinParOut.Val{count} = [LinParOut.Val{count} testVal];
            end
        end
        LinParOut.Label{count}=testVal; %Save label
        
        count=count+1;
    end %endif
end %end while


fclose(fid); %close file

end %end function