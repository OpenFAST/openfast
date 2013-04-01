% Linear2Matlab
% Function for reading linearization config file into a MATLAB struct.
%
%
%This function returns a structure DataOut, which contains
%
%   .Val: A cell array of values
%   .Label: A cell array of matching labels
%
% Jason Laks, Jan 2013

function DataOut = Linear2Matlab(Lin_file,hdrLines)

if nargin < 2
    hdrLines = 0;
end

%----------------------Read FST main file----------------------------------
fid = fopen([Lin_file],'r');
if fid == -1
    Flag = 0;
    error('Lin file could not be found')
end

%skip hdr
for hi = 1:hdrLines
    fgets_result=fgets(fid);    %fgets returns -1 at EOF
end
if ~exist('fgets_result','var')
    fgets_result=1;
end;
if fgets_result==-1
    error('End of linearization template reached during load!');
end;

%Loop through the file line by line, looking for value-label pairs
%Stop once we've reached NInputs, meaning we've reached the list of control
%inputs:
count=1;

past_header=false;  %set true when it appears header has been passed

while fgets_result~=-1 %loop until discovering NInputs, or until EOF
    
    
    % Get the Value, number or string
    testVal=fscanf(fid,'%f',1);  %First check if line begins with a number
    if isempty(testVal)
        testVal=fscanf(fid,'%s',1);  %If not look for a string instead
    elseif (past_header)&&(count==1)
        error('The first item in the Lin file should not be numeric');
    end

    if ~past_header %first item should be a true or false
        past_header=( strcmpi('true',testVal)||strcmpi('false',testVal) );
        if ~past_header
            fgets_result=fgets(fid); %Advanced to the next line
        end;
    end;
    if ischar(testVal)
        not_separator=isempty(strfind(testVal,'--'));
    else
        not_separator=true;
    end
    if past_header && not_separator
        %assign value
        DataOut.Val{count}=testVal; %assign Val

        %Now get the label
        % Now get the label, some looping is necessary because often
        % times, old values for FAST parameters are kept next to new
        % ones seperated by a space and need to be ignored
        test=0;
        while test==0
            testVal=fscanf(fid,'%f',1);
            if isempty(testVal) %if we've reached something besides a number
                testVal=fscanf(fid,'%s',1);
                test=1;
            end
        end
        DataOut.Label{count}=testVal; %Now save label

        fgets_result=fgets(fid); %Advanced to the next line
        count=count+1;

        %check if the last label read is the NInputs
        if strmatch(DataOut.Label{count-1},'NInputs')
            break; %if it does end the loop
        end
    end;
end %end while
if fgets_result==-1
    error('End of linearization template reached during load!');
end;


%% Read in Control Input list
NInputs = GetFastPar(DataOut,'NInputs');
if NInputs>0
    testVal=fscanf(fid,'%f',NInputs);
    if isempty(testVal) %if we've reached something besides a number
        error('Could not find CntrlInpt');
    elseif length(testVal)~=NInputs
        error('Did not find the correct number of CntrlInpt');
    end
    DataOut.Val{count}=testVal';    %Matlab2FASTLin wants row vect
else
    DataOut.Val{count}=1;
end

%now get label
test=0;
while test==0
    testVal=fscanf(fid,'%f',1);
    if isempty(testVal) %if we've reached something besides a number
        testVal=fscanf(fid,'%s',1);
        test=1;
    end
end
DataOut.Label{count}=testVal; %Now save label

fgets_result=fgets(fid); %Advanced to the next line
if fgets_result==-1
    error('End of linearization template reached during load!');
end;
count=count+1;

% Get the NDisturbs
testVal=fscanf(fid,'%f',1);  %First check if line begins with a number
if isempty(testVal)
    error('Could not find NDisturbs');
end

%assign value
DataOut.Val{count}=testVal; %assign Val

% Now get the label, some looping is necessary because often
% times, old values for FAST parameters are kept next to new
% ones seperated by a space and need to be ignored
test=0;
while test==0
    testVal=fscanf(fid,'%f',1);
    if isempty(testVal) %if we've reached something besides a number
        testVal=fscanf(fid,'%s',1);
        test=1;
    end
end
DataOut.Label{count}=testVal; %Now save label

fgets_result=fgets(fid); %Advanced to the next line
if fgets_result==-1
    error('End of linearization template reached during load!');
end;
count=count+1;

%% Read in Disturbance Input list
NDisturbs = GetFastPar(DataOut,'NDisturbs');
if NInputs>0
    testVal=fscanf(fid,'%f',NDisturbs);
    if isempty(testVal) %if we've reached something besides a number
        error('Could not find Disturbnc');
    elseif length(testVal)~=NDisturbs
        error('Did not find the correct number of Disturbnc');
    end
    DataOut.Val{count}=testVal';        %Matlab2FASTLin wants row vect
else
    DataOut.Val{count}=1;
end

% Now get the label, some looping is necessary because often
% times, old values for FAST parameters are kept next to new
% ones seperated by a space and need to be ignored
test=0;
while test==0
    testVal=fscanf(fid,'%f',1);
    if isempty(testVal) %if we've reached something besides a number
        testVal=fscanf(fid,'%s',1);
        test=1;
    end
end
DataOut.Label{count}=testVal; %Now save label

fclose(fid); %close file

end %end function