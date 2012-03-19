function [u, t, txt, units] = loadColumnData(fname, uc, tab, rSkip, tr)
% function [u, t, txt, units] = loadColumnData(fname, uc, tab, rSkip, tr)
% This function reads column output and returns a time vector and matrices
% of that output.
%
% Required Input:
%    fname  - the name of the file to be opened
%
% Optional Input:
%    uc      - the vector of column numbers to be output, default is [2,3,4]
%              it can also contain the string 'all' for all columns
%    tab     - the column delimiter, default is '\t' ( i.e. tab-delimited )
%    rSkip   - a vector containing 
%               (1) the number of rows to skip before reading in numeric
%                   data, default is 8 (as most FAST headers)
%               (2) the number of columns to skip before reading data,
%                   default is 0
%               (3) the last row to read, default is the last row
%               (4) the last column to read, default is the last column
%              rSkip can have a length 1,2,or 4. See dlmread() for details.
%    tr      - the number of rows to skip before reading titles, default is
%               max( rSkip(1) - 2, 0 )
%
% Output:
%    u       - the data from columns uc
%    t       - the time data, found in column 1
%    txt     - a cell array containing the text for columns uc, from row tr+1 
%              (text for time data is in the last output cell)
%    units   - a cell array containing the units for columns uc, from row
%              tr+2 (units for time data is in the last output cell)

num_args   = nargin;
num_argOut = nargout;

switch num_args
    case 1
        uc    = [2:4];
        tab   = '\t';
        rSkip = 8;
        tr    = rSkip-2;
    case 2
        tab   = '\t';
        rSkip = 8;
        tr    = rSkip-2;
    case 3
        rSkip = 8;
        tr    = rSkip-2;
    case 4
        tr    = rSkip(1)-2;
    case 5
    otherwise
        error('Invalid number of input arguments in loadColumnData().');
end  
%%
specificCols = isnumeric(uc);
    
if (specificCols)
    if ( all( size(uc) > 1 ) )
        error('Invalid column numbers, uc, in loadColumnData().');
    end
else
    uc = 1;
end

tc  = 1;                                                % time is always the first column

sz = size(rSkip);
if all( sz > 1 )
    error('Invalid range numbers, rSkip, in loadColumnData().');
end

sz = prod(sz);
if sz > 1
    nColSkip = rSkip(2);
    
    if sz == 4
        M = dlmread(fname, tab, rSkip);                 % assume rSkip is a range of data
    elseif sz == 2
        M = dlmread(fname, tab, rSkip(1), rSkip(2));    % assume rSkip is row and column starting place
    else
        error('Invalid range numbers, rSkip, in loadColumnData().');
    end
else
    nColSkip = 0;
    
    M   = dlmread(fname, tab, rSkip, 0);                % read tab-delimited file, skipping rSkip rows 
end

msz = size(M,2);
if ~specificCols
    uc = 2:msz;    
end

mx  = max( [tc; uc(:)] );

if mx > msz
    error ([ 'The file ', fname, ' does not contain ', int2str(mx), ' columns. Only ', int2str(msz), ' columns were found.']);
end

% Get the requested columns
t  = M(:,tc);
u  = M(:,uc);
%%
if num_argOut > 2
	% Get text strings from column headers
    nCells = length(uc)+1;
    
	txt = cell(nCells,1);
    
    if num_argOut > 3   
        units = cell(nCells,1);
    end
	
	fid = fopen(fname);
	
	if (fid > 0)
        tr = max(0, tr);
        
        tmp = '';
        for n = 1:tr+1
            tmp = fgetl(fid);
        end
        tmpUnit = fgetl(fid);
        
        fclose(fid);

        if length(tmp) < 1
            return;
        end
        
        if tmp(1) == '!' && tmpUnit(1) == '!'
            tmp     = tmp(2:length(tmp));        
            tmpUnit = tmpUnit(2:length(tmpUnit));
        end 
        
        mx = max(uc)+nColSkip;
        
        if isempty(tab)
            tmp1 = textscan( tmp, '%s', mx );
            tmp1 = tmp1{1};
            
            if num_argOut > 3 
                tmp2 = textscan( tmpUnit, '%s', mx );                        
                tmp2 = tmp2{1};
            end            
        else
            if tab(1) == '\'
                if tab(2) == 't'
                    tab = char(9);       % The ASCII tab character
                elseif tab(2) == 'n'
                    tab = char(13);      % The ASCII carriage return character (what is the ASCII newline character?)
                end
            end
            
            n = 1;
            while ( n <= mx )
           
                [tmp1{n}, tmp] = strtok ( tmp, tab );
                
                if num_argOut > 3   
                    [tmp2{n}, tmpUnit] = strtok ( tmpUnit, tab );
                end
                                      
                if isempty(tmp) || isempty(tmpUnit)
                    n = mx + 1;            
                end
                n = n + 1;
                
            end %while
        end 

        n = length(tmp1);
                
        for i = 1:nCells-1           
            NewIX = uc(i) + nColSkip;
            if NewIX <= n            
                txt{i} = strtrim(tmp1{NewIX});
            else
                txt{i} = 'BLANK';
            end
            if num_argOut > 3   
                if NewIX <= length(tmp2) && rSkip(1) > tr + 1
                    units{i} = strtrim(tmp2{NewIX});
                else
                    units{i} = '(-)';
                end
            end
        end
        txt{nCells} = strtrim(tmp1{1});
        
        if num_argOut > 3  
            if rSkip(1) > tr + 1
                units{nCells} = strtrim(tmp2{1});
            else
                units{nCells} = '(-)';
            end                
        end
                       
	end % fid > 0
end % num_argOut > 2
    
return;
