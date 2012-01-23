function [data,time,colTxt,colUnits] = readHHbin(FileName,fileFmt)
%[data,time,colTxt,colUnits] = readHHbin(FileName)
%                       or     readHHbin(FileName,fileFmt)
%--------------------------------------------------------------------------
% Required Input:
%   FileName    - string     - name of the TurbSim HH binary file to read
%   
% Optional Input:
%   fileFmt     - string     - type of values in the binary file.  By default,
%                              this is 'float32', the type of TurbSim data.
% Output:
%   data        - 2D matrix  - This Nx13 matrix contains the 13 data columns
%                              in the file at each of the N time steps.
%   time        - vector     - This column vector of length N contains the
%                              time corresponding to the values in the
%                              "data" matrix.
%   colTxt      - 14 strings - These cells of strings are the text labels
%                              describing the 13 columns in the "data"
%                              matrix and time (in cell 14).  
%                              colTxt{i} describes data(:,i).
%   colUnits    - 14 strings - These cells of strings are the units for 
%                              the corresponding colTxt data.
%--------------------------------------------------------------------------
    if ( nargin < 2 )
        fileFmt = 'float32';
    end

    fid = fopen( FileName, 'r' );

    if fid > 0
        
        
        A = fread(fid, [14, inf], fileFmt);
%         A = fread(fid, [14, 1], fileFmt);

        data = A(2:14,:)';
        time = A(1,:)';

        colTxt  = {'U',    'Uh',   'Ut',   'V',    'W',    'u''',  'v''',  'w''',  'u''w''', 'u''v''', 'v''w''', 'TKE',    'CTKE',   'Time'};
        colUnits= {'(m/s)','(m/s)','(m/s)','(m/s)','(m/s)','(m/s)','(m/s)','(m/s)','(m/s)^2','(m/s)^2','(m/s)^2','(m/s)^2','(m/s)^2','(s)' }; 

        fclose(fid);
    else
        error( [' Unable to open ' FileName '.'] );
    end

return;