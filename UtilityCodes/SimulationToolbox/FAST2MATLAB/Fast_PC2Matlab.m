% Function for reading in FAST/AeroDyn PC file
%
% In:   PC_file    -   Full file name of PC file
%       ADver      -   Aerodyn version (12 or 13)
% Out   -   Data structure containg the [AoA Cl Cd] for all sets in foildata (Reynolds numbers)
%
% Knud A. Kragh

function DataOut=Fast_PC2Matlab(PC_file, ADver)

% clc;clear all; close all;
% PC_file='C:\Documents and Settings\knkr\My Documents\PhD\FAST\Cart3Model\AeroData\C3_05_S818_AD12.dat';
% ADver=12;

fid = fopen([PC_file],'r');
if fid == -1
    Flag = 0;
    disp('  ')
    disp('==============================================================')
    disp(['PC file could not be found'])
    disp('--------------------------------------------------------------')
    return
end

if ADver==13
    
    tline = fgets(fid); % Line n
    tline = fgets(fid); % Line n
    tline = fgets(fid); % Line n
    
    DataOut.nFoils=fscanf(fid,'%f',1);
    tline = fgets(fid); % Line n
    
    DataOut.Data=zeros(2,2,DataOut.nFoils);
    DataOut.ReyNo=[];
    DataOut.ReyNo=[DataOut.ReyNo fscanf(fid,'%f',1)*1e6];
    
    for i=1:DataOut.nFoils
        for ii=1:7
            tline = fgets(fid); % Line n
        end
        %for i=1:DataOut.nFoils
        count=1;
        tline = fgets(fid); % Line n
        
        while DataOut.Data(end,1,i)~=180
            DataOut.Data(count,1,i)=fscanf(fid,'%f',1);
            DataOut.Data(count,2,i)=fscanf(fid,'%f',1);
            DataOut.Data(count,3,i)=fscanf(fid,'%f',1);
            count=count+1;
        end
        tline = fgets(fid); % Line n
        tline = fgets(fid); % Line n
        
        DataOut.ReyNo=[DataOut.ReyNo fscanf(fid,'%f',1)*1e6];
        
    end
    
elseif ADver==12
    tline = fgets(fid); % Line n
    tline = fgets(fid); % Line n
    
    DataOut.nFoils=fscanf(fid,'%f',1);
    
    tline = fgets(fid); % Line n
    
    DataOut.ReyNo=[];
    for i=1:DataOut.nFoils
        DataOut.ReyNo=[DataOut.ReyNo fscanf(fid,'%f',1)*1e6];
    end
    for i=1:11
        tline = fgets(fid); % Line n
    end
    
    count=1;
    temp=zeros(2,2);
    while temp(end,1)<180
        for i=1:DataOut.nFoils*2+1
            test=fscanf(fid,'%f',1);
            
            if isempty(test)
                fscanf(fid,'%s',1);
                temp(count,i)=fscanf(fid,'%f',1);
            else
                temp(count,i)=test;
            end
        end
        tline = fgets(fid); % Line n
        count=count+1;
    end
    
    for i=1:DataOut.nFoils
        DataOut.Data(:,1,i)=temp(:,1);
        DataOut.Data(:,2,i)=temp(:,i*2);
        DataOut.Data(:,3,i)=temp(:,i*2+1);
    end
    
    
    
    
end
DataOut.DataAv=zeros(size(DataOut.Data(:,:,1)));
for i=1:DataOut.nFoils
    DataOut.DataAv=DataOut.DataAv+DataOut.Data(:,:,i);
end
DataOut.DataAv=DataOut.DataAv./DataOut.nFoils;
fclose(fid);
end