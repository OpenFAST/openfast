% Function for writing Hawc2 st file from ST in Hawc2 format stored in
% matlab structure
%
% In:   StData          -   Structure containing structural data for Tower,
%                           shaft, hub, and blade in Hawc2 format
%       OutPath         -   Path to folder containing Hawc2 files
%       OutFileNames    -   Name of output st-file
%
% Out:  Hawc2 .st file placed in OutPath\Data\
%
%
%
% Knud A. Kragh

function st2Hawc2(StData,OutPath,OutFileNames)

mkdir(OutPath,'Data');
% Writing ae file
fid = fopen([OutPath 'Data\' OutFileNames '.st' ], 'w+');

fprintf(fid,'4  Nset \n');
fprintf(fid,'This file is converted from FAST, fit values!! \n');

%% Writing tower data
% Flexible
fprintf(fid,'#1  Tower\n');
fprintf(fid,'r			m		x_cg		y_cg ri_x		ri_y		x_sh		y_sh	E		G	I_x		I_y		I_p	k_x  k_y	A		pitch	x_e	y_e \n');
fprintf(fid,['$1 ' num2str(size(StData.Tower,1))  ' Flexible \n']);

for i=1:size(StData.Tower,1) % Rows
    for ii=1:19 % columns
        fprintf(fid,[num2str((100*StData.Tower(i,ii))/100) '\t \t']);
    end
    fprintf(fid,['\n']);
end

% Stiff
fprintf(fid,['$2 ' num2str(size(StData.Tower,1))  ' Stiff \n']);

for i=1:size(StData.Tower,1) % Rows
    for ii=1:19 % columns
        if ii==9
            fprintf(fid,[num2str((100*StData.Tower(i,ii)*10e5)/100) '\t \t']);
        elseif ii==10
            fprintf(fid,[num2str((100*StData.Tower(i,ii)*10e5)/100) '\t \t']);
        elseif ii==13
            fprintf(fid,[num2str((100*StData.Tower(i,ii)*100)/100) '\t \t']);
        else
            fprintf(fid,[num2str((100*StData.Tower(i,ii))/100) '\t \t']);
        end
    end
    fprintf(fid,['\n']);
end

%% Writing shaft data
% Flexible
fprintf(fid,'#3  Shaft \n');
fprintf(fid,'r			m		x_cg		y_cg ri_x		ri_y		x_sh		y_sh	E		G	I_x		I_y		I_p	k_x  k_y	A		pitch	x_e	y_e \n');
fprintf(fid,['$1 ' num2str(size(StData.Shaft,1))  ' Flexible \n']);

for i=1:size(StData.Shaft,1) % Rows
    for ii=1:19 % columns
        fprintf(fid,[num2str((100*StData.Shaft(i,ii))/100) '\t \t']);
    end
    fprintf(fid,['\n']);
end

% Stiff
fprintf(fid,['$2 ' num2str(size(StData.Shaft,1))  ' Stiff \n']);

for i=1:size(StData.Shaft,1) % Rows
    for ii=1:19 % columns
        if ii==9
            fprintf(fid,[num2str((100*StData.Shaft(i,ii)*10e5)/100) '\t \t']);
        elseif ii==10
            fprintf(fid,[num2str((100*StData.Shaft(i,ii)*10e5)/100) '\t \t']);
        elseif ii==13
            fprintf(fid,[num2str((100*StData.Shaft(i,ii)*100)/100) '\t \t']);
        else
            fprintf(fid,[num2str((100*StData.Shaft(i,ii))/100) '\t \t']);
        end
    end
    fprintf(fid,['\n']);
end




%% Writing hub data
% Flexible
fprintf(fid,'#4  Hub \n');
fprintf(fid,'r			m		x_cg		y_cg ri_x		ri_y		x_sh		y_sh	E		G	I_x		I_y		I_p	k_x  k_y	A		pitch	x_e	y_e \n');
fprintf(fid,['$1 ' num2str(size(StData.Hub,1))  ' Flexible \n']);

for i=1:size(StData.Hub,1) % Rows
    for ii=1:19 % columns
        fprintf(fid,[num2str((100*StData.Hub(i,ii))/100) '\t \t']);
    end
    fprintf(fid,['\n']);
end

% Stiff
fprintf(fid,['$2 ' num2str(size(StData.Hub,1))  ' Stiff \n']);

for i=1:size(StData.Hub,1) % Rows
    for ii=1:19 % columns
        if ii==9
            fprintf(fid,[num2str((100*StData.Hub(i,ii)*10e5)/100) '\t \t']);
        elseif ii==10
            fprintf(fid,[num2str((100*StData.Hub(i,ii)*10e5)/100) '\t \t']);
        elseif ii==13
            fprintf(fid,[num2str((100*StData.Hub(i,ii)*100)/100) '\t \t']);
        else
            fprintf(fid,[num2str((100*StData.Hub(i,ii))/100) '\t \t']);
        end
    end
    fprintf(fid,['\n']);
end

%% Writing blade data
%Flexible
fprintf(fid,'#5  Blade \n');
fprintf(fid,'r			m		x_cg		y_cg ri_x		ri_y		x_sh		y_sh	E		G	I_x		I_y		I_p	k_x  k_y	A		pitch	x_e	y_e \n');
fprintf(fid,['$1 ' num2str(size(StData.Blade,1))  ' Flexible \n']);

for i=1:size(StData.Blade,1) % Rows
    for ii=1:19 % columns
        fprintf(fid,[num2str((100*StData.Blade(i,ii))/100) '\t \t']);
    end
    fprintf(fid,['\n']);
end

% Stiff
fprintf(fid,['$2 ' num2str(size(StData.Blade,1))  ' Stiff \n']);

for i=1:size(StData.Blade,1) % Rows
    for ii=1:19 % columns
        if ii==9
            fprintf(fid,[num2str((100*StData.Blade(i,ii)*10e5)/100) '\t \t']);
        elseif ii==10
            fprintf(fid,[num2str((100*StData.Blade(i,ii)*10e5)/100) '\t \t']);
        elseif ii==13
            fprintf(fid,[num2str((100*StData.Blade(i,ii)*100)/100) '\t \t']);
        else
            fprintf(fid,[num2str((100*StData.Blade(i,ii))/100) '\t \t']);
        end
    end
    fprintf(fid,['\n']);
end

fclose(fid);