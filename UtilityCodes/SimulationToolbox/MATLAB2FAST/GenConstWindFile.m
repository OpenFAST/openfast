% GenConstWindFile
% Function for generating a constant wind speed file given:
% 1) A base filename
% 2) A windspeed
% - The generated file will have a name BaseName_Windspeed_MPS.wnd
%
%In:    baseName          -   directory and basename of generated wind file
%       windSpeed         -   wind speed in mps
%
%Out:   fileName          -   final filename of generated const wind file
%
% Paul Fleming, JUNE 2011
% Using code copied from functions written by Jan-Willem van Wingerden in
% July 2010


function fileName = GenConstWindFile(baseName,windSpeed)

%Build the final filename
intSpeed = num2str(floor(windSpeed));
tenthSpeed = num2str(floor(mod(windSpeed*10,10)));
fileName = [baseName '_' intSpeed '_' tenthSpeed 'MPS.wnd'];


%Now open the file and write the constant wind file
fid= fopen(fileName,'w');

fprintf(fid,'!Constant wind speed file. \n');
fprintf(fid,'!Time  Wind     Wind	Vert.       Horiz.      Vert.       LinV        Gust \n');
fprintf(fid,'!      Speed    Dir    Speed       Shear		Shear       Shear       Speed \n');
fprintf(fid,'%2.2f  %2.2f   %2.2f   %2.2f       %2.2f       %2.2f       %2.2f       %2.2f\n', 0, windSpeed, 0, 0,0, 0,0,0);
fprintf(fid,'%2.2f  %2.2f   %2.2f   %2.2f       %2.2f       %2.2f       %2.2f       %2.2f\n', 0.1, windSpeed, 0, 0,0, 0,0,0);
fprintf(fid,'%2.2f  %2.2f   %2.2f   %2.2f       %2.2f       %2.2f       %2.2f       %2.2f\n', 999, windSpeed, 0, 0,0, 0,0,0);

fclose(fid);

end