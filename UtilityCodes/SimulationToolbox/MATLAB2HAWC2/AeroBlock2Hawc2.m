% Function for writing aero block if htc file
%
% In:   fid             -   File indentifier of file in which to write
%       OutFileNames    -   Output file names, name for ae and pc file
%                           without extension
%       Hawc2Set        -   Hawc2 input settings structure
%
%
%
% Knud A. Kragh

function AeroBlock2Hawc2(fid,OutFileNames,Hawc2Set)

fprintf(fid,['begin aero ;\n']);
fprintf(fid,['  nblades  3;\n']);
fprintf(fid,['  hub_vec shaft -3 ;          vector from hub (normal to rotor plane) directed towards tower top\n']);
fprintf(fid,['  link 1 mbdy_c2_def blade1;\n']);
fprintf(fid,['  link 2 mbdy_c2_def blade2;\n']);
fprintf(fid,['  link 3 mbdy_c2_def blade3;\n']);
fprintf(fid,['  ae_filename        ./data/' OutFileNames '.ae ;\n']);
fprintf(fid,['  pc_filename        ./data/' OutFileNames '.pc ;\n']);
fprintf(fid,['  induction_method   1 ;     0=none, 1=normal\n']);
fprintf(fid,['  aerocalc_method    1 ;     0=ingen aerodynamic, 1=med aerodynamic\n']);
fprintf(fid,['  aero_distribution       ae_file 1 ;\n']);
fprintf(fid,['  ae_sets            1 1 1;\n']);
fprintf(fid,['  tiploss_method     1 ;     0=none, 1=normal\n']);
fprintf(fid,['  dynstall_method    2 ;     0=none, 1=stig øye method,2=mhh method\n']);
fprintf(fid,['end aero ;\n']);
