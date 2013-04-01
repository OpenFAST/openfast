% Script for writing structure block in HAWC2 htc file from FAST data
% stored in matlab
%
%
%
% Knud A. Kragh

fprintf(fid,['begin new_htc_structure; \n']);
% Defining structure analysis output files
fprintf(fid,['  beam_output_file_name  ./logfiles/beam_' OutFileNames '.dat;                    Optional - Calculated beam properties of the bodies are written to file; \n']);
fprintf(fid,['  body_output_file_name  ./logfiles/body_' OutFileNames '.dat;                    Optional - Body initial position and orientation are written to file \n']);
fprintf(fid,[';  body_eigenanalysis_file_name ./eigenfrq/body_' OutFileNames '.dat; \n']);
fprintf(fid,['  structure_eigenanalysis_file_name ./eigenfrq/structure_' OutFileNames '.dat; \n']);
fprintf(fid,['; \n']);

% Defining main bodies
DefTower(FastPar,Hawc2Set,TowerData,OutFileNames,fid);
DefShaft(FastPar,Hawc2Set,OutFileNames,fid);
DefHub(FastPar,Hawc2Set,OutFileNames,fid);
DefBlade(FastPar,Hawc2Set,BladeData,OutFileNames,fid);

% Defining orientations
DefOrient(FastPar,fid)

% Defining constraints
DefConstraint(FastPar,fid)


fprintf(fid,['end new_htc_structure; \n']);