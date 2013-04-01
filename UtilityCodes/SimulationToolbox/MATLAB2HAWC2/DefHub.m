% Function for defining Hub main bodies
%
% In: FastPar       -   FAST input parameters extracted using Fast2Matlab.m
%     Hawc2Set      -   Hawc2 input settings
%     OutFileNames  -   Name of Hawc2 st file
%     fid           -   File identifier of Hawc2 htc file
%
% Knud A. Kragh
function DefHub(FastPar,Hawc2Set,OutFileNames,fid)

fprintf(fid,['  begin main_body; \n']);
fprintf(fid,['    name        hub1 ;  \n']);
fprintf(fid,['    type        timoschenko ; \n']);
fprintf(fid,['    nbodies     ' num2str(Hawc2Set.Shaft.nbodies) ' ; \n']);
fprintf(fid,['    node_distribution     c2_def ; \n']);
fprintf(fid,['    damping   5.0E-02 5.0E-02 8.0E-01 1.0E-03 1.0E-03 4.5E-04 ; \n']);

fprintf(fid,['  begin timoschenko_input; \n']);
fprintf(fid,['    filename ./data/' OutFileNames '.st ; \n']);
fprintf(fid,['    set 4 1 ; \n']);
fprintf(fid,['  end timoschenko_input; \n']);

fprintf(fid,['  begin c2_def; \n']);
fprintf(fid,['    nsec ' num2str(Hawc2Set.Hub.nSec+1) '; \n']);
fprintf(fid,['    sec ' num2str(1) '\t 0.0 0.0 0.0 0.0 ;  x,y,z,twist; \n']);
for i=1:Hawc2Set.Hub.nSec
    fprintf(fid,['    sec ' num2str(i+1) '\t 0.0 0.0 ' num2str(round(100*(GetFastPar(FastPar,'HubRad'))*i*1/Hawc2Set.Hub.nSec)/100) '\t 0.0 ;  x,y,z,twist; \n']);
end
fprintf(fid,['  end c2_def; \n']);

fprintf(fid,['  end main_body; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['  begin main_body; \n']);
fprintf(fid,['  name  hub2  ; \n']);
fprintf(fid,['    copy_main_body  hub1; \n']);
fprintf(fid,['    end main_body; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['  begin main_body; \n']);
fprintf(fid,['    name  hub3  ; \n']);
fprintf(fid,['    copy_main_body  hub1; \n']);
fprintf(fid,['  end main_body; \n']);


fprintf(fid,['; \n']);


end
