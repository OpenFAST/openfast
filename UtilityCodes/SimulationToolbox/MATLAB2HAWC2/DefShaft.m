% Function for defining Shaft main body and adding hub mass
%
% In: FastPar       -   FAST input parameters extracted using Fast2Matlab.m
%     Hawc2Set      -   Hawc2 input settings
%     OutFileNames  -   Name of Hawc2 st file
%     fid           -   File identifier of Hawc2 htc file
%
% Knud A. Kragh
function DefShaft(FastPar,Hawc2Set,OutFileNames,fid)

fprintf(fid,['  begin main_body; \n']);
fprintf(fid,['    name        shaft ;  \n']);
fprintf(fid,['    type        timoschenko ; \n']);
fprintf(fid,['    nbodies     ' num2str(Hawc2Set.Shaft.nbodies) ' ; \n']);
fprintf(fid,['    node_distribution     c2_def ; \n']);
fprintf(fid,['    damping   0.01 0.001 0.001 4.0E-03 1.0E-03 4.5E-04 ; \n']);
fprintf(fid,['    concentrated_mass ' num2str(Hawc2Set.Shaft.nSec+1) ' 0.0 0.0 0.0 ' num2str(GetFastPar(FastPar,'HubMass')) ' 0.0 0.0 ' num2str(GetFastPar(FastPar,'HubIner')) '; \n']);

fprintf(fid,['  begin timoschenko_input; \n']);
fprintf(fid,['    filename ./data/' OutFileNames '.st ; \n']);
fprintf(fid,['    set 3 1 ; \n']);
fprintf(fid,['  end timoschenko_input; \n']);

fprintf(fid,['  begin c2_def; \n']);
fprintf(fid,['    nsec ' num2str(Hawc2Set.Shaft.nSec+1) '; \n']);
fprintf(fid,['    sec ' num2str(1) '\t 0.0 0.0 0.0 0.0 ;  x,y,z,twist; \n']);
for i=1:Hawc2Set.Shaft.nSec
    fprintf(fid,['    sec ' num2str(i+1) '\t 0.0 0.0 ' num2str(-1*round(100*(GetFastPar(FastPar,'OverHang'))*i*1/Hawc2Set.Shaft.nSec)/100) '\t 0.0 ;  x,y,z,twist; \n']);
end
fprintf(fid,['  end c2_def; \n']);

fprintf(fid,['  end main_body; \n']);
fprintf(fid,['; \n']);


end
