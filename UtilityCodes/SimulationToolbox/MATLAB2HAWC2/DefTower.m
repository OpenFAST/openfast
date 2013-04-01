% Function for defining Tower main body in Hawc2 htc file and adding
% nacelle mass
%
% In: FastPar       -   FAST input parameters extracted using Fast2Matlab.m
%     Hawc2Set      -   Hawc2 input settings
%     TowerData     -   Tower data in extracted from FAST tower file
%     OutFileNames  -   Name of Hawc2 st file
%     fid           -   File identifier of Hawc2 htc file
%
% Knud A. Kragh
function DefTower(FastPar,Hawc2Set,TowerData,OutFileNames,fid)

fprintf(fid,['  begin main_body; \n']);
fprintf(fid,['    name        tower ;  \n']);
fprintf(fid,['    type        timoschenko ; \n']);
fprintf(fid,['    nbodies     ' num2str(Hawc2Set.Tower.nbodies) ' ; \n']);
fprintf(fid,['    node_distribution     c2_def ; \n']);
fprintf(fid,['    damping   0.01 0.001 0.001 4.0E-03 1.0E-03 4.5E-04 ; \n']);
fprintf(fid,['    concentrated_mass ' num2str(GetFastPar(TowerData,'NTwInpSt')) ' ' num2str(GetFastPar(FastPar,'NacCMyn')) ' ' num2str(GetFastPar(FastPar,'NacCMxn')) ' 0.0 ' num2str(GetFastPar(FastPar,'NacMass')) ' 0.0 0.0 ' num2str(GetFastPar(FastPar,'NacYIner')) '; \n']);

fprintf(fid,['  begin timoschenko_input; \n']);
fprintf(fid,['    filename ./data/' OutFileNames '.st ; \n']);
fprintf(fid,['    set 1 1 ; \n']);
fprintf(fid,['  end timoschenko_input; \n']);

fprintf(fid,['  begin c2_def; \n']);
fprintf(fid,['    nsec ' num2str(GetFastPar(TowerData,'NTwInpSt')) '; \n']);
for i=1:GetFastPar(TowerData,'NTwInpSt')
    fprintf(fid,['    sec ' num2str(i) '\t 0.0 0.0 ' num2str(round(100*-(GetFastPar(FastPar,'TowerHt')+GetFastPar(FastPar,'Twr2Shft'))*TowerData.TowProp(i,1))/100) '\t 0.0 ;  x,y,z,twist; \n']);
end
fprintf(fid,['  end c2_def; \n']);

fprintf(fid,['  end main_body; \n']);
fprintf(fid,['; \n']);


end
