% Function for defining Blade main bodies in htc file
%
% In: FastPar       -   FAST input parameters extracted using Fast2Matlab.m
%     Hawc2Set      -   Hawc2 input settings
%     BladeData     -   Blade data in extracted from FAST tower file
%     OutFileNames  -   Name of Hawc2 st file
%     fid           -   File identifier of Hawc2 htc file
%
% Knud A. Kragh

function DefBlade(FastPar,Hawc2Set,BladeData,OutFileNames,fid)

fprintf(fid,['  begin main_body; \n']);
fprintf(fid,['    name        blade1 ;  \n']);
fprintf(fid,['    type        timoschenko ; \n']);
fprintf(fid,['    nbodies     ' num2str(Hawc2Set.Blade.nbodies) ' ; \n']);
fprintf(fid,['    node_distribution     c2_def ; \n']);
fprintf(fid,['    damping   5.0E-02 5.0E-02 8.0E-01 1.3E-03 2.0E-03 4.5E-04 ; \n']);

fprintf(fid,['  begin timoschenko_input; \n']);
fprintf(fid,['    filename ./data/' OutFileNames '.st ; \n']);
fprintf(fid,['    set 5 1 ; \n']);
fprintf(fid,['  end timoschenko_input; \n']);

fprintf(fid,['  begin c2_def; \n']);
fprintf(fid,['    nsec ' num2str(GetFastPar(BladeData,'NBlInpSt')) '; \n']);
for i=1:GetFastPar(BladeData,'NBlInpSt')
    fprintf(fid,['    sec ' num2str(i) '\t' num2str(-BladeData.BldProp(i,13)) '\t' num2str(BladeData.BldProp(i,12)) '\t' num2str(round(100*((GetFastPar(FastPar,'TipRad')-GetFastPar(FastPar,'HubRad'))*BladeData.BldProp(i,1)))/100) '\t' num2str(-BladeData.BldProp(i,3)) ';  x,y,z,twist; \n']);
end
fprintf(fid,['  end c2_def; \n']);


fprintf(fid,['  end main_body; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['  begin main_body; \n']);
fprintf(fid,['    name  blade2  ; \n']);
fprintf(fid,['    copy_main_body  blade1; \n']);
fprintf(fid,['  end main_body; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['  begin main_body; \n']);
fprintf(fid,['    name  blade3  ; \n']);
fprintf(fid,['    copy_main_body  blade1; \n']);
fprintf(fid,['  end main_body; \n']);


fprintf(fid,['; \n']);


end
