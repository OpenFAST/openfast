% Function for writing wind block
%
%
%
% Knud A. Kragh

function Wind2Hawc2(fid,windpar)

fprintf(fid,['begin wind ;\n']);
fprintf(fid,['  density                 1.225 ;\n']);
fprintf(fid,['  wsp                     ' num2str(windpar.wsp) '  ;\n']);
fprintf(fid,['  tint                    ' num2str(windpar.tint) ' ;\n']);
fprintf(fid,['  horizontal_input        1     ;\n']);
fprintf(fid,['  windfield_rotations     ' num2str(windpar.wrot) ' 0.0  0.0 ;\n']);
fprintf(fid,['  center_pos0             0.0 0.0 ' num2str(windpar.hubheight) ' ;  \n']);
fprintf(fid,['  shear_format            ' num2str(windpar.sformat) '\t '   num2str(windpar.shearexp) ' ;  0=none, 1=constant, 2=log, 3=power, 4=linear\n']);
fprintf(fid,['  turb_format             ' num2str(windpar.turbform) '       ;  0=none, 1=mann, 2=flex\n']);
fprintf(fid,['  tower_shadow_method     0       ;  0=none, 1=potential flow, 2=jet\n']);
fprintf(fid,[';\n']);
fprintf(fid,['  begin mann;\n']);
fprintf(fid,['    create_turb_parameters 33.6 1.0 3.9 ' num2str(windpar.turbseed) ' 1.0 ;\n']);
fprintf(fid,['    filename_u    .\dturb\du.bin ;      \n']);
fprintf(fid,['    filename_v    .\dturb\dv.bin ;        \n']);
fprintf(fid,['    filename_w    .\dturb\dw.bin ;      \n']);
fprintf(fid,['    box_dim_u   16384 1.031494 ;\n']);
fprintf(fid,['    box_dim_v     32 2.8125 ;\n']);
fprintf(fid,['    box_dim_w     32 2.8125 ;          \n']);
fprintf(fid,['    std_scaling   1.0 0.7 0.5 ;\n']);
fprintf(fid,['  end mann;\n']);
fprintf(fid,[';\n']);
fprintf(fid,['end wind;\n']);