% Define and writing constraints in Hawc2 htc file
%
% In: FastPar       -   FAST input parameters extracted using Fast2Matlab.m
%     fid           -   File identifier of Hawc2 htc file
%
%
% Knud A. Kragh

function DefConstraint(FastPar,fid)

fprintf(fid,['  begin constraint; \n']);

fprintf(fid,['  begin fix0;  fixed to ground in translation and rotation of node 1 \n']);
fprintf(fid,['    body tower; \n']);
fprintf(fid,['  end fix0; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['  begin bearing1;  Free bearing \n']);
fprintf(fid,['    name  shaft_rot ; \n']);
fprintf(fid,['    body1 tower last; \n']);
fprintf(fid,['    body2 shaft 1; \n']);
fprintf(fid,['    bearing_vector 2 0.0 0.0 -1.0; \n']);
fprintf(fid,['  end bearing1; \n']);
fprintf(fid,['; \n']);

fprintf(fid,[';  begin bearing2;  Bearing where rotation is controlled from external dll\n']);
fprintf(fid,[';    name  shaft_rot ; \n']);
fprintf(fid,[';    body1 tower last; \n']);
fprintf(fid,[';    body2 shaft 1; \n']);
fprintf(fid,[';    bearing_vector 2 0.0 0.0 -1.0; \n']);
fprintf(fid,[';  end bearing2; \n']);
fprintf(fid,['; \n']);

fprintf(fid,[';  begin bearing3;  Bearing where rotation is controlled from external dll\n']);
fprintf(fid,[';    name  shaft_rot ; \n']);
fprintf(fid,[';    body1 tower last; \n']);
fprintf(fid,[';    body2 shaft 1; \n']);
fprintf(fid,[';    bearing_vector 2 0.0 0.0 -1.0; \n']);
fprintf(fid,[';    omegas ' num2str(GetFastPar(FastPar,'RotSpeed')/60*2*pi) '; \n']);
fprintf(fid,[';  end bearing3; \n']);
fprintf(fid,['; \n']);
 
fprintf(fid,['  begin fix1;   \n']);
fprintf(fid,['    body1 shaft last; \n']);
fprintf(fid,['    body2 hub1 1; \n']);
fprintf(fid,['  end fix1; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['  begin fix1;   \n']);
fprintf(fid,['    body1 shaft last; \n']);
fprintf(fid,['    body2 hub2 1; \n']);
fprintf(fid,['  end fix1; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['  begin fix1;   \n']);
fprintf(fid,['    body1 shaft last; \n']);
fprintf(fid,['    body2 hub3 1; \n']);
fprintf(fid,['  end fix1; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['  begin bearing2;  Bearing where rotation is controlled from external dll\n']);
fprintf(fid,['    name  pitch1 ; \n']);
fprintf(fid,['    body1 hub1 last; \n']);
fprintf(fid,['    body2 blade1 1; \n']);
fprintf(fid,['    bearing_vector 2 0.0 0.0 -1.0; \n']);
fprintf(fid,['  end bearing2; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['  begin bearing2;  Bearing where rotation is controlled from external dll\n']);
fprintf(fid,['    name  pitch2 ; \n']);
fprintf(fid,['    body1 hub2 last; \n']);
fprintf(fid,['    body2 blade2 1; \n']);
fprintf(fid,['    bearing_vector 2 0.0 0.0 -1.0; \n']);
fprintf(fid,['  end bearing2; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['  begin bearing2;  Bearing where rotation is controlled from external dll\n']);
fprintf(fid,['    name  pitch3 ; \n']);
fprintf(fid,['    body1 hub3 last; \n']);
fprintf(fid,['    body2 blade3 1; \n']);
fprintf(fid,['    bearing_vector 2 0.0 0.0 -1.0; \n']);
fprintf(fid,['  end bearing2; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['  end constraint; \n']);
fprintf(fid,['; \n']);

end