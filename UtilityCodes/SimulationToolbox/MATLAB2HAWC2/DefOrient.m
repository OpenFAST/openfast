% Defining and writting orientations of main bodies in htc file
%
% In: FastPar       -   FAST input parameters extracted using Fast2Matlab.m
%     fid           -   File identifier of Hawc2 htc file
%
% Knud A. Kragh

function DefOrient(FastPar,fid)

fprintf(fid,['  begin orientation; \n']);
fprintf(fid,['    begin base; \n']);
fprintf(fid,['      body   tower;; \n']);
fprintf(fid,['      inipos        0.0 0.0 0.0 ;         initial position of node 1 \n']);
fprintf(fid,['      body_eulerang 0.0 0.0 0.0;         initial position of node 1 \n']);
fprintf(fid,['    end base; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['    begin relative; \n']);
fprintf(fid,['      body1  tower last; \n']);
fprintf(fid,['      body2  shaft 1; \n']);
fprintf(fid,['      body2_eulerang 90.0 0.0 0.0;  \n']);
fprintf(fid,['      body2_eulerang ' num2str(GetFastPar(FastPar,'ShftTilt')) ' 0.0 0.0;     tilt \n']);
fprintf(fid,['      body2_ini_rotvec_d1 0.0 0.0 -1.0 ' num2str(GetFastPar(FastPar,'RotSpeed')/60*2*pi) '; body initial rotation velocity x,y,z,angle velocity[rad/s]  (body 2 coordinates) \n']);
fprintf(fid,['    end relative; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['    begin relative; \n']);
fprintf(fid,['      body1  shaft last; \n']);
fprintf(fid,['      body2  hub1 1; \n']);
fprintf(fid,['      body2_eulerang -90.0 0.0 0.0;  \n']);
fprintf(fid,['      body2_eulerang ' num2str(-GetFastPar(FastPar,'PreCone(1)')) ' 0.0 0.0;     tilt \n']);
fprintf(fid,['    end relative; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['    begin relative; \n']);
fprintf(fid,['      body1  shaft last; \n']);
fprintf(fid,['      body2  hub2 1; \n']);
fprintf(fid,['      body2_eulerang -90.0 0.0 0.0;  \n']);
fprintf(fid,['      body2_eulerang 0.0 -120.0 0.0;  \n']);
fprintf(fid,['      body2_eulerang ' num2str(GetFastPar(FastPar,'PreCone(2)')) ' 0.0 0.0;     tilt \n']);
fprintf(fid,['    end relative; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['    begin relative; \n']);
fprintf(fid,['      body1  shaft last; \n']);
fprintf(fid,['      body2  hub3 1; \n']);
fprintf(fid,['      body2_eulerang -90.0 0.0 0.0;  \n']);
fprintf(fid,['      body2_eulerang 0.0 120.0 0.0;  \n']);
fprintf(fid,['      body2_eulerang ' num2str(GetFastPar(FastPar,'PreCone(2)')) ' 0.0 0.0;     tilt \n']);
fprintf(fid,['    end relative; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['    begin relative; \n']);
fprintf(fid,['      body1  hub1 last; \n']);
fprintf(fid,['      body2  blade1 1; \n']);
fprintf(fid,['      body2_eulerang 0.0 0.0 0.0;  \n']);
fprintf(fid,['    end relative; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['    begin relative; \n']);
fprintf(fid,['      body1  hub2 last; \n']);
fprintf(fid,['      body2  blade2 1; \n']);
fprintf(fid,['      body2_eulerang 0.0 0.0 0.0;  \n']);
fprintf(fid,['    end relative; \n']);
fprintf(fid,['; \n']);

fprintf(fid,['    begin relative; \n']);
fprintf(fid,['      body1  hub3 last; \n']);
fprintf(fid,['      body2  blade3 1; \n']);
fprintf(fid,['      body2_eulerang 0.0 0.0 0.0;  \n']);
fprintf(fid,['    end relative; \n']);
fprintf(fid,['; \n']);

fprintf(fid,[  '  end orientation; \n']);
fprintf(fid,['; \n']);
end
