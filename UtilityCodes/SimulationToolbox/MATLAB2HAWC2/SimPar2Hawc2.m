% Function for writting simulation block in HAWC2 input file
%
% In:   T             -     Simulation time
%       dt            -     Time step
%       OutFileNames  -     Name of animation an logfile (Animation is out commented)
%       fid           -     File identifier for HAWC2 input file
%
% Knud A. Kragh

function SimPar2Hawc2(T,dt,OutFileNames,fid)

fprintf(fid,['begin Simulation; \n']);
fprintf(fid,['  time_stop '   num2str(T)  ' ; \n']);
fprintf(fid,['  solvertype   1 ;    (newmark) \n']);
fprintf(fid,['  on_no_convergence continue ; \n']);
fprintf(fid,['  logfile ./logfiles/' OutFileNames '.log ; \n']);
fprintf(fid,[';  animation ./animation/' OutFileNames '.dat; \n']);
fprintf(fid,['; \n']);
fprintf(fid,['  begin newmark; \n']);
fprintf(fid,['    deltat '    num2str(dt) ' ;  \n']);
fprintf(fid,['  end newmark; \n']);
fprintf(fid,['end Simulation; \n']);
fprintf(fid,['; \n']);


end