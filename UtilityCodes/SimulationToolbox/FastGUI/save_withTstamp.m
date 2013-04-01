function save_file=save_withTstamp(file_name,data_list)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%				Jason Laks
%		University of Colorado at Boulder
%				Nov., 2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_time=fix(clock);
save_time=num2str(save_time(2:end-1));

save_file=[file_name,regexprep(save_time,'\s*','_')];
save_cmd=['save ',save_file,' ',data_list]; %,' save_file'];
try
    evalin('caller',save_cmd);
catch
    err=lasterror;
    disp(err.message)
end;

 