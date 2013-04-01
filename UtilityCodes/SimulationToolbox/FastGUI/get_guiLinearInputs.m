function get_guiLinearInputs(handles)
parentHandles=getappdata(handles.figure1,'parentHandles');
FastData=getappdata(parentHandles.figure1,'FastData');

Controls=[];
LinearInputDesc={};
if get(handles.checkbox_Yaw,'value')==1
    Controls=[Controls,1];
	LinearInputDesc=[LinearInputDesc;'Yaw'];
end;
if get(handles.checkbox_YawRate,'value')==1
    Controls=[Controls,2];
	LinearInputDesc=[LinearInputDesc;'YawRate'];
end;
if get(handles.checkbox_GenTq,'value')==1
    Controls=[Controls,3];
	LinearInputDesc=[LinearInputDesc;'GenTq'];
end;
if get(handles.checkbox_PitchC,'value')==1
    Controls=[Controls,4];
	LinearInputDesc=[LinearInputDesc;'PC'];
end;
if get(handles.checkbox_Pitch1,'value')==1
    Controls=[Controls,5];
	LinearInputDesc=[LinearInputDesc;'Pitch1'];
end;
if get(handles.checkbox_Pitch2,'value')==1
    Controls=[Controls,6];
	LinearInputDesc=[LinearInputDesc;'Pitch2'];
end;
if get(handles.checkbox_Pitch3,'value')==1
    Controls=[Controls,7];
	LinearInputDesc=[LinearInputDesc;'Pitch3'];
end;

FastData.LinearParams=SetFastPar(FastData.LinearParams,'NInputs',length(Controls));
if isempty(Controls)
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'CntrlInpt',1);
else
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'CntrlInpt',Controls);
end

Disturbs=[];
if get(handles.checkbox_HorWind,'value')==1
    Disturbs=[Disturbs,1];
	LinearInputDesc=[LinearInputDesc;'HorWind'];
end;
if get(handles.checkbox_HWdirect,'value')==1
    Disturbs=[Disturbs,2];
	LinearInputDesc=[LinearInputDesc;'HWdirect'];
end;
if get(handles.checkbox_VerWind,'value')==1
    Disturbs=[Disturbs,3];
	LinearInputDesc=[LinearInputDesc;'VerWind'];
end;
if get(handles.checkbox_HorShear,'value')==1
    Disturbs=[Disturbs,4];
	LinearInputDesc=[LinearInputDesc;'HorShear'];
end;
if get(handles.checkbox_PwrVshear,'value')==1
    Disturbs=[Disturbs,5];
	LinearInputDesc=[LinearInputDesc;'PwrVshear'];
end;
if get(handles.checkbox_LinVshear,'value')==1
    Disturbs=[Disturbs,6];
	LinearInputDesc=[LinearInputDesc;'LinVshear'];
end;
if get(handles.checkbox_HHHgust,'value')==1
    Disturbs=[Disturbs,7];
	LinearInputDesc=[LinearInputDesc;'HHHgust'];
end;

FastData.LinearInputDesc=LinearInputDesc;
FastData.LinearParams=SetFastPar(FastData.LinearParams,'NDisturbs',length(Disturbs));
if isempty(Disturbs)
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'Disturbnc',1);
else
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'Disturbnc',Disturbs);
end
set(parentHandles.text_EndPitch,'string','not run','ForegroundColor','k');
set(parentHandles.text_EndTorque,'string','not run','ForegroundColor','k');
setappdata(parentHandles.figure1,'FastData',FastData);

