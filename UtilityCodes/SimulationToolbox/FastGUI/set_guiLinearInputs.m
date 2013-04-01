function set_guiLinearInputs(handles)
parentHandles=getappdata(handles.figure1,'parentHandles');
FastData=getappdata(parentHandles.figure1,'FastData');

%first clear and then set according to fast lin data
set(handles.checkbox_Yaw,'value',0);
set(handles.checkbox_YawRate,'value',0);
set(handles.checkbox_GenTq,'value',0);
set(handles.checkbox_PitchC,'value',0);
set(handles.checkbox_Pitch1,'value',0);
set(handles.checkbox_Pitch2,'value',0);
set(handles.checkbox_Pitch3,'value',0);

NInputs = GetFastPar(FastData.LinearParams,'NInputs');

%set check boxes if NInputs>0
if NInputs>0
    Controls = GetFastPar(FastData.LinearParams,'CntrlInpt');
    for index=1:length(Controls)
        switch Controls(index)
            case 1
                set(handles.checkbox_Yaw,'value',1);
            case 2
                set(handles.checkbox_YawRate,'value',1);
            case 3
                set(handles.checkbox_GenTq,'value',1);
            case 4
                set(handles.checkbox_PitchC,'value',1);
            case 5
                set(handles.checkbox_Pitch1,'value',1);
            case 6
                set(handles.checkbox_Pitch2,'value',1);
            case 7
                set(handles.checkbox_Pitch3,'value',1);
        end %switch
    end;
end;

%now do the same for the disturbance inputs
set(handles.checkbox_HorWind,'value',0);
set(handles.checkbox_HWdirect,'value',0);
set(handles.checkbox_VerWind,'value',0);
set(handles.checkbox_HHHgust,'value',0);
set(handles.checkbox_HorShear,'value',0);
set(handles.checkbox_LinVshear,'value',0);
set(handles.checkbox_PwrVshear,'value',0);

NDisturbs = GetFastPar(FastData.LinearParams,'NDisturbs');

%set check boxes if NDisturbs>0
if NDisturbs>0
    disturbs = GetFastPar(FastData.LinearParams,'Disturbnc');
    for index=1:length(disturbs)
        switch disturbs(index)
            case 1
                set(handles.checkbox_HorWind,'value',1);
            case 2
                set(handles.checkbox_HWdirect,'value',1);
            case 3
                set(handles.checkbox_VerWind,'value',1);
            case 4
                set(handles.checkbox_HorShear,'value',1);
            case 5
                set(handles.checkbox_PwrVshear,'value',1);
            case 6
                set(handles.checkbox_LinVshear,'value',1);
            case 7
                set(handles.checkbox_HHHgust,'value',1);
        end %switch
    end;
end;

