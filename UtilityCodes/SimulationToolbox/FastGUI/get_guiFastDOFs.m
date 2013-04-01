function get_guiFastDOFs(handles)
parentHandles=getappdata(handles.figure1,'parentHandles');
FastData=getappdata(parentHandles.figure1,'FastData');

index=find(strcmp('CompAero',FastData.FastParams.Label));
if get(handles.checkbox_Aero,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('GenDOF',FastData.FastParams.Label));
if get(handles.checkbox_Generator,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('DrTrDOF',FastData.FastParams.Label));
if get(handles.checkbox_DriveTrain,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('FlapDOF1',FastData.FastParams.Label));
if get(handles.checkbox_Flap1,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('FlapDOF2',FastData.FastParams.Label));
if get(handles.checkbox_Flap2,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('EdgeDOF',FastData.FastParams.Label));
if get(handles.checkbox_Edge,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('TwSSDOF1',FastData.FastParams.Label));
if get(handles.checkbox_Side2Side1,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('TwSSDOF2',FastData.FastParams.Label));
if get(handles.checkbox_Side2Side2,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('YawDOF',FastData.FastParams.Label));
if get(handles.checkbox_Yaw,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('TwFADOF1',FastData.FastParams.Label));
if get(handles.checkbox_ForeAft1,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('TwFADOF2',FastData.FastParams.Label));
if get(handles.checkbox_ForeAft2,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('TeetDOF',FastData.FastParams.Label));
if get(handles.checkbox_Teeter,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;

index=find(strcmp('CompNoise',FastData.FastParams.Label));
if get(handles.checkbox_Noise,'value')==1
    FastData.FastParams.Val{index}='True';
else
    FastData.FastParams.Val{index}='False';
end;
set(parentHandles.text_EndPitch,'string','not run','ForegroundColor','k');
set(parentHandles.text_EndTorque,'string','not run','ForegroundColor','k');
setappdata(parentHandles.figure1,'FastData',FastData);

