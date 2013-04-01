function RunLin4GUI(handles)

%make fast command and then execute from working directory
fastexe=get(handles.edit_FastExecutable,'string');
fastcmd=['"',fastexe,'" gui_FastConfig.fst >gui_fastlin_out.txt'];

HHWind=get(handles.edit_HHWind,'string');
NomPitch=get(handles.edit_StartPitch,'string');
NomTorque=get(handles.edit_StartTorque,'string');
RotorSpeed=get(handles.edit_RotorSpeed,'string');
disp('Running Linearization:');
disp(['   Wind Speed:  ',HHWind]);
disp(['   Rotor Speed: ',RotorSpeed]);
disp(['   Nom. Pitch:  ',NomPitch]);
disp(['   Nom. Torque: ',NomTorque]);

set(handles.text_EndPitch,'string','running','ForegroundColor','g');
set(handles.text_EndTorque,'string','running','ForegroundColor','g');

%disable as much of the gui as possible
names=fieldnames(handles);
for index=1:length(names)
    if ~strcmp(names{index},'text_EndPitch')
        if ~strcmp(names{index},'text_EndTorque')
            if ~strcmp(names{index},'figure1')
                %disable object if possible
                handle=getfield(handles,names{index});
                if isfield(get(handle),'Enable')
                    set(handle,'Enable','off');
                end;
            end;
        end
    end
end

system(fastcmd);

%re-enable the gui
for index=1:length(names)
    if ~strcmp(names{index},'text_EndPitch')
        if ~strcmp(names{index},'text_EndTorque')
            if ~strcmp(names{index},'figure1')
                %enable object if possible
                handle=getfield(handles,names{index});
                if isfield(get(handle),'Enable')
                    set(handle,'Enable','on');
                end;
            end;
        end
    end
end

%check to see if linearization was successfull
fid=fopen('gui_fastlin_out.txt'); %fast cmd line outs are piped to this file

if fid>0
    LinConsoleOut_text=fscanf(fid,'%c');
    fclose(fid);
    
    test=regexp(LinConsoleOut_text,'Beginning iteration','once');
    if isempty(test)
        test=regexp(LinConsoleOut_text,'Linearizing FAST model about initial conditions','once');
        if (~isempty(test))
            success=1;
            LinResults=get_Lin4GUI(FASTLinName);
            set(handles.text_EndPitch,'string','IC Lin','ForegroundColor','g');
            set(handles.text_EndTorque,'string','IC Lin','ForegroundColor','g');
			
		    FastData=getappdata(handles.figure1,'FastData');
			FastData.EndPitch_text='IC Lin';
			FastData.EndTorque_text='IC Lin';
        else
            success=0;
            set(handles.text_EndPitch,'string','Failed','ForegroundColor','r');
            set(handles.text_EndTorque,'string','Failed','ForegroundColor','r');
        end;
    else
        test=regexp(LinConsoleOut_text,'does not appear to converge','once');
        if isempty(test)
            success=1;
            [trimPitch,trimGenTq]=get_FastLinPitchTorque(LinConsoleOut_text);
            LinResults=get_Lin4GUI('gui_FastConfig.lin');
            LinResults.trimPitch=trimPitch;
            LinResults.trimGenTq=trimGenTq;
			trimPitch=sprintf('%2.1f',trimPitch);
			trimGenTq=sprintf('%5.1f',trimGenTq);
            set(handles.text_EndPitch,'string',trimPitch,'ForegroundColor','g');
            set(handles.text_EndTorque,'string',trimGenTq,'ForegroundColor','g');
			
		    FastData=getappdata(handles.figure1,'FastData');
			FastData.EndPitch_text=trimPitch;
			FastData.EndTorque_text=trimGenTq;
        else
            success=0;
            set(handles.text_EndPitch,'string','Failed','ForegroundColor','r');
            set(handles.text_EndTorque,'string','Failed','ForegroundColor','r');
        end;
    end;
else %something went very wrong
    disp('System was not able to run:');
    disp(get(handles.edit_FastExecutable,'string'));
    set(handles.text_EndPitch,'string','not run','ForegroundColor','k');
    set(handles.text_EndTorque,'string','not run','ForegroundColor','k');
    success=0;
end;
if success==1
	LinResults.HHWind=FastData.HHWind;
	LinResults.RotorSpeed=FastData.RotorSpeed;
	FastData.LinConsoleOut_text=LinConsoleOut_text;
	fid=fopen('gui_FastConfig.lin'); %fast cmd line outs are piped to this file
	if fid>0
		LinResults_text=fscanf(fid,'%c');
		fclose(fid);
		FastData.LinResults_text=LinResults_text;
	else
		disp('Could not read gui_FastConfig.lin into FastData')
	end;
	W=['_W',sprintf('%2.1f',LinResults.HHWind)];
	W=strrep(W,'.','p');
	Rspd=['_Rspd',sprintf('%2.1f',LinResults.RotorSpeed)];
	Rspd=strrep(Rspd,'.','p');
	P=['_P',sprintf('%2.1f',LinResults.trimPitch)];
	P=strrep(P,'.','p');
	Tq=['_Tq',sprintf('%4.0f',LinResults.trimGenTq)];
	savename=['gui_LinResults',W,Rspd,P,Tq];
    save(savename,'LinResults','FastData');
	setappdata(handles.figure1,'FastData',FastData);
end;
