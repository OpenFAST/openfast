function varargout = FastGUI(varargin)
% FASTGUI MATLAB code for FastGUI.fig
%      FASTGUI, by itself, creates a new FASTGUI or raises the existing
%      singleton*.
%
%      H = FASTGUI returns the handle to a new FASTGUI or the handle to
%      the existing singleton*.
%
%      FASTGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FASTGUI.M with the given input arguments.
%
%      FASTGUI('Property','Value',...) creates a new FASTGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FastGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FastGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FastGUI

% Last Modified by GUIDE v2.5 17-Jan-2013 15:46:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FastGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @FastGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FastGUI is made visible.
function FastGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FastGUI (see VARARGIN)

% Choose default command line output for FastGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

%Try loading in gui data file from current directory
try
    load('FastGUIdata');
catch merr
    disp('Could not load FastGUIdata.mat from current directory');
    disp(' ');
    disp('Try loading:')
    default_guidata='Y:\Wind\Public\Projects\Projects A-E\Controls\Software\MATLAB\SimulationToolbox\FastGUI\FastGUIdata.mat';
    disp(default_guidata);
    try
        load(default_guidata);
    catch merr
        disp(' ');
        disp('Could not load:');
        disp(default_guidata);
    end;
end;
if exist('FastData','var') %then load above worked
    setappdata(handles.figure1,'FastData',FastData);  %initialize empty struct
    set_guidata(handles);
else
    setappdata(handles.figure1,'FastData',struct);  %initialize empty struct
end;


% --- Outputs from this function are returned to the command line.
function varargout = FastGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_WorkingDir_Callback(hObject, eventdata, handles)
% hObject    handle to edit_WorkingDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guidata(handles);


% --- Executes during object creation, after setting all properties.
function edit_WorkingDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_WorkingDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_WorkingDir.
function pushbutton_WorkingDir_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_WorkingDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
path=uigetdir;
if ischar(path)
    set(handles.edit_WorkingDir,'string',path);
    get_guidata(handles);
end;



function edit_FastTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FastTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guidata(handles);


% --- Executes during object creation, after setting all properties.
function edit_FastTemplate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FastTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_FastTemplate.
function pushbutton_FastTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_FastTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path]=uigetfile('*.*');
if ischar(path)
    set(handles.edit_FastTemplate,'string',[path,file]);
    get_guidata(handles);
end;


% --- Executes on button press in pushbutton_LoadTemplates.
function pushbutton_LoadTemplates_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LoadTemplates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
template=get(handles.edit_FastTemplate,'string');
if exist(template,'file')
    try
        FastParams=Fast2Matlab(template);
    catch merr
        disp('Fast2Matlab could not load:');
        disp(template);
        disp(' ');
        disp('Error from Fast2Matlab:')
        disp(merr.message);
    end;
    if exist('FastParams','var')
        fid=fopen(template,'r');
        fast_text=fscanf(fid,'%c');
        fclose(fid);
        %update app data
        FastData=getappdata(handles.figure1,'FastData');
        FastData.FastTemplateText=fast_text;
        FastData.FastParams=FastParams;
        setappdata(handles.figure1,'FastData',FastData);
        set_guidata(handles);
    end;
else
        disp('Cannot find Fast Template:');
        disp(template);
end;

template=get(handles.edit_AeroDynTemplate,'string');
if exist(template,'file')
    try
        AeroDynParams=AeroDyn2Matlab(template);
    catch merr
        disp('Could not load:');
        disp(template);
        disp(' ');
        disp('Error from AeroDyn2Matlab:')
        disp(merr.message);
    end;
    if exist('AeroDynParams','var')
        fid=fopen(template,'r');
        fast_text=fscanf(fid,'%c');
        fclose(fid);
        %update app data
        FastData=getappdata(handles.figure1,'FastData');
        FastData.AeroDynTemplateText=fast_text;
        FastData.AeroDynParams=AeroDynParams;
        setappdata(handles.figure1,'FastData',FastData);
        set_guidata(handles);
    end;
else
        disp('Cannot find Aerodyn Template:');
        disp(template);
end;

template=get(handles.edit_LinearTemplate,'string');
if exist(template,'file')
    try
        LinearParams=Linear2Matlab(template);
    catch merr
        disp('Could not load:');
        disp(template);
        disp(' ');
        disp('Error from Linear2Matlab:')
        disp(merr.message);
    end;
    if exist('LinearParams','var')
        fid=fopen(template,'r');
        fast_text=fscanf(fid,'%c');
        fclose(fid);
        %update app data
        FastData=getappdata(handles.figure1,'FastData');
        FastData.LinearTemplateText=fast_text;
        FastData.LinearParams=LinearParams;
		
		%set linearization input description
		LinearInputDesc={};
		NInputs = GetFastPar(FastData.LinearParams,'NInputs');
		if NInputs>0
			Controls = GetFastPar(FastData.LinearParams,'CntrlInpt');
			for index=1:length(Controls)
				switch Controls(index)
					case 1
						LinearInputDesc=[LinearInputDesc;'Yaw'];
					case 2
						LinearInputDesc=[LinearInputDesc;'YawRate'];
					case 3
						LinearInputDesc=[LinearInputDesc;'GenTq'];
					case 4
						LinearInputDesc=[LinearInputDesc;'PitchC'];
					case 5
						LinearInputDesc=[LinearInputDesc;'Pitch1'];
					case 6
						LinearInputDesc=[LinearInputDesc;'Pitch2'];
					case 7
						LinearInputDesc=[LinearInputDesc;'Pitch3'];
				end %switch
			end;
		end;
		NDisturbs = GetFastPar(FastData.LinearParams,'NDisturbs');
		if NDisturbs>0
			disturbs = GetFastPar(FastData.LinearParams,'Disturbnc');
			for index=1:length(disturbs)
				switch disturbs(index)
					case 1
						LinearInputDesc=[LinearInputDesc;'HorWind'];
					case 2
						LinearInputDesc=[LinearInputDesc;'HWdirect'];
					case 3
						LinearInputDesc=[LinearInputDesc;'VerWind'];
					case 4
						LinearInputDesc=[LinearInputDesc;'HorShear'];
					case 5
						LinearInputDesc=[LinearInputDesc;'PwrVshear'];
					case 6
						LinearInputDesc=[LinearInputDesc;'LinVshear'];
					case 7
						LinearInputDesc=[LinearInputDesc;'HHHgust'];
				end %switch
			end;
		end;
		FastData.LinearInputDesc=LinearInputDesc;
		
        setappdata(handles.figure1,'FastData',FastData);
        set_guidata(handles);
    end;
else
        disp('Cannot find Linearization Template:');
        disp(template);
end;

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guidata(handles);
FastData=getappdata(handles.figure1,'FastData');
curdir=pwd;
try
    cd(FastData.WorkingDir);
catch merr
    disp('Could not save gui data in working dir');
    disp(['Saving in: ',curdir]);
    cd(curdir);
    save('FastGUIdata','FastData');
    return;
end;
save('FastGUIdata','FastData');
cd(curdir);

% --------------------------------------------------------------------
function menu_Load_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path]=uigetfile('.mat');
if isstr(path) %filename is good, load data and update gui
    gui_file=[path,file];
    load(gui_file);
    setappdata(handles.figure1,'FastData',FastData);
    set_guidata(handles);
else
    disp('Nothing loaded.');
end;


% --------------------------------------------------------------------
function menu_File_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_FastExecutable.
function pushbutton_FastExecutable_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_FastExecutable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path]=uigetfile('*.*');
if ischar(path)
    set(handles.edit_FastExecutable,'string',[path,file]);
    get_guidata(handles);
end;

function edit_FastExecutable_Callback(hObject, eventdata, handles)
% hObject    handle to edit_FastExecutable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guidata(handles);

% --- Executes during object creation, after setting all properties.
function edit_FastExecutable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_FastExecutable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_FastDOFs.
function pushbutton_FastDOFs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_FastDOFs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FastDOFs(handles);


function edit_TMax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_TMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_TMax as text
%        str2double(get(hObject,'String')) returns contents of edit_TMax as a double


% --- Executes during object creation, after setting all properties.
function edit_TMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_TMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_DT_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_DT as text
%        str2double(get(hObject,'String')) returns contents of edit_DT as a double


% --- Executes during object creation, after setting all properties.
function edit_DT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_SimulationMode.
function uipanel_SimulationMode_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_SimulationMode 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
FastData=getappdata(handles.figure1,'FastData');

if strcmp(get(eventdata.NewValue,'string'),'TimeMarch')
    FastData.FastParams=SetFastPar(FastData.FastParams,'AnalMode',1);
else
    FastData.FastParams=SetFastPar(FastData.FastParams,'AnalMode',2);
end;
setappdata(handles.figure1,'FastData',FastData);


% --- Executes when selected object is changed in uipanel_PitchControl.
function uipanel_PitchControl_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_PitchControl 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
FastData=getappdata(handles.figure1,'FastData');
if strcmp(get(eventdata.NewValue,'string'),'PitchNone')
    FastData.FastParams=SetFastPar(FastData.FastParams,'PCMode',0);
elseif strcmp(get(eventdata.NewValue,'string'),'PitchCode')
    FastData.FastParams=SetFastPar(FastData.FastParams,'PCMode',0);
else
    FastData.FastParams=SetFastPar(FastData.FastParams,'PCMode',0);
end;
setappdata(handles.figure1,'FastData',FastData);


% --- Executes when selected object is changed in uipanel_YawControl.
function uipanel_YawControl_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_YawControl 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
FastData=getappdata(handles.figure1,'FastData');
if strcmp(get(eventdata.NewValue,'string'),'YawNone')
    FastData.FastParams=SetFastPar(FastData.FastParams,'YCMode',0);
elseif strcmp(get(eventdata.NewValue,'string'),'YawCode')
    FastData.FastParams=SetFastPar(FastData.FastParams,'YCMode',1);
else
    FastData.FastParams=SetFastPar(FastData.FastParams,'YCMode',2);
end;
setappdata(handles.figure1,'FastData',FastData);


% --- Executes when selected object is changed in uipanel_GeneratorControl.
function uipanel_GeneratorControl_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_GeneratorControl 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
FastData=getappdata(handles.figure1,'FastData');
if strcmp(get(eventdata.NewValue,'string'),'GenNone')
    FastData.FastParams=SetFastPar(FastData.FastParams,'VSContrl',0);
elseif strcmp(get(eventdata.NewValue,'string'),'GenSimple')
    FastData.FastParams=SetFastPar(FastData.FastParams,'VSContrl',1);
elseif strcmp(get(eventdata.NewValue,'string'),'GenCode')
    FastData.FastParams=SetFastPar(FastData.FastParams,'VSContrl',2);
else
    FastData.FastParams=SetFastPar(FastData.FastParams,'VSContrl',3);
end;
setappdata(handles.figure1,'FastData',FastData);


function edit_AeroDynTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_AeroDynTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guidata(handles);


% --- Executes during object creation, after setting all properties.
function edit_AeroDynTemplate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_AeroDynTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_AeroDynTemplate.
function pushbutton_AeroDynTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_AeroDynTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path]=uigetfile('*.*');
if ischar(path)
    set(handles.edit_AeroDynTemplate,'string',[path,file]);
    get_guidata(handles);
end;





function edit_StartPitch_Callback(hObject, eventdata, handles)
% hObject    handle to edit_StartPitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_StartPitch as text
%        str2double(get(hObject,'String')) returns contents of edit_StartPitch as a double


% --- Executes during object creation, after setting all properties.
function edit_StartPitch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_StartPitch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_StartTorque_Callback(hObject, eventdata, handles)
% hObject    handle to edit_StartTorque (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_StartTorque as text
%        str2double(get(hObject,'String')) returns contents of edit_StartTorque as a double


% --- Executes during object creation, after setting all properties.
function edit_StartTorque_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_StartTorque (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_HHWind_Callback(hObject, eventdata, handles)
% hObject    handle to edit_HHWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_HHWind as text
%        str2double(get(hObject,'String')) returns contents of edit_HHWind as a double


% --- Executes during object creation, after setting all properties.
function edit_HHWind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_HHWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_RotorSpeed_Callback(hObject, eventdata, handles)
% hObject    handle to edit_RotorSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_RotorSpeed as text
%        str2double(get(hObject,'String')) returns contents of edit_RotorSpeed as a double


% --- Executes during object creation, after setting all properties.
function edit_RotorSpeed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_RotorSpeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_LinearTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_LinearTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guidata(handles);


% --- Executes during object creation, after setting all properties.
function edit_LinearTemplate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_LinearTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_LinearTemplate.
function pushbutton_LinearTemplate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LinearTemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file,path]=uigetfile('*.*');
if ischar(path)
    set(handles.edit_LinearTemplate,'string',[path,file]);
    get_guidata(handles);
end;

function get_guidata(handles)
FastData=getappdata(handles.figure1,'FastData');

%any change in gui settings makes linearization status invalid:
set(handles.text_EndPitch,'string','not run','ForegroundColor','k');
set(handles.text_EndTorque,'string','not run','ForegroundColor','k');
if isfield(FastData,'LinResults_text')
	FastData=rmfield(FastData,'LinResults_text');
end;
if isfield(FastData,'LinConsoleOut_text')
	FastData=rmfield(FastData,'LinConsoleOut_text');
end;
if isfield(FastData,'EndPitch_text')
	FastData=rmfield(FastData,'EndPitch_text');
end;
if isfield(FastData,'EndTorque_text')
	FastData=rmfield(FastData,'EndTorque_text');
end;

%store editable text
FastData.FastExecutable=get(handles.edit_FastExecutable,'string');
FastData.WorkingDir=get(handles.edit_WorkingDir,'string');
FastData.FastTemplate=get(handles.edit_FastTemplate,'string');
FastData.AeroDynTemplate=get(handles.edit_AeroDynTemplate,'string');
FastData.LinearTemplate=get(handles.edit_LinearTemplate,'string');

%store editable numbers
TMax=str2num(get(handles.edit_TMax,'String'));
if isempty(TMax)
    set(handles.edit_TMax,'String','0')
    FastData.FastParams=SetFastPar(FastData.FastParams,'TMax',0);
else
    FastData.FastParams=SetFastPar(FastData.FastParams,'TMax',TMax);
end;

DT=str2num(get(handles.edit_DT,'String'));
if isempty(DT)
    set(handles.edit_DT,'String','0.025')
    FastData.FastParams=SetFastPar(FastData.FastParams,'DT',0.025);
else
    FastData.FastParams=SetFastPar(FastData.FastParams,'DT',DT);
end;

StartPitch=str2num(get(handles.edit_StartPitch,'String'));
if isempty(StartPitch)
    set(handles.edit_StartPitch,'String','10')
    FastData.FastParams=SetFastPar(FastData.FastParams,'BlPitch(1)',10);
    FastData.FastParams=SetFastPar(FastData.FastParams,'BlPitch(2)',10);
    FastData.FastParams=SetFastPar(FastData.FastParams,'BlPitch(3)',10);
	FastData.StartPitch=10;
else
    FastData.FastParams=SetFastPar(FastData.FastParams,'BlPitch(1)',StartPitch);
    FastData.FastParams=SetFastPar(FastData.FastParams,'BlPitch(2)',StartPitch);
    FastData.FastParams=SetFastPar(FastData.FastParams,'BlPitch(3)',StartPitch);
	FastData.StartPitch=StartPitch;
end;

FastData.StartTorque=str2num(get(handles.edit_StartTorque,'String'));
if isempty(FastData.StartTorque)
    set(handles.edit_StartTorque,'String','1000')
    FastData.StartTorque=1000;
end;

DispTol=get(handles.edit_DispTol,'String');
if iscell(DispTol)
    DispTol=DispTol{:};
end;
DispTol=str2num(DispTol);
if isempty(DispTol)
    set(handles.edit_DispTol,'String','0.01')
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'DispTol',0.01);
else
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'DispTol',DispTol);
end;

VelTol=get(handles.edit_VelTol,'String');
if iscell(VelTol)
    VelTol=VelTol{:};
end;
VelTol=str2num(VelTol);
if isempty(VelTol)
    set(handles.edit_VelTol,'String','0.01')
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'VelTol',0.01);
else
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'VelTol',VelTol);
end;

Azimuths=get(handles.edit_Azimuths,'String');
if iscell(Azimuths)
    Azimuths=Azimuths{:};
end;
Azimuths=str2num(Azimuths);
if isempty(Azimuths)
    set(handles.edit_Azimuths,'String','12')
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'NAzimStep',12);
else
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'NAzimStep',Azimuths);
end;

HHWind=get(handles.edit_HHWind,'String');
if iscell(HHWind)
    HHWind=HHWind{:};
end;

FastData.HHWind=str2num(HHWind);
if isempty(FastData.HHWind)
    set(handles.edit_StartTorque,'String','18')
    FastData.HHWind=18;
end;

RotorSpeed=str2num(get(handles.edit_RotorSpeed,'String'));
if isempty(RotorSpeed)
    set(handles.edit_RotorSpeed,'String','30')
    FastData.FastParams=SetFastPar(FastData.FastParams,'RotSpeed',30);
	FastData.RotorSpeed=30;
else
    FastData.FastParams=SetFastPar(FastData.FastParams,'RotSpeed',RotorSpeed);
	FastData.RotorSpeed=RotorSpeed;
end;

%check boxes
if get(handles.checkbox_CalcSteady,'value')==1
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'CalcStdy','True');
else
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'CalcStdy','False');
end;

setappdata(handles.figure1,'FastData',FastData);


function set_guidata(handles)
FastData=getappdata(handles.figure1,'FastData'); 
%set non-editable text
if isfield(FastData,'EndPitch_text')
	set(handles.text_EndPitch,'string',FastData.EndPitch_text,'ForegroundColor','k');
else
	set(handles.text_EndPitch,'string','not run','ForegroundColor','k');
end;
if isfield(FastData,'EndTorque_text')
	set(handles.text_EndTorque,'string',FastData.EndTorque_text,'ForegroundColor','k');
else
	set(handles.text_EndTorque,'string','not run','ForegroundColor','k');
end;

%set editable text
set(handles.edit_FastExecutable,'string',FastData.FastExecutable);
set(handles.edit_WorkingDir,'string',FastData.WorkingDir);
set(handles.edit_FastTemplate,'string',FastData.FastTemplate);
set(handles.edit_AeroDynTemplate,'string',FastData.AeroDynTemplate);
set(handles.edit_LinearTemplate,'string',FastData.LinearTemplate);

%set radio buttons
if GetFastPar(FastData.FastParams,'AnalMode')==1
    set(handles.radiobutton_TimeMarch,'Value',1);
else
    set(handles.radiobutton_PeriodicLin,'Value',1);
end;

if GetFastPar(FastData.FastParams,'YCMode')==0
    set(handles.radiobutton_YawNone,'Value',1);
elseif GetFastPar(FastData.FastParams,'YCMode')==1
    set(handles.radiobutton_YawCode,'Value',1);
else
    set(handles.radiobutton_YawSimulink,'Value',1);
end;

if GetFastPar(FastData.FastParams,'PCMode')==0
    set(handles.radiobutton_PitchNone,'Value',1);
elseif GetFastPar(FastData.FastParams,'PCMode')==1
    set(handles.radiobutton_PitchCode,'Value',1);
else
    set(handles.radiobutton_PitchSimulink,'Value',1);
end;

if GetFastPar(FastData.FastParams,'VSContrl')==0
    set(handles.radiobutton_GenNone,'Value',1);
elseif GetFastPar(FastData.FastParams,'VSContrl')==1
    set(handles.radiobutton_GenSimple,'Value',1);
elseif GetFastPar(FastData.FastParams,'VSContrl')==2
    set(handles.radiobutton_GenCode,'Value',1);
else
    set(handles.radiobutton_GenSimulink,'Value',1);
end;

if GetFastPar(FastData.LinearParams,'TrimCase')==1
    set(handles.radiobutton_TrimYaw,'Value',1);
elseif GetFastPar(FastData.LinearParams,'TrimCase')==2
    set(handles.radiobutton_TrimTorque,'Value',1);
else
    set(handles.radiobutton_TrimPitch,'Value',1);
end;

if GetFastPar(FastData.LinearParams,'MdlOrder')==1
    set(handles.radiobutton_ModOrder1,'Value',1);
else
    set(handles.radiobutton_ModOrder2,'Value',1);
end;

%check boxes
if strcmpi(GetFastPar(FastData.LinearParams,'CalcStdy'),'True')
    set(handles.checkbox_CalcSteady,'Value',1);
else
    set(handles.checkbox_CalcSteady,'Value',0);
end;

%set editable numbers
TMax=sprintf('%3.1f',GetFastPar(FastData.FastParams,'TMax'));
set(handles.edit_TMax,'String',TMax);

DT=sprintf('%1.2d',GetFastPar(FastData.FastParams,'DT'));
set(handles.edit_DT,'String',DT);

StartPitch=GetFastPar(FastData.FastParams,'BlPitch(1)');
StartPitch=StartPitch+GetFastPar(FastData.FastParams,'BlPitch(2)');
StartPitch=StartPitch+GetFastPar(FastData.FastParams,'BlPitch(3)');
StartPitch=sprintf('%2.1f',StartPitch/3);
set(handles.edit_StartPitch,'String',StartPitch);

StartTorque=GetFastPar(FastData.FastParams,'VS_RtTq');
StartTorque=sprintf('%5.2f',StartTorque);
set(handles.edit_StartTorque,'String',StartTorque);

DispTol=GetFastPar(FastData.LinearParams,'DispTol');
set(handles.edit_DispTol,'String',DispTol);

VelTol=GetFastPar(FastData.LinearParams,'VelTol');
set(handles.edit_VelTol,'String',VelTol);

Azimuths=GetFastPar(FastData.LinearParams,'NAzimStep');
set(handles.edit_Azimuths,'String',Azimuths);

HHWind=sprintf('%3.2f',FastData.HHWind);
set(handles.edit_HHWind,'String',HHWind);

RotorSpeed=GetFastPar(FastData.FastParams,'RotSpeed');
RotorSpeed=sprintf('%3.2f',RotorSpeed);
set(handles.edit_RotorSpeed,'String',RotorSpeed);

% --- Executes on button press in pushbutton_LinearInputs.
function pushbutton_LinearInputs_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_LinearInputs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
LinearInputs(handles);


% --- Executes on button press in pushbutton_WriteTemplates.
function pushbutton_WriteTemplates_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_WriteTemplates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%make sure all gui data items have been placed in FastData
get_guidata(handles);

workdir=get(handles.edit_WorkingDir,'string');
if ~exist(workdir,'file')
    error('Not writing files; cannot find working directory.')
end;

currentdir=pwd;
cd(workdir);

FastData=getappdata(handles.figure1,'FastData');

%write out fast config file
if ~isfield(FastData,'FastTemplateText')
    error('Not writing fast file-- source template text is not loaded.')
else
    %point to gui aerodyn file to be written subsequently
    FastData.FastParams=SetFastPar(FastData.FastParams,'ADFile','"gui_AeroDynConfig.ipt"');
    %point to gui linearization file to be written subsequently
    FastData.FastParams=SetFastPar(FastData.FastParams,'LinFile','"gui_LinearConfig.ipt"');
    
    fid=fopen('gui_tempfile','w');
    if fid>=0
        fprintf(fid,'%s',FastData.FastTemplateText);
        fclose(fid);
        Matlab2FAST(FastData.FastParams,'gui_tempfile','gui_FastConfig.fst');
    else
        error('Not writing fast file-- could not open temp file.')
    end;
end;    

%write out AeroDyn config file
if ~isfield(FastData,'AeroDynTemplateText')
    error('Not writing AeroDyn file-- source template text is not loaded.')
else
    %point to gui wind file to be written following writing of aerodyn
    FastData.AeroDynParams=SetFastPar(FastData.AeroDynParams,'WindFile','gui_wind.wnd');
    
    fid=fopen('gui_tempfile','w');
    if fid>=0
        fprintf(fid,'%s',FastData.AeroDynTemplateText);
        fclose(fid);
        Matlab2AD(FastData.AeroDynParams,'gui_tempfile','gui_AeroDynConfig.ipt')
    else
        error('Not writing AeroDyn file-- could not open temp file.')
    end;
end;

%write out lin config file
if ~isfield(FastData,'LinearTemplateText')
    error('Not writing linear config file-- source template text is not loaded.')
else
    
    fid=fopen('gui_tempfile','w');
    if fid>=0
        fprintf(fid,'%s',FastData.LinearTemplateText);
        fclose(fid);
        Matlab2FASTLin(FastData.LinearParams,'gui_tempfile','gui_LinearConfig.ipt')
    else
        error('Not writing Linear Config file-- could not open temp file.')
    end;
end;

%write out wind file
wind=FastData.HHWind;
wind=[wind;wind];  %Constant wind (until gui is enhanced)
try
    writeUniformWindFile(wind,'gui_wind.wnd');
catch merr
    error(['Could not write wind file: ',merr.message]);
end;

cd(currentdir);


% --- Executes on button press in checkbox_CalcSteady.
function checkbox_CalcSteady_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_CalcSteady (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guidata(handles);



function edit_VelTol_Callback(hObject, eventdata, handles)
% hObject    handle to edit_VelTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guidata(handles);


% --- Executes during object creation, after setting all properties.
function edit_VelTol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_VelTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_DispTol_Callback(hObject, eventdata, handles)
% hObject    handle to edit_DispTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guidata(handles);


% --- Executes during object creation, after setting all properties.
function edit_DispTol_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_DispTol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_RunLinearization.
function pushbutton_RunLinearization_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_RunLinearization (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%set simulation type to periodic linearization
FastData=getappdata(handles.figure1,'FastData');
curdir=pwd;
try
    cd(FastData.WorkingDir);
catch merr
    disp('Could not cd to:');
    disp(FastData.WorkingDir);
    disp(' ');
    disp('Not running FAST!');
    return;
end;

set(handles.radiobutton_PeriodicLin,'Value',1);
% % set(handles.radiobutton_YawNone,'Value',1);
% % set(handles.radiobutton_PitchNone,'Value',1);
% % set(handles.radiobutton_GenNone,'Value',1);
get_guidata(handles);

%write out config files
try
    pushbutton_WriteTemplates_Callback(hObject, eventdata, handles);
catch merr
    disp('Could not write config files');
    disp(merr.message);
    return;
end    

RunLin4GUI(handles);
cd(curdir);

function edit_Azimuths_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Azimuths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guidata(handles);


% --- Executes during object creation, after setting all properties.
function edit_Azimuths_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Azimuths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_ModelOrder.
function uipanel_ModelOrder_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_ModelOrder 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
FastData=getappdata(handles.figure1,'FastData');
if strcmp(get(eventdata.NewValue,'string'),'TimeMarch')
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'MdlOrder',1);
else
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'MdlOrder',2);
end;
setappdata(handles.figure1,'FastData',FastData);


% --- Executes when selected object is changed in uipanel_TrimType.
function uipanel_TrimType_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_TrimType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
FastData=getappdata(handles.figure1,'FastData');
if strcmp(get(eventdata.NewValue,'string'),'TrimYaw')
    FastData.FastParams=SetFastPar(FastData.LinearParams,'TrimCase',1);
elseif strcmp(get(eventdata.NewValue,'string'),'TrimTorque')
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'TrimCase',2);
else
    FastData.LinearParams=SetFastPar(FastData.LinearParams,'TrimCase',3);
end;
setappdata(handles.figure1,'FastData',FastData);
