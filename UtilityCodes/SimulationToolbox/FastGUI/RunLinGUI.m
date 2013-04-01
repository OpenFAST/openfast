function varargout = RunLinGUI(varargin)
% RUNLINGUI M-file for RunLinGUI.fig
%      RUNLINGUI, by itself, creates a new RUNLINGUI or raises the existing
%      singleton*.
%
%      H = RUNLINGUI returns the handle to a new RUNLINGUI or the handle to
%      the existing singleton*.
%
%      RUNLINGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RUNLINGUI.M with the given input arguments.
%
%      RUNLINGUI('Property','Value',...) creates a new RUNLINGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before RunLinGUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to RunLinGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help RunLinGUI

% Last Modified by GUIDE v2.5 11-Mar-2009 17:32:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @RunLinGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @RunLinGUI_OutputFcn, ...
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


% --- Executes just before RunLinGUI is made visible.
function RunLinGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to RunLinGUI (see VARARGIN)

% Choose default command line output for RunLinGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

parenthandles=varargin{1};
setappdata(handles.figure1,'parenthandles',parenthandles);

NominalPitch=getappdata(parenthandles.figure1,'NominalPitch');
set(handles.Pitch_EditText,'string',num2str(NominalPitch));

NominalTorque=getappdata(parenthandles.figure1,'NominalTorque');
set(handles.Torque_EditText,'string',num2str(NominalTorque));

NominalSpeed=getappdata(parenthandles.figure1,'NominalSpeed');
set(handles.Speed_EditText,'string',num2str(NominalSpeed));

NominalWind=getappdata(parenthandles.figure1,'NominalWind');
set(handles.Wind_EditText,'string',num2str(NominalWind));

selectTorqueTrim=getappdata(parenthandles.figure1,'TorqueTrim');
set(handles.Torque_radio,'value',selectTorqueTrim);
set(handles.Pitch_radio,'value',~selectTorqueTrim);
set(handles.trimPitchResults_text,'string','not run');
set(handles.trimTorqueResults_text,'string','not run');
set(handles.WptsRun_text,'string','not run');
set(handles.PptsRun_text,'string','not run');
set(handles.LSSVptsRun_text,'string','not run');

RangePmin=getappdata(parenthandles.figure1,'RangePmin');
RangePmax=getappdata(parenthandles.figure1,'RangePmax');
RangePpts=getappdata(parenthandles.figure1,'RangePpts');
RangeWmin=getappdata(parenthandles.figure1,'RangeWmin');
RangeWmax=getappdata(parenthandles.figure1,'RangeWmax');
RangeWpts=getappdata(parenthandles.figure1,'RangeWpts');
RangeLSSVmin=getappdata(parenthandles.figure1,'RangeLSSVmin');
RangeLSSVmax=getappdata(parenthandles.figure1,'RangeLSSVmax');
RangeLSSVpts=getappdata(parenthandles.figure1,'RangeLSSVpts');
TPSopt=getappdata(parenthandles.figure1,'TPSopt');
TipRadius=getappdata(parenthandles.figure1,'TipRadius');

set(handles.Pmin_EditText,'string',num2str(RangePmin));
set(handles.Pmax_EditText,'string',num2str(RangePmax));
set(handles.Ppts_EditText,'string',num2str(RangePpts));
set(handles.Wmin_EditText,'string',num2str(RangeWmin));
set(handles.Wmax_EditText,'string',num2str(RangeWmax));
set(handles.Wpts_EditText,'string',num2str(RangeWpts));
set(handles.LSSVmin_EditText,'string',num2str(RangeLSSVmin));
set(handles.LSSVmax_EditText,'string',num2str(RangeLSSVmax));
set(handles.LSSVpts_EditText,'string',num2str(RangeLSSVpts));

set(handles.Pmin_EditText,'Enable','off');
set(handles.Pmax_EditText,'Enable','off');
set(handles.Ppts_EditText,'Enable','off');
set(handles.Wmin_EditText,'Enable','off');
set(handles.Wmax_EditText,'Enable','off');
set(handles.Wpts_EditText,'Enable','off');
set(handles.LSSVmin_EditText,'Enable','off');
set(handles.LSSVmax_EditText,'Enable','off');
set(handles.LSSVpts_EditText,'Enable','off');
set(handles.TPSopt_EditText,'string',num2str(TPSopt));
set(handles.TipRadius_EditText,'string',num2str(TipRadius));

gridRange=getappdata(parenthandles.figure1,'gridRange');
set(handles.gridrange_radio,'value',gridRange);
set(handles.TPSrange_radio,'value',~gridRange);

set(handles.Rangecheckbox,'Value',0);
% UIWAIT makes RunLinGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = RunLinGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function Pitch_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to Pitch_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NominalPitch=str2double(get(hObject,'String'));
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'NominalPitch',NominalPitch);


% --- Executes during object creation, after setting all properties.
function Pitch_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pitch_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Speed_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to Speed_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NominalSpeed=str2double(get(hObject,'String'));
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'NominalSpeed',NominalSpeed);


% --- Executes during object creation, after setting all properties.
function Speed_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Speed_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Wind_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to Wind_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NominalWind=str2double(get(hObject,'String'));
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'NominalWind',NominalWind);


% --- Executes during object creation, after setting all properties.
function Wind_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Wind_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in StartLinButton.
function StartLinButton_Callback(hObject, eventdata, handles)
% hObject    handle to StartLinButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
success=1;

if ~isappdata(parenthandles.figure1,'FastContentsLin_text')
    disp('FastContentsLin_text not initialized: select and load fast file(s).');
    success=0;
end
if ~isappdata(parenthandles.figure1,'WindLinTemplate_text')
    disp('WindLinTemplate_text not initialized: select and load fast file(s).');
    success=0;
end
if ~isappdata(parenthandles.figure1,'LinFileContentsLin_text')
    disp('LinFileContentsLin_text not initialized: select and load fast file(s).');
    success=0;
end
if ~isappdata(parenthandles.figure1,'ADFileContentsLin_text')
    disp('ADFileContentsLin_text not initialized: select and load fast file(s).');
    success=0;
end

if ~success
    set(handles.trimPitchResults_text,'string','fail');
    set(handles.trimTorqueResults_text,'string','fail');
    set(handles.WptsRun_text,'string','fail');
    set(handles.PptsRun_text,'string','fail');
    set(handles.LSSVptsRun_text,'string','fail');
    return;
end;

NominalPitch=getappdata(parenthandles.figure1,'NominalPitch');
NominalWind=getappdata(parenthandles.figure1,'NominalWind');
NominalTorque=getappdata(parenthandles.figure1,'NominalTorque');
NominalSpeed=getappdata(parenthandles.figure1,'NominalSpeed');

selectTorqueTrim=getappdata(parenthandles.figure1,'TorqueTrim');
gridRange=getappdata(parenthandles.figure1,'gridRange');
TPSopt=getappdata(parenthandles.figure1,'TPSopt');
TipRadius=getappdata(parenthandles.figure1,'TipRadius');

if get(handles.Rangecheckbox,'value')
    if gridRange==1
        RangePmin=getappdata(parenthandles.figure1,'RangePmin');
        RangePmax=getappdata(parenthandles.figure1,'RangePmax');
        RangePpts=getappdata(parenthandles.figure1,'RangePpts');
        RangeWmin=getappdata(parenthandles.figure1,'RangeWmin');
        RangeWmax=getappdata(parenthandles.figure1,'RangeWmax');
        RangeWpts=getappdata(parenthandles.figure1,'RangeWpts');
        RangeLSSVmin=getappdata(parenthandles.figure1,'RangeLSSVmin');
        RangeLSSVmax=getappdata(parenthandles.figure1,'RangeLSSVmax');
        RangeLSSVpts=getappdata(parenthandles.figure1,'RangeLSSVpts');
    else
        RangeWmin=getappdata(parenthandles.figure1,'RangeWmin');
        RangeWmax=getappdata(parenthandles.figure1,'RangeWmax');
        RangeWpts=getappdata(parenthandles.figure1,'RangeWpts');
%         TipSpeedRatio=kron(1./WindsII,TipSpeed);
        RangeLSSVmin=TPSopt*RangeWmin/2/pi/TipRadius*60;
        RangeLSSVmax=TPSopt*RangeWmax/2/pi/TipRadius*60;
        RangeLSSVpts=RangeWpts;
        setappdata(parenthandles.figure1,'RangeLSSVmin',RangeLSSVmin);
        setappdata(parenthandles.figure1,'RangeLSSVmax',RangeLSSVmax);
        setappdata(parenthandles.figure1,'RangeLSSVpts',RangeLSSVpts);
        set(handles.LSSVmin_EditText,'string',(RangeLSSVmin));
        set(handles.LSSVmax_EditText,'string',num2str(RangeLSSVmax));
        set(handles.LSSVpts_EditText,'string',num2str(RangeLSSVpts));

        RangePmin=NominalPitch;
        RangePmax=NominalPitch;
        RangePpts=1;
    end;
else
    RangePmin=NominalPitch;
    RangePmax=NominalPitch;
    RangePpts=1;
    RangeWmin=NominalWind;
    RangeWmax=NominalWind;
    RangeWpts=1;
    RangeLSSVmin=NominalSpeed;
    RangeLSSVmax=NominalSpeed;
    RangeLSSVpts=1;
end;

Pitchs=linspace(RangePmin,RangePmax,RangePpts);
Winds=linspace(RangeWmin,RangeWmax,RangeWpts);
LSSVs=linspace(RangeLSSVmin,RangeLSSVmax,RangeLSSVpts);

if gridRange==1 %doing grid
    N_LSSV=length(LSSVs);
else %LSSV is determined by wind speed and TPSopt
    N_LSSV=1;
end;

%change directory to spec'd FAST data dir
present_dir=pwd;
data_dir=getappdata(parenthandles.figure1,'DataDir');
cd(data_dir);

TRIMPITCHS=zeros(length(Winds),length(Pitchs),length(LSSVs));
TRIMTORQUES=zeros(length(Winds),length(Pitchs),length(LSSVs));
SUCCESSFLAGS=zeros(length(Winds),length(Pitchs),length(LSSVs));
OutputNames=getappdata(parenthandles.figure1,'OutListLin_text');
lin_text=getappdata(parenthandles.figure1,'LinFileContentsLin_text');
[Ninputs,CntrlInputs,Ndisturbs,DistInputs]=getFAST_LinFileInputs(lin_text);

AVG_BDSYS=ss(zeros(length(OutputNames),Ninputs+Ndisturbs,RangeWpts,RangePpts,N_LSSV));
AVG_BDSYS.inputname={CntrlInputs{:},DistInputs{:}};
for w_index=1:length(Winds)
    for p_index=1:length(Pitchs)
        for omega_index=1:N_LSSV
            
            if gridRange==1 %doing grid
                LSSV=LSSVs(omega_index);
                setappdata(parenthandles.figure1,'NominalSpeed',LSSVs(omega_index));
                set(handles.LSSVptsRun_text,'string',sprintf('(%i)%4.2f',omega_index,LSSVs(omega_index)));
            else  %LSSV is determined by wind speed and TPSopt via w_index
                LSSV=LSSVs(w_index);
                setappdata(parenthandles.figure1,'NominalSpeed',LSSVs(w_index));
                set(handles.LSSVptsRun_text,'string',sprintf('(%i)%4.2f',w_index,LSSVs(w_index)));
            end;

            %set op points used by RunLin4GUI
            setappdata(parenthandles.figure1,'NominalPitch',Pitchs(p_index));
            setappdata(parenthandles.figure1,'NominalWind',Winds(w_index));
            
            set(handles.PptsRun_text,'string',sprintf('(%i)%4.2f',p_index,Pitchs(p_index)));
            set(handles.WptsRun_text,'string',sprintf('(%i)%4.2f',w_index,Winds(w_index)));
            
            [success]=RunLin4GUI(handles);
            if ~success
                set(handles.trimPitchResults_text,'string','fail');
                set(handles.trimTorqueResults_text,'string','fail');
                set(handles.WptsRun_text,'string','fail');
                set(handles.PptsRun_text,'string','fail');
                set(handles.LSSVptsRun_text,'string','fail');
            else
                %get final generator trim settings
                fid=fopen('fst_out.txt');
                if fid>0
                    fast_text=fscanf(fid,'%c');
                    fclose(fid);
                    
                    test=regexp(fast_text,'Beginning iteration','once');
                    if isempty(test)
                        test=regexp(fast_text,'Linearizing FAST model about initial conditions','once');
                        if (~isempty(test))
                            success=1;
                            trimPitch=NaN;
                            trimGenTq=NaN;
                            setappdata(parenthandles.figure1,'trimPitch',trimPitch);
                            setappdata(parenthandles.figure1,'trimGenTq',trimGenTq);
                            set(handles.trimPitchResults_text,'string','no trim');
                            set(handles.trimTorqueResults_text,'string','no trim');
                        else
                            success=0;
                        end;
                    else
                        test=regexp(fast_text,'does not appear to converge','once');
                        if isempty(test)
                            success=1;
                            [trimPitch,trimGenTq]=get_FastLinPitchTorque(fast_text);
                            setappdata(parenthandles.figure1,'trimPitch',trimPitch);
                            setappdata(parenthandles.figure1,'trimGenTq',trimGenTq);
                            set(handles.trimPitchResults_text,'string',num2str(trimPitch));
                            set(handles.trimTorqueResults_text,'string',num2str(trimGenTq));
                        else
                            success=0;
                        end;
                    end;
                end;
                if ~success
                    set(handles.trimPitchResults_text,'string','fail');
                    set(handles.trimTorqueResults_text,'string','fail');
                    set(handles.WptsRun_text,'string','fail');
                    set(handles.PptsRun_text,'string','fail');
                    set(handles.LSSVptsRun_text,'string','fail');
                end;
            end;
            
            if ~success  %will put empty mats in results cell arrays
                
                display('Linearization failed.')
                avg_bdsys=ss(zeros(length(OutputNames),Ninputs+Ndisturbs));
                trimPitch=[];
                trimGenTq=[];
                AVG_BDSYS(:,:,w_index,p_index,omega_index)=avg_bdsys;
                TRIMPITCHS(w_index,p_index,omega_index)=NaN;
                TRIMTORQUES(w_index,p_index,omega_index)=NaN;
            else
                display('Linearization Successful.')
                
                %store avg system and pitch/wind settings
                AvgAMat=mean(getappdata(parenthandles.figure1,'AMat'),3);
                AvgBMat=mean(getappdata(parenthandles.figure1,'BMat'),3);
                AvgBdMat=mean(getappdata(parenthandles.figure1,'BdMat'),3);
                AvgCMat=mean(getappdata(parenthandles.figure1,'CMat'),3);
                AvgDMat=mean(getappdata(parenthandles.figure1,'DMat'),3);
                AvgDdMat=mean(getappdata(parenthandles.figure1,'DdMat'),3);
                
                avg_bdsys=ss(AvgAMat,[AvgBMat,AvgBdMat],AvgCMat,[AvgDMat,AvgDdMat]);
                avg_bdsys.OutputName=getappdata(parenthandles.figure1,'OutListLin_text');
                
                AVG_BDSYS(:,:,w_index,p_index,omega_index)=avg_bdsys;
                TRIMPITCHS(w_index,p_index,omega_index)=trimPitch;
                TRIMTORQUES(w_index,p_index,omega_index)=trimGenTq;
            end;
        end;
        
    end;
end;

%store results into gui data
setappdata(parenthandles.figure1,'AVG_BDSYS',AVG_BDSYS);
setappdata(parenthandles.figure1,'TRIMPITCHS',TRIMPITCHS);
setappdata(parenthandles.figure1,'TRIMTORQUES',TRIMTORQUES);
setappdata(parenthandles.figure1,'SUCCESSFLAGS',SUCCESSFLAGS);

%reset nominal settings
setappdata(parenthandles.figure1,'NominalPitch',NominalPitch);
setappdata(parenthandles.figure1,'NominalWind',NominalWind);
setappdata(parenthandles.figure1,'NominalTorque',NominalTorque);

set(handles.PptsRun_text,'string','done');
set(handles.WptsRun_text,'string','done');
set(handles.LSSVptsRun_text,'string','done');

%save results
save_file=fix(clock);
save_file=num2str(save_file(2:end-1));
save_file=['FASTGUILinearization',regexprep(save_file,'\s*','_')];
setappdata(parenthandles.figure1,'save_file',save_file);

Fgui_data=getappdata(parenthandles.figure1);
Fgui_data.Listeners = [];

assignin('base','Fgui_data',Fgui_data);

save(save_file,'Fgui_data');
cd(present_dir);

function Torque_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to Torque_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
NominalTorque=str2double(get(hObject,'String'));
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'NominalTorque',NominalTorque);


% --- Executes during object creation, after setting all properties.
function Torque_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Torque_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
NominalPitch=str2double(get(handles.Pitch_EditText,'String'));
setappdata(parenthandles.figure1,'NominalPitch',NominalPitch);
NominalTorque=str2double(get(handles.Torque_EditText,'String'));
setappdata(parenthandles.figure1,'NominalTorque',NominalTorque);
NominalSpeed=str2double(get(handles.Speed_EditText,'String'));
setappdata(parenthandles.figure1,'NominalSpeed',NominalSpeed);
NominalWind=str2double(get(handles.Wind_EditText,'String'));
setappdata(parenthandles.figure1,'NominalWind',NominalWind);


% --- Executes on button press in Pitch_radio.
function Pitch_radio_Callback(hObject, eventdata, handles)
% hObject    handle to Pitch_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'PitchTrim',1);
setappdata(parenthandles.figure1,'TorqueTrim',0);


% --- Executes on button press in Torque_radio.
function Torque_radio_Callback(hObject, eventdata, handles)
% hObject    handle to Torque_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'TorqueTrim',1);
setappdata(parenthandles.figure1,'PitchTrim',0);



function Pmin_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to Pmin_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'RangePmin',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function Pmin_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pmin_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Ppts_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to Ppts_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'RangePpts',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function Ppts_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ppts_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Wmin_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to Wmin_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'RangeWmin',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function Wmin_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Wmin_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Wpts_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to Wpts_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'RangeWpts',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function Wpts_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Wpts_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Rangecheckbox.
function Rangecheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to Rangecheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
do_range=get(handles.Rangecheckbox,'value');

if do_range
    parenthandles=getappdata(handles.figure1,'parenthandles');
    gridRange=getappdata(parenthandles.figure1,'gridRange');
    
    if gridRange==1
        set(handles.Pmin_EditText,'Enable','on');
        set(handles.Pmax_EditText,'Enable','on');
        set(handles.Ppts_EditText,'Enable','on');
        set(handles.Wmin_EditText,'Enable','on');
        set(handles.Wmax_EditText,'Enable','on');
        set(handles.Wpts_EditText,'Enable','on');
        set(handles.LSSVmin_EditText,'Enable','on');
        set(handles.LSSVmax_EditText,'Enable','on');
        set(handles.LSSVpts_EditText,'Enable','on');
    else
        set(handles.Pmin_EditText,'Enable','off');
        set(handles.Pmax_EditText,'Enable','off');
        set(handles.Ppts_EditText,'Enable','off');
        set(handles.Wmin_EditText,'Enable','on');
        set(handles.Wmax_EditText,'Enable','on');
        set(handles.Wpts_EditText,'Enable','on');
        set(handles.LSSVmin_EditText,'Enable','off');
        set(handles.LSSVmax_EditText,'Enable','off');
        set(handles.LSSVpts_EditText,'Enable','off');
    end;
else
    set(handles.Pmin_EditText,'Enable','off');
    set(handles.Pmax_EditText,'Enable','off');
    set(handles.Ppts_EditText,'Enable','off');
    set(handles.Wmin_EditText,'Enable','off');
    set(handles.Wmax_EditText,'Enable','off');
    set(handles.Wpts_EditText,'Enable','off');
    set(handles.LSSVmin_EditText,'Enable','off');
    set(handles.LSSVmax_EditText,'Enable','off');
    set(handles.LSSVpts_EditText,'Enable','off');
end;



function Pmax_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to Pmax_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'RangePmax',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function Pmax_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pmax_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Wmax_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to Wmax_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'RangeWmax',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function Wmax_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Wmax_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LSSVmin_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to LSSVmin_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'RangeLSSVmin',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function LSSVmin_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LSSVmin_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LSSVpts_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to LSSVpts_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'RangeLSSVpts',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function LSSVpts_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LSSVpts_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LSSVmax_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to LSSVmax_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'RangeLSSVmax',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function LSSVmax_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LSSVmax_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in gridrange_radio.
function gridrange_radio_Callback(hObject, eventdata, handles)
% hObject    handle to gridrange_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'gridRange',1);
setappdata(parenthandles.figure1,'TPSRange',0);
Rangecheckbox_Callback(hObject, eventdata, handles)

% --- Executes on button press in TPSrange_radio.
function TPSrange_radio_Callback(hObject, eventdata, handles)
% hObject    handle to TPSrange_radio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'gridRange',0);
setappdata(parenthandles.figure1,'TPSRange',1);
Rangecheckbox_Callback(hObject, eventdata, handles)





function TPSopt_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to TPSopt_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'TPSopt',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function TPSopt_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TPSopt_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function TipRadius_EditText_Callback(hObject, eventdata, handles)
% hObject    handle to TipRadius_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
parenthandles=getappdata(handles.figure1,'parenthandles');
setappdata(parenthandles.figure1,'TipRadius',str2double(get(hObject,'String')));


% --- Executes during object creation, after setting all properties.
function TipRadius_EditText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TipRadius_EditText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


