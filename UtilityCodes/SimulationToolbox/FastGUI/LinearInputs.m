function varargout = LinearInputs(varargin)
% LINEARINPUTS MATLAB code for LinearInputs.fig
%      LINEARINPUTS, by itself, creates a new LINEARINPUTS or raises the existing
%      singleton*.
%
%      H = LINEARINPUTS returns the handle to a new LINEARINPUTS or the handle to
%      the existing singleton*.
%
%      LINEARINPUTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LINEARINPUTS.M with the given input arguments.
%
%      LINEARINPUTS('Property','Value',...) creates a new LINEARINPUTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LinearInputs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LinearInputs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LinearInputs

% Last Modified by GUIDE v2.5 18-Jan-2013 12:54:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LinearInputs_OpeningFcn, ...
                   'gui_OutputFcn',  @LinearInputs_OutputFcn, ...
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


% --- Executes just before LinearInputs is made visible.
function LinearInputs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LinearInputs (see VARARGIN)

% Choose default command line output for LinearInputs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

parentHandles=varargin{1};
setappdata(handles.figure1,'parentHandles',parentHandles);
set_guiLinearInputs(handles);


% --- Outputs from this function are returned to the command line.
function varargout = LinearInputs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox_Yaw.
function checkbox_Yaw_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Yaw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);

% --- Executes on button press in checkbox_YawRate.
function checkbox_YawRate_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_YawRate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_GenTq.
function checkbox_GenTq_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_GenTq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_HorWind.
function checkbox_HorWind_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_HorWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_VerWind.
function checkbox_VerWind_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_VerWind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_HorShear.
function checkbox_HorShear_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_HorShear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_PwrVshear.
function checkbox_PwrVshear_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PwrVshear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_HHHgust.
function checkbox_HHHgust_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_HHHgust (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_LinVshear.
function checkbox_LinVshear_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_LinVshear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_PitchC.
function checkbox_PitchC_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_PitchC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_Pitch1.
function checkbox_Pitch1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Pitch1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_Pitch2.
function checkbox_Pitch2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Pitch2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_Pitch3.
function checkbox_Pitch3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Pitch3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes on button press in checkbox_HWdirect.
function checkbox_HWdirect_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_HWdirect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiLinearInputs(handles);
