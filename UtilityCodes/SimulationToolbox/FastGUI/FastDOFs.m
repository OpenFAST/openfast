function varargout = FastDOFs(varargin)
% FASTDOFS MATLAB code for FastDOFs.fig
%      FASTDOFS, by itself, creates a new FASTDOFS or raises the existing
%      singleton*.
%
%      H = FASTDOFS returns the handle to a new FASTDOFS or the handle to
%      the existing singleton*.
%
%      FASTDOFS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FASTDOFS.M with the given input arguments.
%
%      FASTDOFS('Property','Value',...) creates a new FASTDOFS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FastDOFs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FastDOFs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FastDOFs

% Last Modified by GUIDE v2.5 15-Jan-2013 10:33:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FastDOFs_OpeningFcn, ...
                   'gui_OutputFcn',  @FastDOFs_OutputFcn, ...
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


% --- Executes just before FastDOFs is made visible.
function FastDOFs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FastDOFs (see VARARGIN)

% Choose default command line output for FastDOFs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

parentHandles=varargin{1};
setappdata(handles.figure1,'parentHandles',parentHandles);
set_guiFastDOFs(handles);

% --- Outputs from this function are returned to the command line.
function varargout = FastDOFs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in checkbox_Flap1.
function checkbox_Flap1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Flap1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes on button press in checkbox_Flap2.
function checkbox_Flap2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Flap2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes on button press in checkbox_Edge.
function checkbox_Edge_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Edge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes on button press in checkbox_Teeter.
function checkbox_Teeter_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Teeter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes on button press in checkbox_Generator.
function checkbox_Generator_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Generator (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes on button press in checkbox_DriveTrain.
function checkbox_DriveTrain_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_DriveTrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes on button press in checkbox_Yaw.
function checkbox_Yaw_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Yaw (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes on button press in checkbox_ForeAft1.
function checkbox_ForeAft1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ForeAft1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes on button press in checkbox_ForeAft2.
function checkbox_ForeAft2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ForeAft2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes on button press in checkbox_Side2Side1.
function checkbox_Side2Side1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Side2Side1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes on button press in checkbox_Side2Side2.
function checkbox_Side2Side2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Side2Side2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes on button press in checkbox_Aero.
function checkbox_Aero_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Aero (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);


% --- Executes on button press in checkbox_Noise.
function checkbox_Noise_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_Noise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
get_guiFastDOFs(handles);
