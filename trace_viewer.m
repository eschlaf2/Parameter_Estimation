function varargout = trace_viewer(varargin)
% TRACE_VIEWER MATLAB code for trace_viewer.fig
%      TRACE_VIEWER, by itself, creates a new TRACE_VIEWER or raises the existing
%      singleton*.
%
%      H = TRACE_VIEWER returns the handle to a new TRACE_VIEWER or the handle to
%      the existing singleton*.
%
%      TRACE_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TRACE_VIEWER.M with the given input arguments.
%
%      TRACE_VIEWER('Property','Value',...) creates a new TRACE_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before trace_viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to trace_viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help trace_viewer

% Last Modified by GUIDE v2.5 15-Sep-2017 19:07:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @trace_viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @trace_viewer_OutputFcn, ...
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


% --- Executes just before trace_viewer is made visible.
function trace_viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to trace_viewer (see VARARGIN)

% Choose default command line output for trace_viewer
handles.output = hObject;
handles.data = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes trace_viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = trace_viewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in series_next.
function series_next_Callback(hObject, eventdata, handles)
% hObject    handle to series_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
series = get_series(handles.series_edit.String) + 1;
series(series > size(handles.data,1)) = size(handles.data,1);
% display(series)
handles.series_edit.String = num2str(series);
update_button_Callback(handles.update_button, eventdata, handles);


% --- Executes on button press in series_prev.
function series_prev_Callback(hObject, eventdata, handles)
% hObject    handle to series_prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
seriesNum = get_series(handles.series_edit.String) - 1;
seriesNum(seriesNum < 1) = 1;
handles.series_edit.String = num2str(seriesNum);
update_button_Callback(handles.update_button, eventdata, handles);


function series_edit_Callback(hObject, eventdata, handles)
% hObject    handle to series_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of series_edit as text
%        str2double(get(hObject,'String')) returns contents of series_edit as a double
seriesNum = str2double(hObject.String);
seriesNum(seriesNum < 1) = 1;
seriesNum(seriesNum > size(handles.data, 1)) = size(handles.data, 1);
hObject.String = num2str(seriesNum);

update_button_Callback(handles.update_button, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function series_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to series_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% hObject.Value = str2double(hObject.String);


% --- Executes on button press in time_prev.
function time_prev_Callback(hObject, eventdata, handles)
% hObject    handle to time_prev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rate = handles.sampleRate_edit.Value;
interval = handles.interval_edit.Value;
step_size = floor((interval * rate) / 2);
display_object = handles.startTime_edit;
display_object.Value = display_object.Value - (step_size / rate);
display_object.String = num2str(display_object.Value);
update_button_Callback(handles.update_button, eventdata, handles);


% --- Executes on button press in time_next.
function time_next_Callback(hObject, eventdata, handles)
% hObject    handle to time_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rate = handles.sampleRate_edit.Value;
interval = handles.interval_edit.Value;
step_size = floor((interval * rate) / 2);
% display_object = handles.startTime_edit;
handles.startTime_edit.Value = handles.startTime_edit.Value + ...
    (step_size / rate);
handles.startTime_edit.String = num2str(handles.startTime_edit.Value);
update_button_Callback(handles.update_button, eventdata, handles);



function interval_edit_Callback(hObject, eventdata, handles)
% hObject    handle to interval_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of interval_edit as text
%        str2double(get(hObject,'String')) returns contents of interval_edit as a double
hObject.Value = str2double(hObject.String);
update_button_Callback(handles.update_button, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function interval_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interval_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
hObject.Value = str2double(hObject.String);



function startTime_edit_Callback(hObject, eventdata, handles)
% hObject    handle to startTime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of startTime_edit as text
%        str2double(get(hObject,'String')) returns contents of startTime_edit as a double
hObject.Value = str2double(hObject.String);
update_button_Callback(handles.update_button, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function startTime_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to startTime_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
hObject.Value = str2double(hObject.String);



function sampleRate_edit_Callback(hObject, eventdata, handles)
% hObject    handle to sampleRate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sampleRate_edit as text
%        str2double(get(hObject,'String')) returns contents of sampleRate_edit as a double
% hObject.Value = str2double(hObject.String);
old_rate = hObject.Value;
new_rate = eval(hObject.String);
new_startTime = convert_units(handles.startTime_edit.Value, ...
    old_rate, new_rate);
new_interval = convert_units(handles.interval_edit.Value, ...
    old_rate, new_rate);

hObject.Value = new_rate;
handles.startTime_edit.Value = new_startTime;
handles.startTime_edit.String = num2str(new_startTime);
handles.interval_edit.Value = new_interval;
handles.interval_edit.String = num2str(new_interval);
update_button_Callback(handles.update_button, eventdata, handles);




% --- Executes during object creation, after setting all properties.
function sampleRate_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sampleRate_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
hObject.Value = str2double(hObject.String);



function units_edit_Callback(hObject, eventdata, handles)
% hObject    handle to units_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of units_edit as text
%        str2double(get(hObject,'String')) returns contents of units_edit as a double
update_button_Callback(handles.update_button, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function units_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to units_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- 
function data_edit_Callback(hObject, eventdata, handles)
% hObject    handle to data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of data_edit as text
%        str2double(get(hObject,'String')) returns contents of data_edit as a double
update_button_Callback(handles.update_button, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function data_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in update_button.
function update_button_Callback(hObject, eventdata, handles)
% hObject    handle to update_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla;
seriesNum = get_series(handles.series_edit.String);
data = evalin('base', handles.data_edit.String);
% display(data)
rate = handles.sampleRate_edit.Value;
desired_interval = rate * handles.interval_edit.Value;
samples_in_interval = round(min(desired_interval, size(data, 2)));
desired_first_sample = get_first_sample(handles.sampleRate_edit.Value, ...
    handles.startTime_edit.Value);
first_sample = round(min(max(desired_first_sample, 1), ...
    size(data, 2) - samples_in_interval + 1));

last_sample = first_sample + samples_in_interval - 1;
% display({series, first_sample, last_sample})

if samples_in_interval ~= desired_interval
    handles.interval_edit.String = num2str(samples_in_interval/rate);
    handles.interval_edit.Value = samples_in_interval;
end
if first_sample ~= desired_first_sample
    handles.startTime_edit.Value = ...
        get_start_time(first_sample, handles.sampleRate_edit.Value);
    handles.startTime_edit.String = num2str(handles.startTime_edit.Value);
end

X = repmat((first_sample:last_sample) / rate, numel(seriesNum));
plot(X', data(seriesNum, (first_sample:last_sample))');
% plot_labels = str2double(get(gca, 'xticklabels')) - 1 + first_sample;
% plot_labels = plot_labels / handles.sampleRate_edit.Value;
% plot_labels = plot_labels + (handles.startTime_edit.Value - plot_labels(1));
% set(gca, 'xticklabels', num2str(plot_labels))
xlabel(['Time (', handles.units_edit.String, ')'])
axis('tight')
ylim([min(data(:)), max(data(:))])

handles.data = data;
guidata(hObject, handles);


% --- internal function
function first_sample = get_first_sample(rate, start_time)
first_sample = rate * start_time + 1;


% --- internal function
function start_time = get_start_time(first_sample, rate)
start_time = (first_sample - 1) / rate;


% --- internal function
function series = get_series(input_string)
try
    series = eval(input_string);
catch
    series = eval(['[', input_string, ']']);
end


% --- Internal
function new_units = convert_units(old_units, old_rate, new_rate)
new_units = (old_units * old_rate)/new_rate;
