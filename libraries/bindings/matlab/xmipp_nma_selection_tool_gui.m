function varargout = xmipp_nma_selection_tool_gui(varargin)
% XMIPP_NMA_SELECTION_TOOL_GUI MATLAB code for xmipp_nma_selection_tool_gui.fig
%      XMIPP_NMA_SELECTION_TOOL_GUI, by itself, creates a new XMIPP_NMA_SELECTION_TOOL_GUI or raises the existing
%      singleton*.
%
%      H = XMIPP_NMA_SELECTION_TOOL_GUI returns the handle to a new XMIPP_NMA_SELECTION_TOOL_GUI or the handle to
%      the existing singleton*.
%
%      XMIPP_NMA_SELECTION_TOOL_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in XMIPP_NMA_SELECTION_TOOL_GUI.M with the given input arguments.
%
%      XMIPP_NMA_SELECTION_TOOL_GUI('Property','Value',...) creates a new XMIPP_NMA_SELECTION_TOOL_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before xmipp_nma_selection_tool_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to xmipp_nma_selection_tool_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help xmipp_nma_selection_tool_gui

% Last Modified by GUIDE v2.5 27-Aug-2013 18:36:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @xmipp_nma_selection_tool_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @xmipp_nma_selection_tool_gui_OutputFcn, ...
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


% --- Executes just before xmipp_nma_selection_tool_gui is made visible.
function xmipp_nma_selection_tool_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to xmipp_nma_selection_tool_gui (see VARARGIN)

% Choose default command line output for xmipp_nma_selection_tool_gui
handles.output = hObject;

handles.rundir=varargin{2};
set(handles.nmaRun,'String',['NMA directory: ' handles.rundir]);
handles.fnProjected=[handles.rundir '/extra/deformationsProjected.txt'];

% Open NMA results

[handles.images, handles.NMAdisplacements, handles.cost]= xmipp_nma_read_alignment(handles.rundir);
%handles.NMAdisplacements=rand(100,4);

%handles.cost=rand(size(handles.NMAdisplacements,1),1);
%for i=1:length(handles.cost)
%    handles.images{i}=['file ' int2str(i) '.xmp'];
%end

if exist(handles.fnProjected,'file') && 0
    handles.NMAdisplacementsProjected=load(handles.fnProjected);
else
    handles.NMAdisplacementsProjected=handles.NMAdisplacements;
end
handles.figHandle=figure();
updateListBox(hObject, handles);
set(handles.listRepresentation,'Value',[1 2])
guidata(hObject,handles)
updatePlot(handles)

% UIWAIT makes xmipp_nma_selection_tool_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = xmipp_nma_selection_tool_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function ndimensions_Callback(hObject, eventdata, handles)
% Do nothing

% --- Executes on selection change in popupProjection.
function popupProjection_Callback(hObject, eventdata, handles)
% hObject    handle to popupProjection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupProjection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupProjection


% --- Executes during object creation, after setting all properties.
function popupProjection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupProjection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonProjection.
function pushbuttonProjection_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonProjection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
projectionOptions=cellstr(get(handles.popupProjection,'String'));
selectedProjection=projectionOptions{get(handles.popupProjection,'Value')};
dout=get(handles.ndimensions,'String');
fnProjector=[handles.rundir '/extra/projector.txt'];

cmd=['xmipp_matrix_dimred -i ' handles.rundir '/extra/deformations.txt --din ' ...
    num2str(size(handles.NMAdisplacements,2)) ' --samples ' ...
    num2str(size(handles.NMAdisplacements,1)) ' -o ' ...
    handles.fnProjected ' --dout ' ...
    dout ' -m '];
if strcmp(selectedProjection,'None')==1
    cmd=['cp ' handles.rundir '/extra/deformations.txt ' handles.fnProjected];
elseif strcmp(selectedProjection,'Principal Component Analysis')==1
    cmd=[cmd 'PCA --saveMapping ' fnProjector];
elseif strcmp(selectedProjection,'Kernel Principal Component Analysis')==1
    cmd=[cmd 'kPCA'];
elseif strcmp(selectedProjection,'Probabilistic Principal Component Analysis')==1
    cmd=[cmd 'pPCA --saveMapping ' fnProjector];
elseif strcmp(selectedProjection,'Local Tangent Space Alignment')==1
    cmd=[cmd 'LTSA'];
elseif strcmp(selectedProjection,'Linear Local Tangent Space Alignment')==1
    cmd=[cmd 'LLTSA --saveMapping ' fnProjector];
elseif strcmp(selectedProjection,'Diffusion Map')==1
    cmd=[cmd 'DM'];
elseif strcmp(selectedProjection,'Linearity Preserving Projection')==1
    cmd=[cmd 'LPP --saveMapping ' fnProjector];
elseif strcmp(selectedProjection,'Laplacian Eigenmap')==1
    cmd=[cmd 'LE'];
elseif strcmp(selectedProjection,'Hessian Locally Linear Embedding')==1
    cmd=[cmd 'HLLE'];
elseif strcmp(selectedProjection,'Stochastic Proximity Embedding')==1
    cmd=[cmd 'SPE'];
elseif strcmp(selectedProjection,'Neighborhood Preserving Embedding')==1
    cmd=[cmd 'NPE --saveMapping ' fnProjector];
end
if exist(fnProjector,'file')==2
    system(['rm -f ' fnProjector]);
end
system(cmd);
handles.NMAdisplacementsProjected=load(handles.fnProjected);
updateListBox(gcbo, handles);
set(handles.listRepresentation,'Value',1:min(2,str2num(dout)))
guidata(gcbo,handles)
updatePlot(handles)

% --- Executes during object creation, after setting all properties.
function listRepresentation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listRepresentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in listRepresentation.
function listRepresentation_Callback(hObject, eventdata, handles)
% hObject    handle to listRepresentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listRepresentation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listRepresentation

% --- Executes on button press in pushbuttonRepresentation.
function pushbuttonRepresentation_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRepresentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updatePlot(handles)

% --- Executes during object creation, after setting all properties.
function ndimensions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ndimensions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function updateListBox(hObject, handles)
    listboxString={};
    for i=1:size(handles.NMAdisplacementsProjected,2)
        listboxString{i}=['X' num2str(i)];
    end
    set(handles.ndimensions,'String',num2str(size(handles.NMAdisplacementsProjected,2)));
    set(handles.listRepresentation,'String',listboxString);
    guidata(hObject,handles)

function updatePlot(handles)
    idx=get(handles.listRepresentation,'Value');
    figure(handles.figHandle)
    if length(idx)==1
        hist(handles.NMAdisplacementsProjected(:,idx),50);
        xlabel(['X' num2str(idx)])
    elseif length(idx)==2
        plot(handles.NMAdisplacementsProjected(:,idx(1)),handles.NMAdisplacementsProjected(:,idx(2)),'o');
        xlabel(['X' num2str(idx(1))])
        ylabel(['X' num2str(idx(2))])
        grid on
        axis square
    elseif length(idx)==3
        plot3(handles.NMAdisplacementsProjected(:,idx(1)),handles.NMAdisplacementsProjected(:,idx(2)),...
            handles.NMAdisplacementsProjected(:,idx(3)),'o');
        xlabel(['X' num2str(idx(1))])
        ylabel(['X' num2str(idx(2))])
        zlabel(['X' num2str(idx(3))])
        grid on
        axis square
    end
