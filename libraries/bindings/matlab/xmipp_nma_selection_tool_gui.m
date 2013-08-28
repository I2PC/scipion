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

% Last Modified by GUIDE v2.5 28-Aug-2013 20:45:33

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
handles.fnState=[handles.rundir '/extra/projectionConfig.mat'];

% Open NMA results

[handles.images, handles.NMAdisplacements, handles.cost]= xmipp_nma_read_alignment(handles.rundir);
%handles.NMAdisplacements=rand(100,4);

%handles.cost=rand(size(handles.NMAdisplacements,1),1);
%for i=1:length(handles.cost)
%    handles.images{i}=['file ' int2str(i) '.xmp'];
%end

if exist(handles.fnProjected,'file')
    handles.NMAdisplacementsProjected=load(handles.fnProjected);
else
    handles.NMAdisplacementsProjected=handles.NMAdisplacements;
end

if exist(handles.fnState,'file')
    handles=loadState(handles);
else
    handles.inCluster=zeros(size(handles.NMAdisplacementsProjected,1),1);
    handles.included=ones(size(handles.NMAdisplacementsProjected,1),1);
    updateListBox(hObject, handles);
    set(handles.listRepresentation,'Value',[1 2])
    saveState(handles);
end
handles.trajectoryX=[];
handles.trajectoryY=[];
handles.figHandle=figure();
guidata(hObject,handles)
updatePlot(handles)

function handlesOut=loadState(handles)
    handlesOut=handles;
    load(handles.fnState);
    set(handlesOut.popupProjection,'Value',projectionMethod);
    set(handlesOut.listRepresentation,'Value',representationIdx);
    set(handlesOut.ndimensions,'String',ndimensions);
    set(handlesOut.listRepresentation,'String',listboxString);
    handlesOut.inCluster=inCluster;
    handlesOut.included=included;

function saveState(handles)
    projectionMethod=get(handles.popupProjection,'Value');
    representationIdx=get(handles.listRepresentation,'Value');
    ndimensions=get(handles.ndimensions,'String');
    listboxString=get(handles.listRepresentation,'String');
    inCluster=handles.inCluster;
    included=handles.included;
    save(handles.fnState,'projectionMethod','representationIdx','ndimensions',...
        'listboxString','inCluster','included');

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
saveState(handles)

% --- Executes on button press in exportWorkspace.
function exportWorkspace_Callback(hObject, eventdata, handles)
% hObject    handle to exportWorkspace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    assignin('base','imgs',handles.images)
    assignin('base','NMAdisplacements',handles.NMAdisplacements)
    assignin('base','NMAdisplacementsProjected',handles.NMAdisplacementsProjected)
    assignin('base','cost',handles.cost)
    disp('The following variables have been created in the workspace')
    disp('   images:                    cell array with the image names')
    disp('   NMAdisplacements:          matrix with the raw displacements')
    disp('   NMAdisplacementsProjected: matrix with the displacements projected onto a lower dimensional space')
    disp('   cost:                      cost of each one of the images (lower costs, better assignments)')
    
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

function conditionExcludeString_Callback(hObject, eventdata, handles)
% hObject    handle to conditionExcludeString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of conditionExcludeString as text
%        str2double(get(hObject,'String')) returns contents of conditionExcludeString as a double
    for i=1:size(handles.NMAdisplacementsProjected,2)
        eval(['X' num2str(i) '=handles.NMAdisplacementsProjected(:,' num2str(i) ');']);
    end
    eval(['idx=' get(handles.conditionExcludeString,'String') ';']);
    handles.included(idx)=0;
    guidata(gcbo,handles)
    updatePlot(handles)

% --- Executes during object creation, after setting all properties.
function conditionExcludeString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to conditionExcludeString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in excludeReset.
function excludeReset_Callback(hObject, eventdata, handles)
% hObject    handle to excludeReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.included=ones(size(handles.included));
    guidata(gcbo,handles)
    updatePlot(handles)

% --- Executes on button press in clusterReset.
function clusterReset_Callback(hObject, eventdata, handles)
% hObject    handle to clusterReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    handles.inCluster=zeros(size(handles.inCluster));
    guidata(gcbo,handles)
    updatePlot(handles)

% --- Executes on button press in clusterLoad.
function clusterLoad_Callback(hObject, eventdata, handles)
% hObject    handle to clusterLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fileList=dir([handles.rundir '/extra/cluster__*.mat']);
    if length(fileList)==0
        return
    end
    listString={};
    for i=1:length(fileList)
        aux=strrep(fileList(i).name,'.mat','');
        listString{i}=strrep(aux,'cluster__','');
    end
    [Selection,ok] = listdlg('ListString',listString,'SelectionMode','single');
    if ok
        load([handles.rundir '/extra/cluster__' listString{Selection} '.mat'])
        handles.inCluster=inCluster;
        guidata(gcbo,handles)
        updatePlot(handles)
    end

% --- Executes on button press in clusterSave.
function clusterSave_Callback(hObject, eventdata, handles)
% hObject    handle to clusterSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    clusterName=inputdlg('Name of the cluster','Save cluster',1);
    if ~isempty(clusterName)
        inCluster=handles.inCluster;
        save([handles.rundir '/extra/cluster__' clusterName{1} '.mat'],'inCluster')
        xmipp_nma_save_cluster(handles.rundir,clusterName{1},inCluster);
    end

function conditionString_Callback(hObject, eventdata, handles)
% hObject    handle to conditionString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of conditionString as text
%        str2double(get(hObject,'String')) returns contents of conditionString as a double
    % Produce Xi variables
    for i=1:size(handles.NMAdisplacementsProjected,2)
        eval(['X' num2str(i) '=handles.NMAdisplacementsProjected(:,' num2str(i) ');']);
    end
    eval(['idx=' get(handles.conditionString,'String') ';']);
    handles.inCluster(idx)=1;
    guidata(gcbo,handles)
    updatePlot(handles)

% --- Executes during object creation, after setting all properties.
function conditionString_CreateFcn(hObject, eventdata, handles)
% hObject    handle to conditionString (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in freehand.
function freehand_Callback(hObject, eventdata, handles)
% hObject    handle to freehand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    idxVars=get(handles.listRepresentation,'Value');
    if length(idxVars)~=2
        return
    end

    k = waitforbuttonpress;
    point1 = get(gca,'CurrentPoint');    % button down detected
    finalRect = rbbox;                   % return figure units
    point2 = get(gca,'CurrentPoint');    % button up detected
    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions

    X=handles.NMAdisplacementsProjected(:,idxVars(1));
    Y=handles.NMAdisplacementsProjected(:,idxVars(2));
    idx = X>=p1(1) & X<=(p1(1)+offset(1)) & Y>=p1(2) & Y<=(p1(2)+offset(2));
    handles.inCluster(idx)=1;
    guidata(gcbo,handles)
    updatePlot(handles)

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
    c=handles.inCluster;
    in=handles.included;
    hold off
    if length(idx)==1
        x=handles.NMAdisplacementsProjected(:,idx);
        m=min(x);
        M=max(x);
        N=50;
        d=(M-m)/N;
        bins=m+d/2+[0:(N-1)]*d;
        hist(x(find(c==0 & in)),bins);
        hold on
        hist(x(find(c==1 & in)),bins);
        h = findobj(gca,'Type','patch');
        set(h(1),'FaceColor','b');
        set(h(2),'FaceColor','r');
        xlabel(['X' num2str(idx)])
    elseif length(idx)==2
        idxc=find(c==0 & in);
        X=handles.NMAdisplacementsProjected(:,idx(1));
        Y=handles.NMAdisplacementsProjected(:,idx(2));
        sizePointsX=(max(abs(X))-min(abs(X)))/15;
        sizePointsY=(max(abs(Y))-min(abs(Y)))/15;
        sizePoints=max(5,min([sizePointsX,sizePointsY]));
        scatter(X(idxc),Y(idxc),sizePoints,handles.cost(idxc),'o');
        hold on
        idxc=find(c==1 & in);
        plot(X(idxc),Y(idxc),'*b');
        h=colorbar;
        ylabel(h,'Error')
        xlabel(['X' num2str(idx(1))])
        ylabel(['X' num2str(idx(2))])
        grid on
        axis square
        
        if length(handles.trajectoryX)>0
            plot(handles.trajectoryX,handles.trajectoryY,'LineWidth',2);
        end
    elseif length(idx)==3
        idxc=find(c==0 & in);
        X=handles.NMAdisplacementsProjected(:,idx(1));
        Y=handles.NMAdisplacementsProjected(:,idx(2));
        Z=handles.NMAdisplacementsProjected(:,idx(3));
        sizePointsX=(max(abs(X))-min(abs(X)))/15;
        sizePointsY=(max(abs(Y))-min(abs(Y)))/15;
        sizePointsZ=(max(abs(Z))-min(abs(Z)))/15;
        sizePoints=max(5,min([sizePointsX,sizePointsY,sizePointsZ]));
        scatter3(X(idxc),Y(idxc),Z(idxc),sizePoints,handles.cost(idxc),'o');
        hold on
        idxc=find(c==1 & in);
        plot3(X(idxc),Y(idxc),Z(idxc),'*b');
        xlabel(['X' num2str(idx(1))])
        ylabel(['X' num2str(idx(2))])
        zlabel(['X' num2str(idx(3))])
        h=colorbar;
        ylabel(h,'Error')
        grid on
        axis square
    end
    set(handles.clusterSizeText,'String',[num2str(sum(c)) '/' num2str(length(c)) ' images'])
    saveState(handles)


% --- Executes on button press in drawTrajectory.
function drawTrajectory_Callback(hObject, eventdata, handles)
% hObject    handle to drawTrajectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    idxVars=get(handles.listRepresentation,'Value');
    if length(idxVars)~=2
        return
    end

    k = waitforbuttonpress;
    [handles.trajectoryX,handles.trajectoryY]=getline();
    guidata(gcbo,handles);
    updatePlot(handles);

% --- Executes on button press in generateAnimation.
function generateAnimation_Callback(hObject, eventdata, handles)
% hObject    handle to generateAnimation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in loadTrajectory.
function loadTrajectory_Callback(hObject, eventdata, handles)
% hObject    handle to loadTrajectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    fileList=dir([handles.rundir '/extra/trajectory__*.mat']);
    if length(fileList)==0
        return
    end
    listString={};
    for i=1:length(fileList)
        aux=strrep(fileList(i).name,'.mat','');
        listString{i}=strrep(aux,'trajectory__','');
    end
    [Selection,ok] = listdlg('ListString',listString,'SelectionMode','single');
    if ok
        load([handles.rundir '/extra/trajectory__' listString{Selection} '.mat'])
        handles.trajectoryX=trajectoryX;
        handles.trajectoryY=trajectoryY;
        guidata(gcbo,handles)
        updatePlot(handles)
    end

% --- Executes on button press in saveTrajectory.
function saveTrajectory_Callback(hObject, eventdata, handles)
% hObject    handle to saveTrajectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    trajectoryName=inputdlg('Name of the trajectory','Save trajectory',1);
    if ~isempty(trajectoryName)
        trajectoryX=handles.trajectoryX;
        trajectoryY=handles.trajectoryY;
        save([handles.rundir '/extra/trajectory__' trajectoryName{1} '.mat'],'trajectoryX','trajectoryY')
    end

% --- Executes on button press in resetTrajectory.
function resetTrajectory_Callback(hObject, eventdata, handles)
% hObject    handle to resetTrajectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        handles.trajectoryX=[];
        handles.trajectoryY=[];
        guidata(gcbo,handles)
        updatePlot(handles)
