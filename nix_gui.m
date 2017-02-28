function varargout = nix_gui(varargin)
% NIX_GUI MATLAB code for nix_gui.fig
%      NIX_GUI, by itself, creates a new NIX_GUI or raises the existing
%      singleton*.
%
%      H = NIX_GUI returns the handle to a new NIX_GUI or the handle to
%      the existing singleton*.
%
%      NIX_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NIX_GUI.M with the given input arguments.
%
%      NIX_GUI('Property','Value',...) creates a new NIX_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before nix_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to nix_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help nix_gui

% Last Modified by GUIDE v2.5 25-Jan-2017 16:44:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @nix_gui_OpeningFcn, ...
    'gui_OutputFcn',  @nix_gui_OutputFcn, ...
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
if isempty(varargin),
    
    global LUE
    
    if ispc, start_stats = [5,300,140,20];
        elseif ismac, start_stats = [5,285,140,20];
        else,    start_stats = [10,320,160,20]; end;
    
    if LUE.type == 2,
        
        LUE.handles.more(1) = uicontrol('Parent',LUE.handles.gui,'Style','Text','Position',start_stats+[start_stats(3),0,0,0],'String',sprintf('Factor-Name'),'HorizontalAlignment','center','BackgroundColor',get(LUE.handles.gui,'color'));
        LUE.handles.more(2) = uicontrol('Parent',LUE.handles.gui,'Style','Text','Position',start_stats+[start_stats(3)*2,0,0,0],'String',sprintf('Factor-Levels'),'HorizontalAlignment','center','BackgroundColor',get(LUE.handles.gui,'color'));
        
        if LUE.amountfactors(1) == 0, LUE.amountfactors(1) = 1; later = 1; else, later = 0; end;
        
        for i = 1:LUE.amountfactors(1),
            LUE.handles.bt(i).text  = uicontrol('Parent',LUE.handles.gui,'Style','Text','Position',start_stats+[0,0-i*20,0,0],'String',sprintf('Between-Factor #% 2d:',i),'HorizontalAlignment','right','BackgroundColor',get(LUE.handles.gui,'color'));
            LUE.handles.bt(i).edit1 = uicontrol('Parent',LUE.handles.gui,'Style','Edit','Position',start_stats+[start_stats(3),0-i*20,0,0],'String',sprintf('BetweenFactorNr%d',i),'HorizontalAlignment','right','BackgroundColor',get(LUE.handles.gui,'color'));
            LUE.handles.bt(i).edit2 = uicontrol('Parent',LUE.handles.gui,'Style','Edit','Position',start_stats+[start_stats(3)*2,0-i*20,0,0],'String','2','HorizontalAlignment','right','BackgroundColor',get(LUE.handles.gui,'Color'));
        end;
        
        if later,
            set(LUE.handles.bt(1).text,'Visible','off');
            set(LUE.handles.bt(1).edit1,'String','AllSubjects','Enable','off','Visible','off');
            set(LUE.handles.bt(1).edit2,'String','1','Enable','off','Visible','off');
        end;
        
        start_stats = start_stats + [0,0-(i+1)*20,0,0];
        for i = 1:LUE.amountfactors(2),
            LUE.handles.wt(i).text  = uicontrol('Parent',LUE.handles.gui,'Style','Text','Position',start_stats+[0,0-i*20,0,0],'String',sprintf('Within-Factor #% 2d:',i),'HorizontalAlignment','right','BackgroundColor',get(LUE.handles.gui,'color'));
            LUE.handles.wt(i).edit1 = uicontrol('Parent',LUE.handles.gui,'Style','Edit','Position',start_stats+[start_stats(3),0-i*20,0,0],'String',sprintf('WithinNr%d',i),'HorizontalAlignment','right','BackgroundColor',get(LUE.handles.gui,'color'));
            LUE.handles.wt(i).edit2 = uicontrol('Parent',LUE.handles.gui,'Style','Edit','Position',start_stats+[start_stats(3)*2,0-i*20,0,0],'String','2','HorizontalAlignment','right','BackgroundColor',get(LUE.handles.gui,'Color'));
        end;
        
        %try, addpath(genpath(fullfile(pwd,'spm_subs'))); end;
        %maindir bestimmen

    elseif LUE.type == 1,
        set(LUE.handles.gui,'visible','on');
        
        LUE.handles.more(1) = uicontrol('Parent',LUE.handles.gui,'Style','Text','Position',start_stats+[start_stats(3),0,0,0],'String',sprintf('Dimension-Name'),'HorizontalAlignment','center','BackgroundColor',get(LUE.handles.gui,'color'));
        LUE.handles.more(2) = uicontrol('Parent',LUE.handles.gui,'Style','Text','Position',start_stats+[start_stats(3)*2,0,0,0],'String',sprintf('Dimension-Levels'),'HorizontalAlignment','center','BackgroundColor',get(LUE.handles.gui,'color'));
        
        for i = 1:LUE.amountfactors,
            LUE.handles.bt(i).text  = uicontrol('Parent',LUE.handles.gui,'Style','Text','Position',start_stats+[0,0-i*20,0,0],'String',sprintf('Dimension #% 2d:',i),'HorizontalAlignment','right','BackgroundColor',get(LUE.handles.gui,'color'));
            LUE.handles.bt(i).edit1 = uicontrol('Parent',LUE.handles.gui,'Style','Edit','Position',start_stats+[start_stats(3),0-i*20,0,0],'String',sprintf('Dim#%d',i),'HorizontalAlignment','right','BackgroundColor',get(LUE.handles.gui,'color'));
            LUE.handles.bt(i).edit2 = uicontrol('Parent',LUE.handles.gui,'Style','Edit','Position',start_stats+[start_stats(3)*2,0-i*20,0,0],'String','2','HorizontalAlignment','right','BackgroundColor',get(LUE.handles.gui,'Color'));
        end;
        
    end;
    a = which('nix');
    b = strfind(a,filesep);
    LUE.maindir = a(1:b(end)-1);
end;

% --- Executes just before nix_gui is made visible.
function nix_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to nix_gui (see VARARGIN)

% Choose default command line output for nix_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
clearvars -global LUE
global LUE
LUE = [];
LUE.groups = []; LUE.aim = []; LUE.within = [];  LUE.wtnames = []; LUE.handles.gui = hObject; LUE.wtnamesxls = []; % LUE.maindir = pwd; 
LUE.type   = nix_menu2('Kind of Statistics','Contingency Table','Repeated non-parametric data');

if LUE.type == 2,
    checksum = 0;
    while checksum==0,
        checksum = 1;
        amountfactors = inputdlg({'Number of Between Factors','Number of Within Factors'},'Insert Number of Factors',1,{'1','1'});
        if (str2num(amountfactors{1})<0) | (str2num(amountfactors{1})>1)
            fprintf('The number of Between-factors have to be 0 or 1\n');
            checksum = 0;
        else,
            LUE.amountfactors(1) = str2num(amountfactors{1});
        end;
        if (str2num(amountfactors{2})<1) | (str2num(amountfactors{2})>2),
            fprintf('The number of Within-factors have to be 1 or 2\n');
            checksum = 0;
        else,
            LUE.amountfactors(2) = str2num(amountfactors{2});
        end;
        %     for i = 1:length(amountfactors),
        %         if str2num(amountfactors{i})<1,
        %             fprintf('You can''t choose anything under 0\n');
        %             checksum = 0; break;
        %         end;
        %         LUE.amountfactors(i) = str2num(amountfactors{i});
        %     end;
    end;
elseif LUE.type == 1,
    
    checksum = 0;
    while checksum==0,
    
        amountfactors = inputdlg({'Number of Dimensions (besides lesion data)'},'Choose the Dimensions',1,{'2'});
        LUE.amountfactors = str2num(amountfactors{1});
        
        if LUE.amountfactors > 1,
            checksum = 1;
        else,
            fprintf('The number of dimensions have to be greater than 1\n');
        end;
        
    end;
    
    LUE.within = 1;
    
%     %Name of Dimensions
%     astring = 'LUE.factorsnames = inputdlg({'; astring2 = '';
%     for i = 1 : LUE.amountfactors,
%         astring  = [astring,sprintf('''%d. Dimension'',',i)];
%         astring2 = [astring2, sprintf('''Dim%d'',',i)];
%     end; 
%     astring = [astring(1:end-1), sprintf('},''Enter the Name for each Dimension'',1,{%s});',astring2(1:end-1))];
%     eval(astring);
%     
%     
%     %Number of Dimensions
%     astring = 'LUE.amountfactors = inputdlg({'; astring2 = '';
%     for i = 1 : LUE.amountfactors,
%         astring  = [astring,sprintf('''%s'',',LUE.factorsnames{i})];
%         astring2 = [astring2, '''2'','];
%     end; 
%     astring = [astring(1:end-1), sprintf('},''Size of Dimensions'',1,{%s});',astring2(1:end-1))];
%     eval(astring);
%     LUE.amountfactors = str2num(cell2mat(LUE.amountfactors))';
   
end;

% UIWAIT makes nix_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = nix_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
nix_zusammenstellen(handles);
try, set(handles.gen,'String',num2str(str2num(get(handles.gen,'string'))+1)); end;


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isequal(get(handles.listbox1,'visible'),'on'),
    afiles = strvcat(' ',spm_select([0,Inf],'image'));
    set(handles.listbox1,'String',afiles);
elseif isequal(get(handles.uitable1,'visible'),'on'),
    [afile,apath] = uigetfile('*','Choose File with Within-Data');
    if (afile),
        [columnnames,data] = nix_txteinlesen(fullfile(apath,afile));
        
        anzahl = size(get(handles.uitable1,'data'),1);
        if ~isequal(size(data,1),anzahl),
            msgbox(sprintf('%d datasets found. %d were expected. Check your Data-File or your chosen Files for the group.',size(data,1),anzahl),'Wrong Number of Data','warn');
        else
            set(handles.uitable1,'ColumnName',columnnames);
            set(handles.uitable1,'Data',data);
        end;
    end;
end;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
timepoint_reshuffle(handles);


% --- Executes on button press in check_multicore.
function check_multicore_Callback(hObject, eventdata, handles)
% hObject    handle to check_multicore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of check_multicore
