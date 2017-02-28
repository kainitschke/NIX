function nix_show_results(varargin)%(effectpath),

global hans
try, delete(hans.fh); end;
clearvars -global hans LUE
global LUE hans inn

%% check for SPM
if isempty(which('spm')),
    uiwait(warndlg(sprintf('SPM was not found. Some subroutines of SPM are needed in the NIX-toolbox. SPM is freeware by the Wellcome Trust Centre for Neuroimaging (www.fil.ion.ucl.ac.uk/spm/).\nYou will now be asked to locate SPM on your system. Without SPM the NIX-toolbox will not work.\nTo avoid this message add SPM permanently to your start paths of Matlab.'),'SPM missing','modal'));
    spm_path = uigetdir(pwd,'Select SPM directory');
    addpath(genpath(spm_path));
    try, spm('fMRI'); catch, errordlg('There is something wrong with your SPM. Without it, the NIX-tool is not gonna work. Please correct these mistakes. Unfortunately, I cannot do anything, because the problem occurs in SPM. I am afraid, you have to go to their community for this problem.','SPM ERROR'); end;
else,
    spm_path = which('spm.m');
    weg      = strfind(spm_path,filesep);
    spm_path = spm_path(1:weg(end)-1);
end;

%% 
pverh = 2/3; pvertausgleich = 90;

if length(varargin)<1,
    hans.effectpath = uigetdir(pwd,'Select Original Result Directory');
else,
    hans.effectpath = varargin{1};
end;

%% choose underlay
hans.underlay = [];
if (length(varargin)<2) & (exist(fullfile(hans.effectpath,'underlay.nii'))), % detect previous underlay
    weg = questdlg('Previous underlay picture found.', 'Underlay', 'Previous', 'Other (choose myself)', 'Previous');
    if isequal(weg,'Previous'), 
        preunderlay = 1; 
    else, 
        preunderlay = 0; 
    end;
else,
    preunderlay = 0;
end;

if (length(varargin)<2) & (preunderlay == 0),
    [ufile,upath] = uigetfile({'*.nii','.hdr'},'Select Underlying Picture');
    hans.underlayfile = fullfile(upath,ufile);
elseif preunderlay == 1,
    hans.underlayfile = fullfile(hans.effectpath,'underlay.nii');
else,
    hans.underlayfile = varargin{2};
end;

%% after chosing, so you can lean back
inn = matfile(fullfile(hans.effectpath,'for_result_tool.mat'));
LUE = inn.LUE;

LUE.resultdir = hans.effectpath;
LUE.maindir   = which('nix_show_results.m'); weg = strfind(LUE.maindir,filesep); LUE.maindir = LUE.maindir(1:weg(end)-1);

sample_file = spm_select('List',LUE.resultdir,'^p.*.nii');; %p_ImgData.nii
try, if iscell(sample_file), sample_file = sample_file{1}; end; end;
sample_file = deblank(sample_file(1,:));


%% get Atlas and reshape | preatlas = [PreviousInFolder  Choose Standard]
if (~exist(fullfile(LUE.resultdir,'atlas.nii'))) | (~exist(fullfile(LUE.resultdir,'atlas.txt'))), % is there an atlas file in the result folder
    preatlas(1)    = 0;
    atlasfile      = 'atlas71.nii';
    hans.atlasfile = spm_select('FPListRec',spm_path,atlasfile);
    if iscell(hans.atlasfile), hans.atlasfile = hans.atlasfile{1}; end;
    try, atlaspath      = strfind(hans.atlasfile,filesep); atlaspath = hans.atlasfile(1:atlaspath(end)-1); end;
else,
    preatlas(1)    = 1;
    atlasfile      = 'atlas.nii';
    hans.atlasfile = fullfile(LUE.resultdir,'atlas.nii');
end;

if ~isempty(hans.atlasfile),
    weg = questdlg('Which atlas file should be used?', 'Atlas', 'Standard/Previous', 'Other (choose myself)', 'Standard/Previous');
    if isequal(weg,'Standard/Previous'), 
        preatlas(2) = 1; 
    else, 
        preatlas(2) = 0; 
    end;
else,
    preatlas(2) = 0;
end;

if preatlas(2) == 0,
    [atlasfile,atlaspath] = uigetfile({'*.nii;*.hdr','Imaging-File (*.nii, *.hdr)'},'Did not find standard atlas (atlas71.nii). Please another.',spm_path);
    hans.atlasfile        = fullfile(atlaspath, atlasfile);
end;

% if isempty(hans.atlasfile) & (preatlas==0),
%     [atlasfile,atlaspath] = uigetfile({'*.nii;*.hdr','Imaging-File (*.nii, *.hdr)'},'Did not find standard atlas (atlas71.nii). Please another.',spm_path);
%     hans.atlasfile = fullfile(atlaspath, atlasfile);
% else,
%     hans.atlasfile = deblank(hans.atlasfile(1,1:end));
%     atlaspath = strfind(hans.atlasfile,filesep); atlaspath = hans.atlasfile(1:atlaspath(end)-1); 
% end;

if ~exist(fullfile(hans.atlasfile)),
    atlaspath      = LUE.resultdir; atlasfile = 'atlas.nii';
    fid            = fopen(fullfile(atlaspath,[atlasfile(1:end-4),'.txt']),'w+'); fwrite(fid,sprintf('1\t Emtpy\n')); fclose(fid);
    Vl.fname       = fullfile(atlaspath,atlasfile);
    Vl.dim         = [50,50,50];
    Vl.mat         = eye(4); Vl.mat(1:3,end) = [-25,-25,-25];
    Vl             = spm_create_vol(Vl);
    spm_write_vol(Vl,zeros(50,50,50));
    hans.atlasfile = fullfile(atlaspath,atlasfile);
    %preatlas       = [1 1];
end;

if (preatlas(1) == 0) | (preatlas(2) == 0),
    weg = strfind(atlasfile,'.'); atlasfile = atlasfile(1:weg(end)-1);
    if ~exist(fullfile(atlaspath,[atlasfile,'.txt'])), errordlg('The text file corresponding to the atlasfile was not found. Please search/download manually and restart the nix_show_results'); end;
    if ~isequal(fullfile(atlaspath,[atlasfile,'.txt']),fullfile(LUE.resultdir,'atlas.txt'))
        copyfile(fullfile(atlaspath,[atlasfile,'.txt']),fullfile(LUE.resultdir,'atlas.txt'));
    end;
    
    spm_jobman('serial','batch_imcalc_nearest.m','', ...
        {fullfile(LUE.resultdir,sample_file); hans.atlasfile}, ... File i1 for space; i2 is transformed
        'atlas.nii', ... save name
        {LUE.resultdir}); ...  result dir
    
    hans.atlasfile = fullfile(LUE.resultdir,'atlas.nii');
end;
hans.atlas     = [];

%% Underlay reshape
if preunderlay == 0,
    spm_jobman('serial','batch_imcalc_trilin.m','', ...
        {fullfile(LUE.resultdir,sample_file); hans.underlayfile}, ... File i1 for space; i2 is transformed
        'underlay.nii', ... save name
        {LUE.resultdir}); ...  result dir
        hans.underlayfile = fullfile(LUE.resultdir,'underlay.nii');
end;

%% 
if ~isfield(LUE,'type'), LUE.type = 2; end; %notwendig fuer alte Version

if LUE.type == 1,
    LUE.namesfinal{end+1} = 'ImgData';
    LUE.fac               = [LUE.fac, LUE.lesiontypes];
    LUE.levelnames{end+1} = num2cell(LUE.lesiontypec)';
    for i = 1:length(LUE.levelnames{end}),
        LUE.levelnames{end}{i} = num2str(LUE.levelnames{end}{i});
    end;
end;

addpath(fullfile(LUE.maindir,'ranktest'));
addpath(fullfile(LUE.maindir,'src'));

%% Figure Größe anpassen
hans.fh = figure;
set(hans.fh,'Resize','Off','MenuBar','none','Color',[1,1,1]);
res = nix_res(pvertausgleich,pverh);
if (diff(res([1,3]))+1)*pverh >= diff(res([2,4]))+1,%res(1)*pverh >= res(2),
%     psize(1) = round(res(2)/pverh);
%     psize(2) = res(2);
%     set(hans.fh,'Position',[round((res(1)-psize(1))/2),round(pvertausgleich/2),psize(1),psize(2)]);
    psize(1) = round((diff(res([2,4]))+1) / pverh ) ;%- pvertausgleich;
    psize(2) = diff(res([2,4]))+1 - pvertausgleich;
    set(hans.fh,'Position',[res(1) + round((diff(res([1,3]))+1-psize(1))/2) + 1, round(pvertausgleich/2) + 1, psize(1), psize(2)]);
else,
%     psize(1) = res(1);
%     psize(2) = round(res(1)*pverh);
%     set(hans.fh,'Position',[1,round((res(2)-psize(2))/2),psize(1),psize(2)]);
    psize(1) = diff(res([1,3])) + 1 ;%- pvertausgleich ;
    psize(2) = round((diff(res([1,3]))+1) * pverh);%- pvertausgleich;
    set(hans.fh,'Position',[res(1)+pvertausgleich/2,  round((diff(res([2,4]))+1-psize(2))/2),   psize(1),  psize(2)]);
end;

%% Alles Grundlegende anlegen0
%Axes für Bilder
weg = axes('position',[0,.5,1,.5],'Color',[0,0,0],'XAxisLocation','top'); 
for i = 1:3,
    hans.und(i) = axes('position',[1/3*(i-1),.5,1/3,.5]);
    hans.brah(i) = axes('position',[1/3*(i-1),.5,1/3,.5],'Color','none');
    hans.cross(i) = axes('position',[1/3*(i-1),.5,1/3,.5],'Color','none');
end;
hans.raw = axes('Position',[1/4+.05,0+.05,3/4-.07,.5-.1]);

uisizes = [10,round(psize(2)/2)-140,60,20];

%edits für x,y,z
vstr = {'z','y','x'};
uicontrol('Parent',hans.fh,'style','text','String','Coordinates','position',(uisizes)+[0,uisizes(4)*(4.1),20,0],'HorizontalAlignment','left','BackgroundColor',[1,1,1]);
for i = 3:-1:1,
    uicontrol('Parent',hans.fh,'style','text','String',vstr{i},'position',(uisizes-[6,0,0,0])+[0,uisizes(4)*(i),0,0],'HorizontalAlignment','left','BackgroundColor',[1,1,1]);
    hans.edit(4-i) = uicontrol('parent',hans.fh,'style','edit','string','0','Position',uisizes+[6,uisizes(4)*(i),0,0],'BackgroundColor',[1,1,1]);
    set(hans.edit(4-i),'KeyPressFcn',@(h_obj,evt) nix_pressenter(evt.Key));
end;

% if LUE.type == 2,
    hans.intbutton = uicontrol('Parent',hans.fh,'style','pushbutton','String','Mult. Coord.','Position',uisizes+[6,uisizes(4)*(0),20,0]);
    set(hans.intbutton,'Callback','nix_pinteraction(2);');
% end;

%Checkbox für Post-Hoc
%hans.postcheck = uicontrol('Parent',hans.fh,'Style','Checkbox','String','Permutation Tests','BackgroundColor',[1,1,1],'position',[uisizes]+[0,0,135,0],'Callback','global hans; if get(hans.postcheck,''Value''), set(hans.welchcheck,''Value'',1); end;');
%hans.welchcheck = uicontrol('Parent',hans.fh,'Style','Checkbox','String','Welch Tests','BackgroundColor',[1,1,1],'position',[uisizes]+[0,-20,135,0],'Callback','global hans; if get(hans.welchcheck,''Value'')==0, set(hans.postcheck,''Value'',0); end;');

%Einstellungen für Overlay
vstr = {'max Transp','Max','Min'}; vor_wert = [.8,10,1.5];
rvers = 185;
uicontrol('Parent',hans.fh,'style','text','String','Threshold','position',(uisizes)+[rvers-50+20,uisizes(4)*(4.1),20,0],'HorizontalAlignment','left','BackgroundColor',[1,1,1]);
for i = 3:-1:1,
    uicontrol('Parent',hans.fh,'style','text','String',vstr{i},'position',(uisizes-[12+105,0,-0,0])+[0+rvers,uisizes(4)*(i),20,0],'HorizontalAlignment','right','BackgroundColor',[1,1,1]);
    hans.thres(4-i) = uicontrol('parent',hans.fh,'style','edit','string',num2str(vor_wert(i)),'Position',uisizes+[0+rvers-50+20,uisizes(4)*(i),0,0],'BackgroundColor',[1,1,1]);
    set(hans.thres(4-i),'KeyPressFcn',@(h_obj,evt) nix_pressenter(evt.Key));
end;

%Refresh Button
hans.refreshbutton = uicontrol('parent',hans.fh,'style','pushbutton','string','Refresh','Position',uisizes+[rvers-50+20,0,0,0]);
set(hans.refreshbutton,'Callback','nix_Refresh;');

%Einstellung für Cluster
vstr = {'all n >=','','FDR p-thres','Uncorr p-thres'}; % {'all n >=','Connections','FDR p-thres','Uncorr p-thres'};
vval = {'2','26','.05','.001'};
uicontrol('Parent',hans.fh,'style','text','String','Contours and Cluster','position',(uisizes)+[rvers*2-50-72+20,uisizes(4)*(4.1),90,0],'HorizontalAlignment','left','BackgroundColor',[1,1,1]);
for i = 4:-1:1%3:-1:1,
    hans.clustui(5-i) = uicontrol('parent',hans.fh,'style','edit','string',vval{i},'Position',uisizes+[0+rvers*2-50+20,uisizes(4)*(i-1),0,0],'BackgroundColor',[1,1,1]);
%     if (i==1)&(LUE.type==1),
%         set(hans.clustui(5-i),'Visible','off');
%     else,
        uicontrol('Parent',hans.fh,'style','text','String',vstr{i},'position',(uisizes-[12+160,0,-50,0])+[0+rvers*2+20,uisizes(4)*(i-1),0,0],'HorizontalAlignment','right','BackgroundColor',[1,1,1]);
%     end;
    set(hans.clustui(5-i),'KeyPressFcn',@(h_obj,evt) nix_pressenter(evt.Key));
end;
set(hans.clustui(3),'Visible','off');

%Listmenu für Effekte
hans.effectpopup = uicontrol('Parent',hans.fh,'style','listbox','position',[uisizes(1),20+(uisizes(2)+[uisizes(4)*(1)]-60)/2,round(psize(1)/4-uisizes(1)*2),(uisizes(2)+[uisizes(4)*(1)]-60)/2]);%uisizes(2)+[uisizes(4)*(1)]-60]);

hans.efiles = dir(fullfile(hans.effectpath,'F_*')); hans.efiles = {hans.efiles.name}; for i = 1:length(hans.efiles), hans.efiles{i} = hans.efiles{i}(3:end-4); end; %_F_Files
holder = dir(fullfile(hans.effectpath,'NoVar*')); hans.efiles = [{holder(1).name(1:end-4)},hans.efiles];
hans.efilestype = ones(length(hans.efiles),1);
hans.maxinteraction = [0,0]; for i = 1:length(hans.efiles), if length(strfind(hans.efiles{i},'X'))>hans.maxinteraction(2), hans.maxinteraction = [i,length(strfind(hans.efiles{i},'X'))]; end; end; 
if hans.maxinteraction(2) > 0, hans.maxinteraction = hans.maxinteraction(1); else, hans.maxinteraction = 1; end;
holder = dir(fullfile(hans.effectpath,'Diff_*')); holder = {holder.name}; for i = 1:length(holder), holder{i} = holder{i}(6:end-4); end; hans.efiles = [hans.efiles,holder]; hans.efilestype = [hans.efilestype; ones(length(holder),1)+1]; %_Post-Hoc_files
set(hans.effectpopup,'String',hans.efiles,'Callback','nix_Refresh;');
hans.alreadyloaded = zeros(length(hans.efiles),1);
% schwarze balken für später heraussuchen
if LUE.type == 2,
for i = 1:length(hans.efiles), if hans.efilestype(i)==2,
        if LUE.btall > 1, 
            j = 1; 
            apoint = strfind(hans.efiles{i},LUE.btnames{j});
            runner=0; ok = 1; 
                while ok == 1,
                    try, 
                        hol = str2num(hans.efiles{i}(apoint+length(LUE.btnames{j}):apoint+length(LUE.btnames{j})+runner)); 
                        if ~isempty(hol), hans.phausp{i,2} = hol; else, ok = 0; end;
                        runner = runner + 1;
                    catch, 
                        ok = 0; 
                    end;
            end;
        else, hans.phausp{i,2} = 1;
        end;
        for j = 1:length(LUE.wtnames),
            apoint = strfind(hans.efiles{i},LUE.wtnames{j});
            for l = 1:2,runner=0; ok = 1; while ok == 1,
                    try, hans.phausp{i,1}(l,j) = str2num(hans.efiles{i}(apoint(l)+length(LUE.wtnames{j}):apoint(l)+length(LUE.wtnames{j})+runner)); runner = runner + 1;
                    catch, ok = 0; end;
end;    end;end;    end;                  end;
end;

for i = 1:length(hans.efiles), hans.clust(i).conn = 0;  hans.clust(i).clust = []; hans.clust(i).out = []; end;

%Listmenu für Peaks
hans.clusterpopup = uicontrol('Parent',hans.fh,'style','listbox','position',[uisizes(1),20,round(psize(1)/4-uisizes(1)*2),(uisizes(2)+[uisizes(4)*(1)]-60)/2-5]);%uisizes(2)+[uisizes(4)*(1)]-60]);
set(hans.clusterpopup,'Callback','nix_Refresh(4);');

%Dinge für Statistik vorbereiten
%hans.novar = spm_vol(fullfile(hans.effectpath,'NoVar.nii'));  
%hans.nogroupvar = spm_vol(fullfile(hans.effectpath,'NoGroupVar.nii')); 

%Buttons fuer rausschreiben
hans.SaveCluster = uicontrol('parent',hans.fh,'style','pushbutton','string','Save Cluster Info','Position',[uisizes(1) 0 round(psize(1)/4-uisizes(1)*2)/3  uisizes(4)]);
set(hans.SaveCluster,'Callback','nix_write_cluster;');
if LUE.type == 2,
    hans.SavePosthocData = uicontrol('parent',hans.fh,'style','pushbutton','string','Save Post Hoc','Position',[uisizes(1)+round(psize(1)/4-uisizes(1)*2)/3   0   round(psize(1)/4-uisizes(1)*2)/3    uisizes(4)]);
    set(hans.SavePosthocData,'Callback','nix_write_posthoc;');
end;
hans.SaveRawData = uicontrol('parent',hans.fh,'style','pushbutton','string','Save Raw Data','Position',[uisizes(1)+(round(psize(1)/4-uisizes(1)*2)/3)*2   0   round(psize(1)/4-uisizes(1)*2)/3    uisizes(4)]);
set(hans.SaveRawData,'Callback','nix_write_rawdata;');

hans.SaveMultCond = [];

%Colourbar mit eigener Unterlage anlegen
eint = 35;
hans.colourbar = axes('position',[(eint-1)/eint, .5+1/75 ,1/eint, .5-2/75],'Color',[0.5,0.5,0.5],'XAxisLocation','top'); 
hans.last = -1;

%t  = spm_vol(fullfile(hans.effectpath,['F_',hans.efiles{1},'.nii']));
t  = spm_vol(fullfile(hans.effectpath,[hans.efiles{1},'.nii']));
[~,aXYZ]   = spm_read_vols(t);
aXYZ       = round(aXYZ*100)/100;
hans.aXYZ  = aXYZ;

acur = diff(aXYZ');
for i = 1:3
    af = find(acur(:,i)~=0);
    hans.aXYZdir(i) = mode(acur(af,i));
end;

hans.rfilter = spm_vol(fullfile(hans.effectpath,'MinGroupSize.nii'));
hans.rfilter = round(spm_read_vols(hans.rfilter));
hans.nfilter = 0;
hans.puval   = str2num(get(hans.clustui(1),'String'));

% % check ob übereinstimmung underlay, ggf. ändern


%% load data
inn         = load(fullfile(hans.effectpath,LUE.niicontfile));
LUE.imgperm = unique(inn.working);

%% additional voxel analysis
% if LUE.type == 2,
    nix_pinteraction(0);
% end;

%% for closing
set(hans.fh,'CloseRequestFcn','nix_result_close;');

%% First Refresh and picturing
nix_Refresh;