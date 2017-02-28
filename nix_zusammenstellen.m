function nix_zusammenstellen(hand),
% collect of relevant data by GUI
global LUE

%% check whether "ImgData" as name or any double name
ok = 1;
for i = 1:length(LUE.handles.bt), if isequal(get(LUE.handles.bt(i).edit1,'String'),'ImgData'), ok = 0; end; end;
if LUE.type == 2,  for i = 1:length(LUE.handles.wt), if isequal(get(LUE.handles.wt(i).edit1,'String'),'ImgData'), ok = 0; end; end; end;
for i = 1:length(LUE.handles.bt),
    for j = i+1:length(LUE.handles.bt), if isequal(get(LUE.handles.bt(i).edit1,'String'),get(LUE.handles.bt(j).edit1,'String')), ok = -1; end; end;
    if LUE.type == 2, for j = 1:length(LUE.handles.wt), if isequal(get(LUE.handles.bt(i).edit1,'String'),get(LUE.handles.wt(j).edit1,'String')), ok = -1; end; end; end;
end;
if LUE.type == 2, for i = 1:length(LUE.handles.wt), 
    for j = i+1:length(LUE.handles.wt),if isequal(get(LUE.handles.wt(i).edit1,'String'),get(LUE.handles.wt(j).edit1,'String')), ok = -1; end; end;
end; end;

if ~(ok == 1),
    set(hand.gen,'string',num2str(str2num(get(hand.gen,'string')) - 1));
    if ok == -1, warndlg('You must not use the same factor name more than once.','Forbidden Names');
    else,        warndlg('You must not name any factor "ImgData".','Forbidden Names');
    end;
    return;
end;

%% weiter

    
set(hand.listbox1,'Visible','off');
set(hand.pushbutton2,'Visible','off');
set(hand.uitable1,'Visible','off');

runde = str2num(get(hand.gen,'string'));

if LUE.type == 2,
    
switch runde,
    case 0,
        set(hand.uitable1,'ColumnEditable',true);
        set(LUE.handles.more(1),'Visible','off');
        set(LUE.handles.more(2),'Visible','off');
        for i = 1:length(LUE.handles.bt),
            LUE.bt(i) = str2num(get(LUE.handles.bt(i).edit2,'String'));
            LUE.btnames{i} = get(LUE.handles.bt(i).edit1,'String');
            set(LUE.handles.bt(i).text,'Visible','Off');
            set(LUE.handles.bt(i).edit1,'Visible','Off');
            set(LUE.handles.bt(i).edit2,'Visible','Off');
        end;
        
        for i = 1:length(LUE.handles.wt),
            LUE.wt(i) = str2num(get(LUE.handles.wt(i).edit2,'String'));
            LUE.wtnames{i} = get(LUE.handles.wt(i).edit1,'String');
            set(LUE.handles.wt(i).text,'Visible','Off');
            set(LUE.handles.wt(i).edit1,'Visible','Off');
            set(LUE.handles.wt(i).edit2,'Visible','Off');
        end;

        LUE.btall = prod(LUE.bt);
        LUE.wtall = prod(LUE.wt);
        LUE.wtgroups = [];
            
        klrunde = [0:LUE.wtall-1]';
        LUE.wtgroups = [];
        for i = 1:length(LUE.wt),
            LUE.wtgroups = [LUE.wtgroups, mod(floor((klrunde)/prod(LUE.wt(1:i-1))),LUE.wt(i))+1];
        end;
        klrunde = [0:LUE.btall-1]';
        LUE.btgroups = [];
        for i = 1:length(LUE.bt),
            LUE.btgroups = [LUE.btgroups, mod(floor((klrunde)/prod(LUE.bt(1:i-1))),LUE.bt(i))+1];
        end;
        
        if isequal(LUE.bt,1),
            LUE.levelnames{1}{1} = 'AllSubjects';
        else,
            for i = 1 : size(LUE.bt,2),
                for j = 1 : LUE.bt(i),
                    ggc = inputdlg(sprintf('Name for factor level\n%s: %d',LUE.btnames{i},j),'Naming of the factor levels',1,{''});
                    LUE.levelnames{i}{j} = ggc{1};
                end;
            end;
        end;

    case LUE.btall*2+1%LUE.btall+LUE.wtall+1,
        LUE.final = get(hand.uitable1,'Data');
        LUE.namesfinal = get(hand.uitable1,'ColumnName');
        
    otherwise,
        if runde <= LUE.btall,
            %group = mod(runde-1,LUE.btall)+1;
            group = [];
            for i = 1:length(LUE.bt);
                group = [group, mod(floor((runde-1)/prod(LUE.bt(1:i-1))),LUE.bt(i))+1];
            end;
            aim = get(hand.listbox1,'String'); aim = aim(2:end,:);
            if ~isempty(aim),
                LUE.aim = strvcat(LUE.aim, aim);
                LUE.groups = [LUE.groups; repmat(group,size(aim,1),1)];
            else,
                set(hand.gen,'String',num2str(str2num(get(hand.gen,'string'))-1));
                runde = runde - 1;
                msgbox('You have to specify image-files');
            end;
            set(hand.listbox1,'String','');
        else
            bt = runde-LUE.bt;
            data = get(hand.uitable1,'Data');
            names = get(hand.uitable1,'ColumnName');
            LUE.within = [LUE.within; data];
            if isempty(LUE.wtnamesxls), LUE.wtnamesxls = [LUE.wtnamesxls; names']; end;
        end;
        

end;    

%% set next visible objects for input
if (runde)== LUE.btall*(1+sign(LUE.wtall)) % alles eingelesen, bestätigen Seite
    if LUE.type == 2, set(hand.check_multicore,'visible','on'); end;
    set(hand.pushbutton1,'String','Perform Analysis');
    set(hand.infotext,'string','Ready to Perform Analysis. Attend that this will take some time. Don''t close emerging windows!');
    %set(hand.pushbutton3,'Visible','on');
    %aholder = questdlg('Should the values be z-standardized (recommended)?','Z-Standardization','Yes','No','Yes');
    %if isequal(aholder,'Yes'), LUE.Z = 1; else, LUE.Z=0; end;
    LUE.Z = 0;
    data_zurechtsortieren;
    set(hand.uitable1,'Visible','on','Data',LUE.final,'ColumnName',LUE.namesfinal);

elseif (runde) > (LUE.btall*(1+sign(LUE.wtall))), % Performing startet
    acors = get(hand.check_multicore,'Value');
    set(hand.gen,'String',num2str(str2num(get(hand.gen,'string'))-1));
    delete(hand.figure1);

    LUE.av  = 'Value'; pause(.001);
    czc     = feature('numCores'); if czc > 12, czc = 12; end;
    if LUE.type == 2,
        %opft    = inputdlg(sprintf(sprintf('On how many cpus should be calculated on?\nMax of "%d" cores detected on your system.\nPlease be cautious with your selection. If you\nare not sure chose "1".\n',czc),czc),'Multi-Core-Calculation',1,{'1'});
        if acors,
            opft = czc;
        else,
            opft = 1;
        end;
    else,
        opft     = 1;
    end;
    %fcorr  = inputdlg(sprintf(sprintf('Choose a threshold for a FDR-Correction.\nAll results will be saved uncorrected\nand with the inserted correction threshold.\n',feature('numCores')),feature('numCores')),'FDR-Correction',1,{'.05'});
    
    %fdrtype = questdlg('What FDR method should be used?','FDR method','conventional (conservative)','step-up (liberal)','conventional (conservative)');
    %if isequal(fdrtype,'conventional (conservative)'), LUE.fdrtype = 1; else, LUE.fdrtype = 2; end;
    LUE.fdrtype = 1;
    
    %fdrcheck= questdlg(sprintf('What FDR correction should be performed?\nThe standard FDR correction will be too conservative but will be\ncalculated faster. The adapted FDR correction will result in more\naccurate results but needs longer to be calculated.\nFor details please see the manual.'),'FDR method','Standard','Adapted','Adapted');
    LUE.fdrcheck = 2;%if isequal(fdrtype,'Standard'), LUE.fdrcheck = 1; else, LUE.fdrcheck = 2; end;
    
    fcorr   = inputdlg(sprintf(sprintf('Choose a nominal threshold for a FDR-Correction.\nAll results will be saved uncorrected\nand with the inserted correction threshold.\nFor more information please see the manual.\n',feature('numCores')),feature('numCores')),'FDR-Correction',1,{'5'});
    if str2num(fcorr{1})<2, warndlg('The nominal threshold has to be at least 2. It is set to now.','Attention'); fcorr = {'2'}; end;
    
    %aholder = questdlg('What statistic test should be performed?','Statistic Test','Anova-type','Wald-type','Anova-type');
    aholder = 'anova';
    if isequal(aholder,'Wald-type'), LUE.statkind = 'wald'; else, LUE.statkind = 'anova'; end;
    
    %aholder = questdlg('Should extensive Posthoc Tests be performed?','Extensive Posthoc Tests','Yes','No','Yes');
    aholder = 'Yes';
    if isequal(aholder,'Yes'), LUE.praepost = 1; else, LUE.praepost = 0; end;
    
    LUE.now       = datestr(now,'yyyymmdd_HHMMSS');
    LUE.resultdir = uigetdir(pwd,'Choose a directory for results');
    LUE.resultdir = fullfile(LUE.resultdir,LUE.now); 
    LUE.core      = opft;%str2num(opft{1});
    LUE.fcorr     = str2num(fcorr{1});
    nix_perform;

elseif runde >= LUE.btall,  % within tabelle zeigen
    set(hand.uitable1,'Visible','on');
    klrunde = runde-LUE.btall;
    astr    = [];
    for i = 1:length(LUE.bt),
        %astr = [astr, sprintf('%s:%d, ',LUE.btnames{i}, mod(floor((klrunde)/prod(LUE.bt(1:i-1))),LUE.bt(i))+1)];
        astr     = [astr, sprintf('%s:%s, ',LUE.btnames{i}, LUE.levelnames{i}{mod(floor((klrunde)/prod(LUE.bt(1:i-1))),LUE.bt(i))+1})];
        group(i) = mod(floor((klrunde)/prod(LUE.bt(1:i-1))),LUE.bt(i))+1;
    end;astr     = astr(1:end-2);
    set(hand.infotext,'String',sprintf('Insert Within-Data of WT-Varibles for %s',astr));
    
    arows = LUE.groups==repmat(group,size(LUE.groups,1),1); arows = find(sum(arows',1)'==size(group,2));
    astr  = [];
    for i = 1:size(LUE.wtgroups,1),
        holder      = sprintf('%d,',LUE.wtgroups(i,:));
        astr{end+1} = sprintf('(%s)',holder(1:end-1));
    end;
    set(hand.uitable1,'Data',repmat({''},length(arows),LUE.wtall),'ColumnName',astr);
    set(hand.pushbutton2,'Visible','on');

elseif runde < LUE.btall, % between-listbox zeigen
    set(hand.listbox1,'Visible','on');
    if ispc, set(hand.infotext,'Position',[1.5,22,81.833,5.2143]); set(hand.listbox1,'Position',[1.6667,4.5,81.667, 17.5]); end;
    set(hand.pushbutton2,'Visible','on');
    %set(hand.infotext,'String',sprintf('Choose Files that contain dependent measures for group %d',runde+1));
    astr     = [];
    for i = 1:length(LUE.bt),
        %astr = [astr, sprintf('%s:%d, ',LUE.btnames{i}, mod(floor((runde)/prod(LUE.bt(1:i-1))),LUE.bt(i))+1)];
        astr = [astr, sprintf('%s:%s, ',LUE.btnames{i}, LUE.levelnames{i}{mod(floor((runde)/prod(LUE.bt(1:i-1))),LUE.bt(i))+1})];
    end;astr = astr(1:end-2);
    set(hand.infotext,'String',sprintf('Choose Files that contain dependent measures for %s',astr));
end;

elseif LUE.type == 1,
      
    if runde > 0; runde = runde - 1; end;
    
    if runde == 0,
        
        set(LUE.handles.more(1),'Visible','off');
        set(LUE.handles.more(2),'Visible','off');
        
        for i = 1:size(LUE.handles.bt,2),
            set(LUE.handles.bt(i).text, 'Visible','off');
        end;
        
        for i = 1:size(LUE.handles.bt,2),
            LUE.namesfinal{i} = get(LUE.handles.bt(i).edit1,'String');
            set(LUE.handles.bt(i).edit1,'Visible','off');
        end;
        
        for i = 1:size(LUE.handles.bt,2),
            LUE.amountfactors(i) = str2num(get(LUE.handles.bt(i).edit2,'String'));
            set(LUE.handles.bt(i).edit2,'Visible','off');
        end;
        
        LUE.files = cell(LUE.amountfactors);
        
        set(hand.figure1,'visible','off');
        
        naming    = nix_menu2('Name the factor levels?','Rename the factor levels','Automatic level naming');
        
        
            for i = 1 : size(LUE.amountfactors,2),
                for j = 1 : LUE.amountfactors(i),
                    if naming == 1,
                        ggc = inputdlg(sprintf('Name for factor level\n%s: %d',LUE.namesfinal{i},j),'Naming of the factor levels',1,{''});
                        LUE.levelnames{i}{j} = ggc{1};
                    else,
                        LUE.levelnames{i}{j} = num2str(j);
                    end;
                end;
            end;
        
        set(hand.figure1,'visible','on');
        
        LUE.final = cell(0,length(LUE.amountfactors)+1);
        
        %LUE.amountfactors = LUE.amountfactors(end : -1 : 1);
        %LUE.namesfinal    = LUE.namesfinal(end : -1 : 1);
        
    end;
    
    runde = runde + 1;
    
%     for i = 1 : numel(LUE.amountfactors),
%         fac(i) = ceil(runde/prod(LUE.amountfactors(i+1:end)));
%     end;     astring = [];
%     fac = fliplr(fac);

%     astring = [];
% 
%     fac = fullfact(LUE.amountfactors);
%     fac = fac(runde,:);
%     
%     for i = 1: numel(LUE.amountfactors),
%         %fac(i) = rem(fac(i)-1,LUE.amountfactors(i)) + 1;
%         astring = [astring, sprintf('%s: %s | ',LUE.namesfinal{i},LUE.levelnames{i}{fac(i)})];
%     end; astring = astring(1:end-3);
    
    
    if runde > 1,
        
        holder   = {get(hand.listbox1,'String')};
        for i = size(holder{1},1) : -1 : 1, if isempty(deblank(holder{1}(i,:))), holder{1}(i,:) = []; end; end;
        astring2 = 'LUE.files(';
        for i = 1 : numel(LUE.amountfactors),
            %astring2 = sprintf('%s%d,',astring2,LUE.fac(i));
            astring2 = sprintf('%s%d,',astring2,LUE.fac(i));
        end;
        astring2 = sprintf('%s) = holder;',astring2(1:end-1));
        eval(astring2);
        
        set(hand.listbox1,'String','');
        
        for i = 1 : size(holder{1},1),
            LUE.final(end+1,1) = {LUE.fac(1)};
            for j = 2 : size(LUE.fac,2),
                LUE.final(end,j) = {LUE.fac(j)};
            end;
            LUE.final(end,end) = {holder{1}(i,:)};
        end;
        %LUE.final{end+1:size(holder{1},1),1:length(LUE.amountfactors)} = {repmat(LUE.fac,size(holder{1},1),1)};
        
        
    end;
    
    if runde <= prod(LUE.amountfactors),
        
    astring = [];

    fac = fullfact(LUE.amountfactors);
    fac = fac(runde,:);
    
    for i = 1: numel(LUE.amountfactors),
        %fac(i) = rem(fac(i)-1,LUE.amountfactors(i)) + 1;
        astring = [astring, sprintf('%s: %s | ',LUE.namesfinal{i},LUE.levelnames{i}{fac(i)})];
    end; astring = astring(1:end-3);
        
        set(hand.listbox1,'Visible','on');
        set(hand.pushbutton2,'Visible','on');
        set(hand.infotext,'String',['Nifti-Files for: ',astring]);
        
    elseif runde == prod(LUE.amountfactors) + 1,
        
        set(hand.infotext,'String','Overview');
        set(hand.uitable1,'Visible','on');
        set(hand.pushbutton1,'String','Perform Analysis');
        
        mytable = reshape(LUE.files,prod(LUE.amountfactors(1:ceil(end/2))),prod(LUE.amountfactors(ceil(end/2)+1:end)));
        
        for i = 1:size(mytable,1),
            for j = 1:size(mytable,2),
                mytables(i,j) = size(mytable{i,j},1);
            end;
        end; 
        
        set(hand.uitable1,'Data',mytables);
               
        fullc   = fullfact(LUE.amountfactors);
                
        repfullc = [];
        for i = 1 : size(fullc,1),
            repfullc{i} = fullc(i,:);
        end;
        astring4 = 'repfullc = reshape(repfullc,';
        for i = 1:numel(LUE.amountfactors),
            astring4 = sprintf('%s%d,',astring4,LUE.amountfactors(i));
        end; astring4 = sprintf('%s);',astring4(1:end-1));
        eval(astring4);
        repfullc = reshape(repfullc,prod(LUE.amountfactors(1:ceil(end/2))),prod(LUE.amountfactors(ceil(end/2)+1:end)));
        
        rownames = cell(1,size(mytables,1));
        colnames = cell(1,size(mytables,2));
        
        % rows
        vergl = cell2mat(repfullc(:,1));
        holder = [];
        for j = 1 : size(vergl,2),
            if length(unique(vergl(:,j))) > 1
                holder = [holder, j];
            end;
        end;
        
        rws = fullfact(LUE.amountfactors(holder));
        for h = 1:length(rownames),
            for j = 1 : length(holder),
                rownames{h} = sprintf('%s%s: %s | ',rownames{h},LUE.namesfinal{j}, LUE.levelnames{holder(j)}{rws(h,j)}); %rws(h,j));
            end;
            rownames{h} = rownames{h}(1:end-3);
        end;
        
        % cols
        vergl = cell2mat(repfullc(1,:)');
        holder = [];
        for j = 1 : size(vergl,2),
            if length(unique(vergl(:,j))) > 1
                holder = [holder, j];
            end;
        end;
        
        cws = fullfact(LUE.amountfactors(holder));
        for h = 1:length(colnames),
            for j = 1 : length(holder),
                colnames{h} = sprintf('%s%s: %s | ',colnames{h},LUE.namesfinal{holder(j)}, LUE.levelnames{holder(j)}{rws(h,j)});%cws(h,(j)));
            end;
            colnames{h} = colnames{h}(1:end-3);
        end;
        
        set(hand.uitable1,'columnName',colnames);
        set(hand.uitable1,'rowName',rownames);
                
    else,
        
        delete(hand.figure1);
        
        czc           = feature('numCores'); if czc > 12, czc = 12; end;
        %opft          = inputdlg(sprintf(sprintf('On how many cpus should be calculated on?\nMax of "%d" cores detected on your system.\nPlease be cautious with your selection. If you\nare not sure chose "1".\n',czc),czc),'Multi-Core-Calculation',1,{'1'});
        opft          = 1;
        
        %fcorr         = inputdlg(sprintf(sprintf('Choose a threshold for a FDR-Correction.\nAll results will be saved uncorrected\nand with the inserted correction threshold.\n',feature('numCores')),feature('numCores')),'FDR-Correction',1,{'.05'});
        fcorr         = inputdlg(sprintf(sprintf('Choose a nominal threshold for a FDR-Correction.\nAll results will be saved uncorrected\nand with the inserted correction threshold.\nFor more information please see the manual.\n',feature('numCores')),feature('numCores')),'FDR-Correction',1,{'5'});
        
        %fdrtype = questdlg('What FDR method should be used?','FDR method','conventional (conservative)','step-up (liberal)','conventional (conservative)');
        %if isequal(fdrtype,'conventional (conservative)'), LUE.fdrtype = 1; else, LUE.fdrtype = 2; end;
        LUE.fdrtype   = 1;
        
        %fdrcheck= questdlg(sprintf('What FDR correction should be performed?\nThe standard FDR correction will be too conservative but will be\ncalculated faster. The adapted FDR correction will result in more\naccurate results but needs longer to be calculated.\nFor details please see the manual.'),'FDR method','Standard','Adapted','Adapted');
        LUE.fdrcheck = 2;%if isequal(fdrtype,'Standard'), LUE.fdrcheck = 1; else, LUE.fdrcheck = 2; end;
        
        %if str2num(fcorr{1})<2, warndlg('The nominal threshold has to be at least 2. It is set to now.','Attention'); fcorr = {'2'}; end;
        
        %aholder       = questdlg('Should extensive Posthoc Tests be performed?','Extensive Posthoc Tests','Yes','No','Yes');
        aholder      = 'No';
        if isequal(aholder,'Yes'), LUE.praepost = 1; else, LUE.praepost = 0; end;
            
        LUE.now       = datestr(now,'yyyymmdd_HHMMSS');
        LUE.resultdir = uigetdir(pwd,'Choose a directory for results');
        LUE.resultdir = fullfile(LUE.resultdir,LUE.now); 
        LUE.core      = opft; %str2num(opft{1});
        LUE.fcorr     = str2num(fcorr{1});
        
        nix_perform;
    end;
    
    try, set(hand.gen,'string',num2str(runde));end;
    try, LUE.fac = fac; end;
end;

function data_zurechtsortieren,
global LUE

LUE.final = cell(0,0);
for i = 1 : size(LUE.within,1),
    LUE.final(i,:) = [num2cell([LUE.groups(i,:), LUE.within(i,:)]), {LUE.aim(i,:)}];
end;
LUE.namesfinal = {LUE.btnames{:},LUE.wtnamesxls{:},'Files'};

