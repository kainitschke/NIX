function nix_pinteraction(lauf,acoords,varargin)
% function for statistically test interaction between regions, lesion and
% within tests
% kai.nitschke@uniklinik-freiburg.de

if ~isempty(varargin)&(lauf~=5 & lauf ~= 6),
    lauf = varargin{1};
    acoords = varargin{2};
end;

%% GUI erstellen 1. GUI
if lauf == 0,
    
    global hans LUE
    
    fest = [600, 300];
    
    if LUE.type == 1,
        acol = prod(LUE.fac(1:end-1))+length(LUE.fac)-1;
        if (acol+2)*30 > fest(2),   fest(2) = (acol+2)*30;          end;
    elseif LUE.type == 2,
        if LUE.btall > 1, fest(2) = fest(2) + (LUE.btall - 2) * 30; end;
        if LUE.wtall > 2, fest(2) = fest(2) + (LUE.wtall - 2) * 30; end;
    end;
    
    res  = nix_res(100);
    
    %hans.voxelan1.fh = figure('Position',[res(1)/2-fest(1)/2 res(2)/2-fest(2)/2 fest(1) fest(2)],'Resize','off','MenuBar','none','Color',[.9 .9 .9],'Name','Specify for Interaction test','NumberTitle','off','visible','off');
    hans.voxelan1.fh = figure('Position',[res(1)+(diff(res([1,3]))+1)/2-fest(1)/2,  res(2) + (diff(res([2,4]))+1)/2-fest(2)/2,  fest(1),  fest(2)],'Resize','off','MenuBar','none','Color',[.9 .9 .9],'Name','Specify for Interaction test','NumberTitle','off','visible','off');
    set(hans.voxelan1.fh,'CloseRequestFcn','global hans; set(hans.voxelan1.fh,''Visible'',''off'')'); % schließen verhindern
    hans.voxelan1.text1 = uicontrol('Parent',hans.voxelan1.fh,'style','text','String','How many voxels:','position',[10 fest(2)-50 400 20],'HorizontalAlignment','left','BackgroundColor',[.9 .9 .9]);
    hans.voxelan1.editvoxel = uicontrol('parent',hans.voxelan1.fh,'style','edit','string','2','Position',[160 fest(2)-48 80 20],'BackgroundColor',[1,1,1]);
    set(hans.voxelan1.editvoxel,'Callback','pause(.0001); nix_pinteraction(3);');
    
    
    if LUE.type == 1,
        
        hans.voxelan1.text2 = uicontrol('Parent',hans.voxelan1.fh,'style','text','String','Factors and levels for analysis:','position',[10 fest(2)-80 400 20],'HorizontalAlignment','left','BackgroundColor',[.9 .9 .9]);
        
        for i = 1 : length(LUE.namesfinal)-1,
            
            eval(sprintf('hans.voxelan1.check%d = uicontrol(''Parent'',hans.voxelan1.fh,''style'',''checkbox'',''String'',''%s'',''position'',[20 fest(2)-50-30*(1+%d) 400 20],''HorizontalAlignment'',''left'',''BackgroundColor'',[.9 .9 .9],''Value'',1);',1+sum(LUE.fac(1:i-1))+length(LUE.fac(1:i-1)),LUE.namesfinal{i},1+sum(LUE.fac(1:i-1))+length(LUE.fac(1:i-1))));
            
            for j = 1:LUE.fac(i),
                eval(sprintf('hans.voxelan1.check%d = uicontrol(''Parent'',hans.voxelan1.fh,''style'',''checkbox'',''String'',''%s'',''position'',[20+20 fest(2)-50-30*(1+%d) 400 20],''HorizontalAlignment'',''left'',''BackgroundColor'',[.9 .9 .9],''Value'',1);',1+j+sum(LUE.fac(1:i-1))+length(LUE.fac(1:i-1)),LUE.levelnames{i}{j},1+j+sum(LUE.fac(1:i-1))+length(LUE.fac(1:i-1))));
            end;
        end;
        
    elseif LUE.type == 2
        
        if (LUE.btall > 1) | (LUE.wtall > 2),
            hans.voxelan1.text2 = uicontrol('Parent',hans.voxelan1.fh,'style','text','String','(Additional) Factors for analysis:','position',[10 fest(2)-80 400 20],'HorizontalAlignment','left','BackgroundColor',[.9 .9 .9]);
        end;
        i = 0;
        if LUE.btall > 1,
            for i = 1 : length(LUE.btnames),
                eval(sprintf('hans.voxelan1.check%d = uicontrol(''Parent'',hans.voxelan1.fh,''style'',''checkbox'',''String'',''%s'',''position'',[20 fest(2)-50-30*(1+%d) 400 20],''HorizontalAlignment'',''left'',''BackgroundColor'',[.9 .9 .9],''Value'',1);',i,LUE.btnames{i},i));
                for j = 1 : length(LUE.levelnames{1}),
                    eval(sprintf('hans.voxelan1.check%d = uicontrol(''Parent'',hans.voxelan1.fh,''style'',''checkbox'',''String'',''%s'',''position'',[20+20 fest(2)-50-30*(1+%d) 400 20],''HorizontalAlignment'',''left'',''BackgroundColor'',[.9 .9 .9],''Value'',1);',i+j,LUE.levelnames{1}{j},i+j));
                end;
            end;
        end;
        if i > 0, i = i + j; end;
        
        if (length(LUE.wtnames) > 1) | (LUE.wtall > 2),
            hans.voxelan1.text3 = uicontrol('Parent',hans.voxelan1.fh,'style','text','String','Which Within Tests:','position',[10 fest(2)-50-30*(2+i) 400 20],'HorizontalAlignment','left','BackgroundColor',[.9 .9 .9]);
            runner = 0;
            if length(LUE.wtnames) > 1, % multiple WTs
                for m = 1:LUE.wt(2),
                    for k = 1:LUE.wt(1),
                        runner = runner + 1;
                        eval(sprintf('hans.voxelan1.check%d%d = uicontrol(''Parent'',hans.voxelan1.fh,''style'',''checkbox'',''String'',''%s'',''position'',[20+80*(m-1) fest(2)-50-30*(2+i+k) 400 20],''HorizontalAlignment'',''left'',''BackgroundColor'',[.9 .9 .9],''Value'',1);',k,m,LUE.wtnamesxls{runner}));
                    end;
                end;
            else,                      % only one WT
                for m = 1:LUE.wt(1),
                    runner = runner + 1;
                    eval(sprintf('hans.voxelan1.check%d = uicontrol(''Parent'',hans.voxelan1.fh,''style'',''checkbox'',''String'',''%s'',''position'',[20 fest(2)-50-30*(2+i+m) 400 20],''HorizontalAlignment'',''left'',''BackgroundColor'',[.9 .9 .9],''Value'',1);',m+i,LUE.wtnamesxls{runner}));
                end;
            end;
        end;
        
    end;
    
    % Coordinate Table
    e = {'    ','    ','    ';'    ','    ','    '};
    hans.voxelan1.table    = uitable(hans.voxelan1.fh,'Data',e,'Columnwidth',{60},'ColumnName',{'x','y','z'},'ColumnEditable',[true,true,true],'Position',[320, 60, fest(1)/2-20*2, fest(2)-20*4]);
    hans.voxelan1.val = 2;
    
    % Calc Button
    hans.voxelan1.continue = uicontrol('Parent',hans.voxelan1.fh,'style','pushbutton','string','Calculate','position',[320 + fest(1)/2-20*2 - 80, 20, 80, 25]);
    set(hans.voxelan1.continue, 'Callback','nix_pinteraction(4);');
    
    % Write Button
    hans.voxelan1.write    = uicontrol('Parent',hans.voxelan1.fh,'style','pushbutton','string','Save Stats','position',[320, 20, 80, 25]);
    set(hans.voxelan1.write, 'Callback','nix_write_mult_coord;');
    
    
    %% initialer aufruf uns sichbarmachung 1. GUI
elseif lauf == 2,
    
    global hans
    set(hans.voxelan1.fh,'Visible','on');
    
    %% sichtbarmachung 2. GUI
elseif lauf == 3,
    
    global hans
    
    drawnow
    aimrows = str2num(get(hans.voxelan1.editvoxel,'String'));
    
    if ~isempty(aimrows),
        dat     = get(hans.voxelan1.table, 'Data');
        if ~isequal(aimrows, size(dat,1)),
            if aimrows > size(dat,1),
                for i = size(dat,1) + 1 : aimrows,
                    dat(i,:) = {'    '};
                end;
            else,
                dat = dat(1:aimrows,:);
            end;
            set(hans.voxelan1.table, 'Data', dat);
        end;
    end;
    
    %% Daten einlesen
elseif lauf == 4,
    
    global hans LUE
    
    bt = 0;
    wt = 0;
    
    dat = get(hans.voxelan1.table, 'Data');
    
    for i = 1:size(dat,1),
        for j = 1:size(dat,2),
            if strfind(dat{i,j},','),
                weg = strfind(dat{i,j},',');
                dat{i,j}(weg) = '.';
            end;
        end;
    end;
    
    try,
        dat = str2double(dat);
    catch,
        errordlg('You made a mistake entering the data','An error occured');
        return;
    end;
    
    if any(any(isnan(dat))),
        errordlg('No coordinate is allowed emtpy','An error occured');
        return;
    end;
    
    if LUE.type == 1,
        
       grps = ones(size(LUE.faccomb,1),1);
        
       for i = 1 : length(LUE.namesfinal)-1,
            
            eval(sprintf('weg = get(hans.voxelan1.check%d, ''Value'');',1+sum(LUE.fac(1:i-1))+length(LUE.fac(1:i-1))));
            if weg == 1,   afactor(i) = 1;
            else,          afactor(i) = 0;  end;
            
            weg  = [];
%             weg2 = zeros(size(LUE.faccomb,1),1);
            for j = 1:LUE.fac(i),
                eval(sprintf('weg(j) = get(hans.voxelan1.check%d, ''Value'');',1+j+sum(LUE.fac(1:i-1))+length(LUE.fac(1:i-1))));
                if ~weg(j),
                     %weg2 = weg2 + ([LUE.faccomb(:,i)]' == j);
                     weg2       = find(LUE.faccomb(:,i) == j);
                     grps(weg2) = 0;
                end;
            end;
            if (afactor(i)==1)&(sum(weg) < 2), 
                errordlg('You chose at least one factor with fewer than 2 levels. That does not make sense. Either choose at least 2 levels for every factor you want to model or uncheck the factor itself.','Error in Input');
                return;
            end;
            %grps = grps .* weg2;
        end; 
        
        %adata = LUE.final(find(grps),1:end-1);

        nix_pinteraction(6,dat,grps,afactor);
        
    elseif LUE.type == 2,
        
        grps = [];
        
        if LUE.btall > 1,
            try, if get(hans.voxelan1.check1,'Value'),
                    bt = 1;
                    for j = 1:length(LUE.levelnames{1}),
                        eval(sprintf('grps(j) = get(hans.voxelan1.check%d,''Value'');',1+j));
                    end;
                    %grps = find(grps);
                else,
                    for j = 1:length(LUE.levelnames{1}),
                        eval(sprintf('set(hans.voxelan1.check%d,''Value'',0);',1+j));
                    end;
                end; end;
        end;
        
        if (length(LUE.wtnames) > 1) | (LUE.wtall > 2),
            wtgrid = [];
            if length(LUE.wtnames) > 1,
                wt = 2;
                for k = 1:LUE.wt(1),
                    for m = 1:LUE.wt(2),
                        eval(sprintf('wtgrid(k,m) = get(hans.voxelan1.check%d%d, ''Value'');',k,m));
                    end;
                end;
                
                if sum(sum(wtgrid)) < 2,
                    errordlg('You have to choose at least two WT-columns','An error occured');
                    return;
                end;
                
                geg = zeros(LUE.wt(1),1);
                for i = 1:LUE.wt(2),
                    if sum(geg)==0,
                        geg = wtgrid(:,i);
                    elseif ~isequal(geg,wtgrid(:,i))&(sum(wtgrid(:,i))~=0),
                        errordlg('The chosen within tests are invalid. Each row has either to be completely unticked or has to be the same tick pattern as all other "ticked" rows.','An error occured');
                        return;
                    end;
                end;
                
            else,
                if LUE.btall > 1, i = 1+length(LUE.levelnames{1}); else, i = 0; end;
                for m = 1:LUE.wt(1),
                    eval(sprintf('wtgrid(m) = get(hans.voxelan1.check%d, ''Value'');',m+i));
                end;
            end;
            %     elseif LUE.wtall > 2,
            %         wt = 1;
            %         disp('EINLESEN FEHLT');
        else,
            wtgrid = [1 1];
        end;
        
        nix_pinteraction(5,dat,bt,wt,wtgrid,grps);
    end;
    
    
    
    %% Rechnen
elseif lauf == 5, % Posthoc for ordinal data
    
    global inn LUE hans
    
    hans.SaveMultCond = [];
    
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n')];
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n_______________________________________________________')];
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n')];
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n  Calculating Interaction effects for multiple voxels')];
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n_______________________________________________________')];
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n\n\n')];
    
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('------------- Voxel Comparison -------------\n')];
    
    %% change in real coords
    for i = 1:size(acoords,1),
        [acoords(i,:), aact(i)] = nix_dist(acoords(i,1),acoords(i,2),acoords(i,3));
        bts(:,i)                = inn.working(:,aact(i));
        aclustch(i)             = 1;
    end;
    nix_calc_pint_results(varargin{1},varargin{2},varargin{3},varargin{4},bts,acoords,aclustch);
    
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n\n\n------------- Cluster Comparison -------------\n')];
    
    cbts = []; aclustch = [];
    for i = 1:size(acoords,1),
        [acoords(i,:), aact(i)] = nix_dist(acoords(i,1),acoords(i,2),acoords(i,3));
        
        actcluster = hans.clust(hans.sel).clust(aact(i));
        if ~isequal(actcluster,0),
            weg         = find(hans.clust(hans.sel).clust == actcluster);
            cbts(:,i)   = mode(inn.working(:,weg)*-1,2)*-1;
            aclustch(i) = 1;
            
        else,
            cbts(:,i)   = bts(:,i);
            aclustch(i) = 0;
        end;
    end;
    nix_calc_pint_results(varargin{1},varargin{2},varargin{3},varargin{4},cbts,acoords,aclustch);
    
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n%s\n',repmat('-',120,1))];
    
    clc;
    disp(LUE.saveposthoc);
    disp(hans.SaveMultCond);
    
elseif lauf == 6, %post hoc for contigency table
    
    global LUE inn hans
    
    hans.SaveMultCond = [];
    clc;
    
    grps     = varargin{1};
    afactor  = varargin{2};
    
%     afaccomb = unique([LUE.faccomb(:,find(afactor)), LUE.faccomb(:,end)], 'rows'); %repmat(LUE.faccomb(:,end),1, size(acoords,1))];
%     for i = 2:size(acoords,1),
%         afaccomb = [repmat(afaccomb,LUE.lesiontypes,1),   reshape(repmat([1:LUE.lesiontypes],size(afaccomb,1),1),size(afaccomb,1)*LUE.lesiontypes,1)];
%     end;
%     
%     fdata = cell2mat(LUE.final(:,find(afactor)));
    
    
%     afaccomb = [LUE.faccomb(:,find(afactor)),  LUE.faccomb(:,end)];
%     weg      = afaccomb; %unique(afaccomb,'rows');
%     while size(weg,1)>0,
%         sucher = ones
%         
%     end;
    
    coltable = []; savetables = [];
    for i = 1 : size(acoords,1), % get tables for voxels and discard non-chosen groups
        [~, aact]     = nix_dist(acoords(i,1),acoords(i,2),acoords(i,3));
        atable        = inn.working(:,aact);
        weg           = find(grps == 0);
        atable(weg)   = 0;
        savetables{i} = atable;
        
        mytable       = reshape(atable,[LUE.fac]); % ,LUE.lesiontypes
        for j = length(LUE.fac)-1 : -1 : 1,
            if ~any(find(afactor)==j), % not modelled dimsension 
                mytable = sum(mytable,j);
            end
        end;
        mytable       = reshape(mytable,[LUE.fac(find(afactor)), LUE.lesiontypes]);
        %savetables{i} = mytable;
        
        trenner       = repmat(':,',1,length(find(afactor))+1);
        eval(sprintf('coltable(%s,%d) = mytable;',trenner(1:end-1),i))
    end;
    
    aeffects = ff2n(length(size(coltable)));
    weg      = find(sum(aeffects,2) > 1);
    aeffects = aeffects(weg,:);
    
    reslog   = nix_contingency(coltable,aeffects);
    
    %% Output
    
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n')];
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n_______________________________________________________')];
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n')];
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n  Calculating Interaction effects for multiple voxels')];
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n_______________________________________________________')];
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n\n\n')];
    
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('------------- Voxel Comparison -------------\n')];
    for i = 1:size(acoords,1),
        hans.SaveMultCond = [hans.SaveMultCond, nix_cont_desc(acoords(i,:),savetables{i},afactor,LUE.namesfinal,LUE.levelnames,LUE.faccomb,LUE.fac,aeffects,LUE.lesiontypes)];
    end;
    
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n\n------------- Statistics -------------\n\n')];
    hans.SaveMultCond = [hans.SaveMultCond, nix_cont_res(reslog,afactor,LUE.namesfinal,aeffects)];
    
    fprintf('%s\n',hans.SaveMultCond);
    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nix_calc_pint_results(varargin)

global LUE hans

cbt      = varargin{1};
cwt      = varargin{2};

bts      = varargin{5};
acoords  = varargin{6};
aclustch = varargin{7};

%% group size evaluation
if cbt == 1,
    bts    = [bts, LUE.groups];
    grsize = length(LUE.btnames{1}); for i = 1:LUE.bt, if length(LUE.levelnames{1}{i}) > grsize, grsize = length(LUE.levelnames{1}{i}); end; end;% for max size later
    agr    = LUE.bt(1);
else,
    bts = [bts, ones(size(bts,1),1)];
    agr = 1;
end;

wtgrid = varargin{3};

if cwt == 2,
    adim(1) = max(sum(wtgrid,1));
    adim(2) = max(sum(wtgrid,2));
elseif cwt == 1,
    disp('FEHLT weiter');
else,
    adim(1) = sum(wtgrid);
end;

combs  = unique(bts(:,1:end-1),'rows');

if ((length(LUE.imgperm) ^ size(acoords,1)) > size(combs,1)) | size(acoords,1) > 2 | cbt == 1, % not enough groups
    testmode = [2,1]; % 1 = keine überschneidungen
elseif (length(adim)>1) & (min(adim)>1),
    testmode = [2,3]; % zwei within-variablen gewählt
elseif length(adim)>1,
    if adim(1) > 1,
        testmode = [1,1];
    else,
        testmode = [1,2];
    end;
else, % enough groups
    testmode = [1,0];
end;
btn = zeros(size(bts,1),1);

for i = 1:size(combs,1),
    agroup = ones(size(bts,1),1);
    for j = 1:size(combs,2),
        agroup = agroup .* (combs(i,j) == bts(:,j));
    end;
    btn(find(agroup)) = i;
    for k = 1:LUE.bt,
        agroupk = agroup .* (LUE.groups == k);
        n(i,k) = length(find(agroupk));
    end;
end;

if cbt == 1, ncollapse = n; else, ncollapse = sum(n,2); end;

if (testmode(1)==1)&(any(ncollapse < 2)),
    testmode = [2,2]; % 2 = zu wenig überschneidung
end;

%% collect infos and calculate

nraus = [];

if testmode(1) == 2,
    if cbt == 1,
        adata = [btn, LUE.groups, LUE.within(:,find(wtgrid))];
    else,
        adata = [btn, ones(size(btn,1),1), LUE.within(:,find(wtgrid))];
    end;
    
    if any(any(ncollapse < 2)),
        for i = 1:size(ncollapse,1),
            if any(any(ncollapse(i,:)<2)),
                adata(find((adata(:,1) == i)),:) = [];
                nraus = [nraus, i];
            end;
        end;
    end;
    
    if cbt == 1,
        grps = varargin{4};
        for j = 1:length(grps),
            if ~grps(j),
                weg = find(adata(:,2)==j);
                adata(weg,:) = [];
            end;
        end;
    end;
    
    switch 1,
        case cwt < 2,% | ((cwt == 2) & (any(adim) < 2)),
            try, result = nix_f2ldf1(adata); end;
            
        case cwt == 2,
            try, result = nix_f1ldf2([adata(:,1),adata(:,3:end)],adim(1),adim(2)); end;
            
        otherwise,
            error(sprintf('''This shouldn''t happen. Strange ...'''));
    end;
    
elseif testmode(1) == 1,
    adata = [bts(:,1:2), LUE.within(:,find(wtgrid))];
    result = nix_f2ldf1(adata);
end;

%% print results
abstand = 6;
if testmode(1) == 1,
    abstand = abstand + 2 + 8 + 1;
end;
if cbt == 1,
    abstand = abstand + length(LUE.btnames{1}) + 1;
end;
if cwt == 0 & length(LUE.wt < 2),
    abstand = abstand + length(LUE.wtnames{1}) + 1;
elseif min(adim) == 1,
    if adim(1) ~= 1,
        abstand = abstand + 1 + length(LUE.wtnames{1});
    else,
        abstand = abstand + 1 + length(LUE.wtnames{2});
    end;
else,
    abstand = abstand + 1 + length(LUE.wtnames{1}) + length(LUE.wtnames{2});
end;

%fprintf('\n'); for i = 1:size(acoords,1), fprintf('Voxel#%d:',i); fprintf('%7.1f',acoords(i,:)); fprintf('  |  '); end; fprintf('\b\b\b\n');
hans.SaveMultCond = [hans.SaveMultCond, sprintf('\nCOORDINATES:')];
hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n%14s%7s%7s','x','y','z')];
for i = 1:size(acoords,1),
    hans.SaveMultCond =  [hans.SaveMultCond, sprintf('\nVoxel#%02d:',i)];
    hans.SaveMultCond =  [hans.SaveMultCond, sprintf('%7.1f',acoords(i,:))];
    if aclustch(i) == 0, hans.SaveMultCond = [hans.SaveMultCond, sprintf(' # ATTENTION: No cluster found; direct voxel used')]; end;
end;
hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n')];

hans.SaveMultCond = [hans.SaveMultCond, sprintf('\nDESCRIPTIVES:')];
for i = 1 : length(LUE.imgperm) ^ size(acoords,1),
    for j = 1:size(acoords,1),
        ac(i,j) = mod(floor((i-1)/length(LUE.imgperm)^(size(acoords,1)-j)),length(LUE.imgperm))+1;
    end;
end;
hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n')];
if cbt==1,
    hans.SaveMultCond = [hans.SaveMultCond, sprintf(sprintf('%%%ds     ',grsize),LUE.btnames{1})];
end;
for i = 1:size(acoords,1),
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('   Voxel#%02d',i)];
end;
hans.SaveMultCond = [hans.SaveMultCond, sprintf('     Number')];
hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n%s',repmat('-',1,3+8*size(acoords,1)+3*(size(acoords,1)-1)+11))];
if cbt==1,
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('%s',repmat('-',1,grsize+5))];
end;
numbcollect = [];
for g = 1:agr,
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n')];
    if cbt==1,
        hans.SaveMultCond = [hans.SaveMultCond, sprintf(sprintf('%%%ds     ',grsize),LUE.levelnames{1}{g})];
    end;
    for i = 1:size(ac,1),
        if (cbt==1) & (i~=1),
            hans.SaveMultCond = [hans.SaveMultCond, sprintf(sprintf('%%%ds',grsize+5),'')];
        end;
        numbs = ones(size(combs,1),1);
        for j = 1 : size(ac,2),
            hans.SaveMultCond = [hans.SaveMultCond, sprintf('%11d',round(LUE.imgperm(ac(i,j))))];
            numbs = numbs .* (LUE.imgperm(ac(i,j)) == combs(:,j));
        end;
        
        vstat = 1;
        for j = 1:length(nraus),
            if isequal(LUE.imgperm(ac(i,:)),combs(nraus(j),:)'),
                vstat = 0;
            end;
        end; % check whether excluded
        
        numbs = find(numbs);
        if isempty(numbs), numbs = 0; else, numbs = ncollapse(numbs,g); end;
        hans.SaveMultCond = [hans.SaveMultCond, sprintf('% 11d',numbs)];
        %if any(ncollapse(:,g)) < 2, fprintf('   # excluded'); end;
        if (numbs<2)|(vstat==0),
            hans.SaveMultCond = [hans.SaveMultCond, sprintf('   # excluded')];
        elseif cbt == 1,
            if any(ismember(find(~grps),g)),
                hans.SaveMultCond = [hans.SaveMultCond, sprintf('   # excluded')];
            end;
        end;
        hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n')];
        numbcollect       = [numbcollect, numbs];
    end;
end;

hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n')];
hans.SaveMultCond = [hans.SaveMultCond, sprintf('\nTEST DESCRIPTION')];
if any(aclustch == 0),
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\nATTENTION: Some or all voxels did not belong to any cluster. Therefore, the direct voxels were used.')];
end;
if testmode(1) == 1,
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\nVoxels were modelled as two independent factors')];
else,
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\nVoxels were modelled only as one factor because:')];
    if size(acoords,1) > 2,   hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n- you entered more than two voxels')]; end;
    if any(any(ncollapse<2)), hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n- at least one group had fewer than 2 people')]; end;
    if cbt==1,                hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n- you also wanted a between factor to be considered')]; end;
    if cwt==2,                hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n- you also wanted a second within factor to be considered')]; end;
    if any(numbcollect<2),    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n- There were voxels that had to be excluded due to insufficient sample sizes, therefore, no full crossed design was available')]; end;
end;
hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n')];

hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n')];
hans.SaveMultCond = [hans.SaveMultCond, sprintf('RESULTS')];
nix_print_results([],[],[],[],2,abstand);
nix_print_results([],[],[],[],0,abstand); % Header drucken
nix_print_results([],[],[],[],2,abstand);

if testmode(1) == 2,
    nix_print_results('Voxel',result.anova.Fgr1,result.anova.dfgr1,result.anova.pgr1,1,abstand);
    
    if length(LUE.wtnames) == 1,
        nix_print_results(LUE.wtnames{1},result.anova.Ft1,result.anova.dft1,result.anova.pt1,1,abstand);
        nix_print_results(sprintf('Voxel:%s',LUE.wtnames{1}),result.anova.Fgr1t1,result.anova.dfgr1t1,result.anova.pgr1t1,1,abstand);
    else,
        if (length(adim)>1) & (min(adim)>1),
            nix_print_results(sprintf('%s',LUE.wtnames{1}),result.anova.Ft1,result.anova.dft1,result.anova.pt1,1,abstand);
            nix_print_results(sprintf('%s',LUE.wtnames{2}),result.anova.Ft2,result.anova.dft2,result.anova.pt2,1,abstand);
            nix_print_results(sprintf('Voxel:%s',LUE.wtnames{1}),result.anova.Fgr1t1,result.anova.dfgr1t1,result.anova.pgr1t1,1,abstand);
            nix_print_results(sprintf('Voxel:%s',LUE.wtnames{2}),result.anova.Fgr1t2,result.anova.dfgr1t2,result.anova.pgr1t2,1,abstand);
            nix_print_results(sprintf('%s:%s',LUE.wtnames{1},LUE.wtnames{2}),result.anova.Ft1t2,result.anova.dft1t2,result.anova.pt1t2,1,abstand);
            nix_print_results(sprintf('Voxel:%s:%s',LUE.wtnames{1},LUE.wtnames{2}),result.anova.Fgr1t1t2,result.anova.dfgr1t1t2,result.anova.pgr1t1t2,1,abstand);
            
        elseif adim(1)>1,
            nix_print_results(sprintf('%s',LUE.wtnames{1}),result.anova.Ft1,result.anova.dft1,result.anova.pt1,1,abstand);
            nix_print_results(sprintf('Voxel:%s',LUE.wtnames{1}),result.anova.Fgr1t1,result.anova.dfgr1t1,result.anova.pgr1t1,1,abstand);
            
        elseif adim(2)>1
            nix_print_results(sprintf('%s',LUE.wtnames{2}),result.anova.Ft2,result.anova.dft2,result.anova.pt2,1,abstand);
            nix_print_results(sprintf('Voxel:%s',LUE.wtnames{2}),result.anova.Fgr1t2,result.anova.dfgr1t2,result.anova.pgr1t2,1,abstand);
        else,
            disp('SOLLTE NICHT PASSIEREN');
            
        end;
    end;
    
    if cbt,
        nix_print_results(LUE.btnames{1},result.anova.Fgr2,result.anova.dfgr2,result.anova.pgr2,1,abstand);
        nix_print_results(sprintf('Voxel:%s',LUE.btnames{1}),result.anova.Fgr1gr2,result.anova.dfgr1gr2,result.anova.pgr1gr2,1,abstand);
        nix_print_results(sprintf('%s:%s',LUE.btnames{1},LUE.wtnames{1}),result.anova.Fgr2t1,result.anova.dfgr2t1,result.anova.pgr2t1,1,abstand);
        nix_print_results(sprintf('Voxel:%s:%s',LUE.btnames{1},LUE.wtnames{1}),result.anova.Fgr1gr2t1,result.anova.dfgr1gr2t1,result.anova.pgr1gr2t1,1,abstand);
    end;
    
elseif testmode(1) == 1,
    astr = LUE.wtnames{1};
    try, if adim(1) == 1, astr = LUE.wtnames{2}; end; end;
    nix_print_results('Voxel#01',result.anova.Fgr1,result.anova.dfgr1,result.anova.pgr1,1,abstand);
    nix_print_results('Voxel#02',result.anova.Fgr2,result.anova.dfgr2,result.anova.pgr2,1,abstand);
    nix_print_results(sprintf('%s',astr),result.anova.Ft1,result.anova.dft1,result.anova.pt1,1,abstand);
    nix_print_results('Voxel#01:Voxel#02',result.anova.Fgr1gr2,result.anova.dfgr1gr2,result.anova.pgr1gr2,1,abstand);
    nix_print_results(sprintf('Voxel#01:%s',astr),result.anova.Fgr1t1,result.anova.dfgr1t1,result.anova.pgr1t1,1,abstand);
    nix_print_results(sprintf('Voxel#02:%s',astr),result.anova.Fgr2t1,result.anova.dfgr2t1,result.anova.pgr2t1,1,abstand);
    nix_print_results(sprintf('Voxel#01:Voxel#02:%s',astr),result.anova.Fgr1gr2t1,result.anova.dfgr1gr2t1,result.anova.pgr1gr2t1,1,abstand);
end;

hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n_______________________________________________________')];

function nix_print_results(rowname,F,df,p,mode,abstand),
global hans
laengen = [10 8 10];
if     mode == 0,
    %eval(sprintf('fprintf(''\\n%%%ds%%%ds%%%ds%%%ds'',''Effect'',''F-value'',''df'',''p-value'');',abstand,laengen(1),laengen(2),laengen(3)));
    hans.SaveMultCond = [hans.SaveMultCond, eval(sprintf('sprintf(''\\n%%%ds%%%ds%%%ds%%%ds'',''Effect'',''F-value'',''df'',''p-value'');',abstand,laengen(1),laengen(2),laengen(3)))];
elseif mode == 1,
    %eval(sprintf('fprintf(''\\n%%%ds%%%d.3f%%%d.2f%%%d.5f'',rowname,F,df,p);',abstand,laengen(1),laengen(2),laengen(3)));
    hans.SaveMultCond = [hans.SaveMultCond, eval(sprintf('sprintf(''\\n%%%ds%%%d.3f%%%d.2f%%%d.5f'',rowname,F,df,p);',abstand,laengen(1),laengen(2),laengen(3)))];
    if (p < .1)&(p >= .05), hans.SaveMultCond = [hans.SaveMultCond, sprintf(' .')]; end;
    if p < .05,  hans.SaveMultCond = [hans.SaveMultCond, sprintf(' *')]; end;
    if p < .01,  hans.SaveMultCond = [hans.SaveMultCond, sprintf('*')];  end;
    if p < .001, hans.SaveMultCond = [hans.SaveMultCond, sprintf('*')];  end;
elseif mode == 2,
    hans.SaveMultCond = [hans.SaveMultCond, sprintf('\n%s',repmat('-',1,abstand + sum(laengen)))];
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function resultstr = nix_cont_desc(acoords,mytable,afactor,namesfinal,levelnames,faccomb,fac,aeffects,lesiontypes),

resultstr = sprintf('\n### Voxel:%7.1f% 7.1f% 7.1f\n',acoords);
resultstr = [resultstr, sprintf('\n')];

maxlength = 5; % find max length for output
for i = 1:length(namesfinal), 
    if length(namesfinal{i}) > maxlength, maxlength = length(namesfinal{i}); end;
    for j = 1:length(levelnames{i}),
        if length(levelnames{i}{j}) > maxlength, maxlength = length(levelnames{i}{j}); end;
    end;
end;
maxlength = maxlength + 2;

for i = 1 : length(namesfinal)-1,
    if afactor(i) == 1,
        %resultdir = sprintf('maxlength')
        eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',%s)];',maxlength,sprintf('namesfinal{%d}',i)));
    end;
end;

eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',maxlength,sprintf('ImgData')));
%eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',maxlength,sprintf('Voxel')));
eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',maxlength,sprintf('Number')));
resultstr = [resultstr, sprintf('\n%s\n',repmat('-',1,(length(find(afactor))+2)*maxlength))];

widx = [find(afactor), length(afactor)+1];
weg  = [fac(find(afactor)), lesiontypes];
allvalues           = zeros(prod(weg),1); %weg wieder frei
rowcomb             = unique(faccomb(:,widx),'rows');
for j = 1 : size(rowcomb,1),
    holder          = ones(size(faccomb,1),1);
    for i = 1 : length(widx),
        holder      = holder .* (faccomb(:,widx(i)) == rowcomb(j,i));
    end;
    allvalues(j)    = sum(mytable(find(holder)));
    
end; N = allvalues;

for i = 1:size(rowcomb,1),
    for j = 1:size(rowcomb,2),
        eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',maxlength,sprintf('%s',levelnames{widx(j)}{rowcomb(i,j)})));
    end;
    eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',maxlength,sprintf('%d',N(i))));
    resultstr = [resultstr, sprintf('\n')];
end;
resultstr = [resultstr, sprintf('%s\n',repmat('-',1,(length(find(afactor))+2)*maxlength))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resultstr = nix_cont_res(reslog,afactor,namesfinal,aeffect)

statlength = [10,6,8];

afactor    = [afactor, 1, 1];
namesfinal{end+1} = 'Voxel';

afac       = find(afactor);
maxlength  = length([namesfinal{afac}])+(length(afac)-1)*3;

resultstr = [];
resultstr = [resultstr, sprintf('\n%s\n',repmat('-',1,maxlength + sum(statlength)))];
eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',maxlength,sprintf('Effect')));
eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',statlength(1),sprintf('lrt')));
eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',statlength(2),sprintf('df')));
eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',statlength(3),sprintf('p')));
resultstr = [resultstr, sprintf('\n%s\n',repmat('-',1,maxlength + sum(statlength)))];

for i = 1 : length(reslog),
    aname = [];
    for j = 1 : size(aeffect,2),
        if aeffect(i,j) == 1,
            aname = sprintf('%s%s X ',aname,namesfinal{afac(j)});
        end;
    end;
    aname = aname(1:end-3);
    eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',maxlength,sprintf('%s',aname)));
    eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',statlength(1),sprintf('%1.3f',reslog(i).lrt)));
    eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',statlength(2),sprintf('%d',reslog(i).df)));
    eval(sprintf('resultstr = [resultstr, sprintf(''%%%ds'',''%s'')];',statlength(3),sprintf('%0.3f',reslog(i).p)));
    resultstr = [resultstr, sprintf('\n')];
end;
resultstr = [resultstr, sprintf('%s\n\n',repmat('-',1,maxlength + sum(statlength)))];