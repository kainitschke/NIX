function nix_result_posthoc(varargin),

global LUE hans inn

%% get varargins
amode   = varargin{1}; % 1 = descriptives | 2 = post-hoc
aeffect = varargin{2};
aact    = varargin{3};
sel     = varargin{4};
if LUE.type == 2,  adata   = varargin{5}; end;


%% Collecting

if LUE.type == 1,
    
    anames = {}; widx = []; nvals = cell(0,0);
    for i = 1 : length(LUE.namesfinal),
        if any(strfind(aeffect,LUE.namesfinal{i})),
            anames{end+1} = LUE.namesfinal{i};
            widx          = [widx, find(ismember(LUE.namesfinal,anames{end}))];
        end;
    end;
    weg = LUE.fac(widx); %start weg
    
    %             if any(strfind(aeffect,'ImgData')),
    %                 anames{end+1} = 'ImgData';
    %                 widx          = [widx, length(LUE.namesfinal)+1];
    %                 weg           = [weg,  LUE.lesiontypes];
    %             end;
    
    allvalues           = zeros(prod(weg),1); %weg wieder frei
    rowcomb             = unique(LUE.faccomb(:,widx),'rows');
    for j = 1 : size(rowcomb,1),
        holder          = ones(size(LUE.faccomb,1),1);
        for i = 1 : length(widx),
            holder      = holder .* (LUE.faccomb(:,widx(i)) == rowcomb(j,i));
        end;
        %if isempty(inn.working), inn = matfile(fullfile(hans.effectpath,LUE.niicontfile)); end;
        allvalues(j)    = sum(inn.working(find(holder),aact));
        
        for i = 1 : size(rowcomb,2),
            anvals{i,j} = sprintf('%s: %s',LUE.namesfinal{widx(i)},LUE.levelnames{widx(i)}{rowcomb(j,i)});
        end; fest = i;
        
        weg = find(~ismember([1:length(LUE.namesfinal)],widx));
        for i = 1 : length(weg),
            weg2        = sprintf('%s: ',LUE.namesfinal{weg(i)});
            weg2        = [weg2, sprintf('%s, ',LUE.levelnames{weg(i)}{:})];
            anvals{fest+i,j} = weg2(1:end-2);
        end;
        
    end; N = allvalues;
    
    
elseif LUE.type == 2,
    
    av            = LUE.ranks;
    av2           = LUE.within;
    
    avcols = cell(0,0);
    if any(strfind(aeffect,'ImgData')),   anames        = {'ImgData'}; else, anames=[]; end; %anames        = {'ImgData'};
    if LUE.bt ~= 1, anames{end+1} = LUE.btnames{1}; end;
    anames{end+1} = 'Within';
    
    awt = LUE.wt; awtnames = LUE.wtnames; if length(LUE.wt)==1, awt = [awt,1]; awtnames{2} = ''; end;
    %disp('HIER ACHTUNG!!!');
    if (~isempty(awtnames{2})) & any(strfind(aeffect,awtnames{2}))
        if any(strfind(aeffect,awtnames{1})),
            for nn = 1:LUE.wtall,
                avcols{nn} = nn;
            end;
        else,
            for nn = 1:awt(2),
                avcols{nn} = [(nn-1)*awt(2)+1 : nn*awt(2)];
            end;
        end;
    else,
        if any(strfind(aeffect,awtnames{1})),
            for mm = 1:awt(1),
                avcols{mm} = [mm : awt(1): LUE.wtall];
            end;
        else,
            avcols{1} = [1:LUE.wtall];
        end;
    end;
    
    if ~isequal(LUE.bt,1), if any(strfind(aeffect,LUE.btnames{1})),
            adata = [adata,LUE.groups];
        end; end;
    
    if isempty(adata), adata = ones(size(LUE.groups,1),1); end;
    
    [~,~,idx] = unique(adata,'rows');
    
    % Zusammenstellen
    hans.mittel = []; hans.SEs = []; hans.werte = []; hans.werteraw = []; hans.ausp = []; anvals = []; N = [];
    for j = 1:length(avcols),
        for i = unique(idx)',
            hans.werte{end+1}  = av(find(idx==i),avcols{j}); hans.werte{end} = reshape(hans.werte{end},size(hans.werte{end},1)*size(hans.werte{end},2),1);
            hans.werteraw{end+1}  = av2(find(idx==i),avcols{j}); hans.werteraw{end} = reshape(hans.werteraw{end},size(hans.werteraw{end},1)*size(hans.werteraw{end},2),1);
            
            nn = median(adata(find(idx==i),:));
            if any(strfind(aeffect,LUE.btnames{1})) & (any(strfind(aeffect,'ImgData'))),
                anvals{1,end+1} = sprintf('%s: %d',anames{1},round(nn(1)));
                anvals{2,end} = sprintf('%s: %s',LUE.btnames{1},LUE.levelnames{1}{round(nn(2))});
                
            elseif any(strfind(aeffect,LUE.btnames{1})),
                anvals{1,end+1} = sprintf('%s: %s',LUE.btnames{1},LUE.levelnames{1}{round(nn(1))});
            else,
                anvals{1,end+1} = sprintf('%s: %d',anames{1},round(nn(1)));
            end;
            hans.ausp{end+1,1} = nn(1);
            
            if length(nn) > 1,
                if ~isequal(LUE.btall,1),
                    %anvals{2,end} = sprintf('%s: %d',anames{2},nn(2));
                    anvals{2,end} = sprintf('%s: %s',anames{2},LUE.levelnames{1}{nn(2)});%nn(2));
                    hans.ausp{end,2} = nn(2);
                end;
                
            end;
            nm = sprintf('%s, ',LUE.wtnamesxls{avcols{j}}); nm = nm(1:end-2);
            anvals{length(nn)+1,end} = sprintf('%s: %s',anames{end},nm);
            N(end+1) = length(find(idx==i));
            hans.ausp{end,length(nn)+1} = avcols{j};
            hans.ausp{end,length(nn)+2} = median(hans.werte{end});%mean(adata(find(idx==i),:));
            hans.ausp{end,length(nn)+3} = N(end);
        end;
    end;
end;

%% barplots

if amode == 1,
    if LUE.type == 1,
        
        bar(allvalues,'FaceColor',[.8,.8,.8],'Parent',hans.raw);
        
        avor = 'lrt';  ylse = [];
        ylse = get(hans.raw,'YLim'); ylse = [0,ylse(2)+diff(ylse)/20]; %ylse = [0-diff(ylse)/20,ylse(2)];
        xlse = get(hans.raw,'XLim');
        
        set(hans.raw,'XtickLabel',[],'YLim',ylse);
        
        text(mean(xlse), ylse(2) + diff(ylse)/20,sprintf('%s = %1.3f, p_u_n_c_o_r_r = %0.3f',avor,hans.resultinf(sel).private.dat(aact),hans.plink(sel).private.dat(aact)),'Parent',hans.raw,'HorizontalAlignment','center','FontSize',16);
        text(xlse(1)-diff(xlse)/17,mean(ylse),'Number per Group','FontSize',18,'HorizontalAlignment','Center','Parent',hans.raw,'Rotation',90);
        
        LUE.astat = sprintf('%s = %1.3f, puncorr = %0.3f',avor,hans.resultinf(sel).private.dat(aact),hans.plink(sel).private.dat(aact));% for save Raw data
        
        %X-Achsebeschriftung
        for i = 1:size(anvals,2)
            for j = 1:size(anvals,1)
                text(i,ylse(1)*1.0-diff(ylse*1.15)*(1/800+(1/45)*(j-1)),anvals{j,i},'FontSize',9-1,'HorizontalAlignment','Center','Parent',hans.raw);
            end;
            t = text(i,ylse(1)+diff(ylse*1.15)/40,sprintf('n = %d',N(i)), 'Fontsize',8,'Parent',hans.raw,'HorizontalAlignment','Center','VerticalAlignment','middle');
        end;
        
    elseif LUE.type == 2,
        allvalues = [];
        for i = 1:length(hans.werte),
            allvalues = [allvalues; hans.werte{i}, repmat(i,length(hans.werte{i}),1)];
        end;
        
        ylse = [0,0]; % enthält die Maximalwerte für Achnenlimits
        lst = zeros(1,length(hans.mittel));
        
        cgroup = ones(1,length(hans.werte));
        if hans.efilestype(sel) == 2,
            fact = find(ismember(LUE.postnames(:,1),hans.efiles{sel}(13:end)));
            for i = 1:length(hans.werte),
                switch 1,
                    case ~isequal(LUE.btall,1),
                        if (hans.ausp{i,2} == LUE.postnames{fact,4}) & ((LUE.postnames{fact,3}==hans.ausp{i,3}) | (LUE.postnames{fact,2}==hans.ausp{i,3})) & ((LUE.postnames{fact,5}==hans.ausp{i,1})|(LUE.postnames{fact,6}==hans.ausp{i,1}))
                            cgroup(i) = 2;
                        end;
                        
                    otherwise,
                        if ((LUE.postnames{fact,3}==hans.ausp{i,2}) | (LUE.postnames{fact,2}==hans.ausp{i,2}))
                            cgroup(i) = 2;
                        end;
                end;
            end;
        end;
        
        b = boxplot(allvalues(:,1),allvalues(:,2),'colors',[0,0,0;1,0,0],'colorgroup',cgroup,'Parent',hans.raw);
        
        ylse = get(hans.raw,'YLim'); ylse = [0-diff(ylse)/20,ylse(2)];
        xlse = get(hans.raw,'XLim');
        set(hans.raw,'XTickLabel',' ','XTick',[1:length(hans.werte)],'TickDir','out','YGrid','on','box','off','YLim',ylse);
        % Beschriftung Graph
        if isequal(LUE.statkind,'wald'), avor = 'Q'; elseif isequal(LUE.statkind,'anova'), avor = 'F'; end;
        if any(ismember(inn.fdrvar,aact)),
            text(mean(xlse), ylse(2) + diff(ylse)/20,sprintf('%s = %1.3f, p_u_n_c_o_r_r = %0.6f, p_F_D_R_a_d_j = %1.6f',avor,abs(hans.resultinf(sel).private.dat(aact)),hans.plink(sel).private.dat(aact),hans.padjlink(sel).private.dat(aact)),'Parent',hans.raw,'HorizontalAlignment','center','FontSize',16);
            LUE.astat = sprintf('%s = %1.3f, puncorr = %0.6f, pFDRadj = %1.6f',avor,abs(hans.resultinf(sel).private.dat(aact)),hans.plink(sel).private.dat(aact),hans.padjlink(sel).private.dat(aact));% for save Raw data
        else,
            text(mean(xlse), ylse(2) + diff(ylse)/20,sprintf('%s = %1.3f, p_u_n_c_o_r_r = %0.6f, p_F_D_R_a_d_j = n.c.',avor,abs(hans.resultinf(sel).private.dat(aact)),hans.plink(sel).private.dat(aact)),'Parent',hans.raw,'HorizontalAlignment','center','FontSize',16);
            LUE.astat = sprintf('%s = %1.3f, puncorr = %0.6f, pFDRadj = n.c.',avor,abs(hans.resultinf(sel).private.dat(aact)),hans.plink(sel).private.dat(aact));% for save Raw data
        end;
        
        text(xlse(1)-diff(xlse)/17,mean(ylse),'Ranks','FontSize',18,'HorizontalAlignment','Center','Parent',hans.raw,'Rotation',90);
        
        %X-Achsebeschriftung
        for i = 1:size(anvals,2)
            for j = 1:size(anvals,1)
                text(i,ylse(1)*1.15-diff(ylse*1.15)*(1/35+(1/35)*(j-1)),anvals{j,i},'FontSize',9-1,'HorizontalAlignment','Center','Parent',hans.raw);
            end;
            t = text(i,ylse(1)+diff(ylse*1.15)/120,sprintf('n = %d',N(i)), 'Fontsize',8,'Parent',hans.raw,'HorizontalAlignment','Center','VerticalAlignment','middle');
        end;
    end;
end;

%% Descriptives

if amode == 2,
    if LUE.type == 2,
        
        %LUE.saveposthoc = [LUE.saveposthoc, sprintf('----------------------------------------\n')];
        LUE.saveposthoc = [LUE.saveposthoc, sprintf('------------- Descriptives -------------\n\n\n')];
        %LUE.saveposthoc = [LUE.saveposthoc, sprintf('----------------------------------------\n\n\n')];
        laengst = 0;
        for i = 1:length(anames), if length(anames{i})>laengst, laengst = length(anames{i}); end; end;
        for i = 1:length(LUE.wtnamesxls), if length(LUE.wtnamesxls{i})>laengst, laengst = length(LUE.wtnamesxls{i}); end; end;
        laengst = laengst + 4;
        LUE.saveposthoc = [LUE.saveposthoc, sprintf('%s\n',repmat('-',1,(laengst+1)*(length(anames)+4)))];
        %eval(sprintf('fprintf(''%% %ds '',''#'')',laengst));
        LUE.saveposthoc = [LUE.saveposthoc, eval(sprintf('sprintf(''%% %ds '',''#'')',laengst))];
        for j = 1:length(anames),
            LUE.saveposthoc = [LUE.saveposthoc, eval(sprintf('sprintf(''%% %ds '',anames{j})',laengst))];
        end;
        LUE.saveposthoc = [LUE.saveposthoc, eval(sprintf('sprintf(''%% %ds %% %ds %% %ds'',''Median'',''N'')',laengst,laengst,laengst))];
        LUE.saveposthoc = [LUE.saveposthoc, sprintf('\n%s\n',repmat('-',1,(laengst+1)*(length(anames)+4)))];
        for i = 1:size(hans.ausp,1),
            LUE.saveposthoc = [LUE.saveposthoc, eval(sprintf('sprintf(''%% %dd '',i)',laengst))];
            for j = 1:size(hans.ausp,2),
                avalx = sprintf('%d,',round(hans.ausp{i,j}));
                %try, if isequal('Within',anames{j}), avalx = sprintf('%s,',LUE.wtnamesxls{round(hans.ausp{i,j})}); end; end;
                try, if isequal('Within',anames{j}), avalx = sprintf('%s,',LUE.wtnamesxls{round(hans.ausp{i,j})}); end; end;
                LUE.saveposthoc = [LUE.saveposthoc, eval(sprintf('sprintf(''%% %ds '',avalx(1:end-1))',laengst))];
            end;
            if hans.efilestype(sel) == 2,
                fact = find(ismember(LUE.postnames(:,1),hans.efiles{sel}(13:end)));
            end;
            if (hans.efilestype(sel)==2)
                switch 1,
                    case ~isequal(LUE.btall,1),
                        if (hans.ausp{i,2} == LUE.postnames{fact,4}) & ((LUE.postnames{fact,3}==hans.ausp{i,3}) | (LUE.postnames{fact,2}==hans.ausp{i,3})), LUE.saveposthoc = [LUE.saveposthoc, sprintf('*')]; end;
                    otherwise,
                        if ((LUE.postnames{fact,3}==hans.ausp{i,2}) | (LUE.postnames{fact,2}==hans.ausp{i,2})), LUE.saveposthoc = [LUE.saveposthoc, sprintf('*')]; end;
                end;
                
                %fprintf('*');
            end;
            LUE.saveposthoc = [LUE.saveposthoc, sprintf('\n')];
        end;
        LUE.saveposthoc = [LUE.saveposthoc, sprintf('%s\n\n\n\n',repmat('-',1,(laengst+1)*(length(anames)+4)))];

%% Post-Hoc
        
        pval = nan(length(hans.werte),length(hans.werte));
        for i = 1:length(hans.werte) - 1,
            for j = i + 1 : length(hans.werte),
                if isequal(hans.ausp(i,1:size(adata,2)),hans.ausp(j,1:size(adata,2)))
                    % within
                    pt = nix_ldf1([hans.werteraw{i},hans.werteraw{j}]);
                    if isnan(eval(sprintf('pt.%s.pt1',LUE.statkind))),
                        pval(i,j) = -99;
                    else,
                        eval(sprintf('pval(%d,%d) = pt.%s.pt1;',i,j,LUE.statkind));
                    end;
                else,
                    % between
                    [~,~,pt] = nix_brunner_munzel(hans.werteraw{i},hans.werteraw{j});
                    if isnan(pt),
                        pval(j,i) = 0;
                    else
                        pval(j,i) = pt;
                    end;
                end;
            end;
        end;
        
        %LUE.saveposthoc = [LUE.saveposthoc, sprintf('-----------------------------------------\n')];
        LUE.saveposthoc = [LUE.saveposthoc, sprintf('------------- Posthoc Tests -------------\n\n\n')];
        %LUE.saveposthoc = [LUE.saveposthoc, sprintf('-----------------------------------------\n\n\n')];
        
        LUE.saveposthoc = [LUE.saveposthoc, eval(sprintf('sprintf(''%% %ds\\n'',''Within Test'')',15*(size(pval,1)+1)-6))];
        LUE.saveposthoc = [LUE.saveposthoc, sprintf('     %s\n',repmat('-',1,15*(size(pval,1)+1)-11))];
        for j = 1:size(pval,2),
            LUE.saveposthoc = [LUE.saveposthoc, sprintf('% 15d',j)];
        end;
        LUE.saveposthoc = [LUE.saveposthoc, sprintf('\n')];
        LUE.saveposthoc = [LUE.saveposthoc, sprintf('     %s\n',repmat('-',1,15*(size(pval,1)+1)-11))];
        for i = 1:size(pval,1),
            LUE.saveposthoc = [LUE.saveposthoc, sprintf('% 3d |',i)];
            for j = 1:size(pval,2),
                
                if i==j,
                    LUE.saveposthoc = [LUE.saveposthoc, sprintf('% 15s','-----------')];
                elseif ~isnan(pval(i,j)),
                    if pval(i,j) ~= -99,
                        aast = sprintf('%1.3f',pval(i,j));
                        if pval(i,j)< .05, aast=[aast,'*']; end;
                        if pval(i,j)< .01, aast=[aast,'*']; end;
                        if pval(i,j)<.001, aast=[aast,'*']; end;
                    else,
                        aast = sprintf('%d',NaN);
                    end;
                    LUE.saveposthoc = [LUE.saveposthoc, sprintf('% 15s',aast)];
                else,
                    LUE.saveposthoc = [LUE.saveposthoc, sprintf('% 15s',sprintf('%s','.'))];
                end;
            end;
            LUE.saveposthoc = [LUE.saveposthoc, sprintf('\n')];
        end;
        LUE.saveposthoc = [LUE.saveposthoc, sprintf('     %s\n',repmat('-',1,15*(size(pval,1)+1)-11))];
        LUE.saveposthoc = [LUE.saveposthoc, sprintf('     %s\n\n\n','Between Test')];
        
    end;
end;