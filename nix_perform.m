function nix_perform,

warning off

global LUE

addpath(fullfile(LUE.maindir,'ranktest'));
if LUE.type == 1, addpath(fullfile(LUE.maindir,'src')); end;

if ~exist(LUE.resultdir), mkdir(LUE.resultdir); end;
LUE.cluster = [];

sub_perform;

fprintf('\n\nResults stored in %s',LUE.resultdir);
fprintf('\n\n--------------- Finished ---------------\n\n');

fprintf('- Starting Result-Tool ... ');
nix_show_results(LUE.resultdir);


function sub_perform,
global LUE

%% Intro
LUE.niicontfile = 'for_result_tool.mat'; %'niicont.mat';

fprintf('\n\n- Scanning Image Data\n\n');

V = spm_vol(char(LUE.final(:,end)));
if LUE.type == 1, % check how many different kinds of lesions
    holder = [];
    for i = 1 : numel(V)
        holder = [holder; unique(V(i).private.dat(:,:,:))];
    end;
    holder = unique(holder);
    LUE.lesiontypes = numel(holder);
    LUE.lesiontypec = holder;
    LUE.typesc      = [1:numel(holder); holder'];
end;

fprintf('- Creating Vis. File\n\n');

if any(exist(fullfile(LUE.resultdir,LUE.niicontfile))), delete(fullfile(LUE.resultdir,LUE.niicontfile)); end;
w = matfile(fullfile(LUE.resultdir,LUE.niicontfile),'Writable',true);
if LUE.type == 1,
    
    %LUE.max     = 0;
    
    LUE.effects = fullfact(2*ones(1,length(LUE.amountfactors)+1))-1; 
    LUE.effects = LUE.effects(sum(LUE.effects')>1,:);
    %LUE.effects = LUE.effects(2:end,:);
    
    LUE.faccomb = fullfact([LUE.amountfactors,LUE.lesiontypes]);
    
    w.working   = zeros(size(LUE.faccomb,1),prod(V(1).dim));
    
elseif LUE.type == 2,
    w.working   = zeros(size(LUE.final,1),prod(V(1).dim));
end;

genvar  = []; novar = []; nogroupvar = []; fdrvar = [];
ugroups = unique(LUE.groups,'rows'); 
mat     = V(1).dim(1) * V(1).dim(2);

if LUE.type == 1,
    
    holder = 0; for i = 2:size(LUE.faccomb,2), holder = holder + factorial(size(LUE.faccomb,2)) / (factorial(i) * factorial(size(LUE.faccomb,2)-i)); end;
    LUE.fa = holder;
    
    result.lrt = zeros(size(LUE.effects,1),prod(V(1).dim));
    result.df  = zeros(size(LUE.effects,1),prod(V(1).dim));
    result.p   = ones(size(LUE.effects,1), prod(V(1).dim));
    
    novar      = zeros(prod(V(1).dim),1);
    
    %Post-Hoc anlegen
    if LUE.praepost == 1,
        %pmytable     = zeros(2,2,2);
        postcomb     = size(LUE.effects(:,1:end-1),2);
        LUE.postcomb = ff2n(postcomb);
        
        for i = size(LUE.postcomb,1) : -1 : 1, % alle außer 2 rausschmeißen
            if sum(LUE.postcomb(i,:)) ~= 2,
                LUE.postcomb(i,:) = [];
            end;
        end;
        
        if length(LUE.fac) > 2, %anlegen für alle Kombinationen, die nicht aktuell getestet werden, aber kontrolliert werden
            
            postcombside = cell(size(LUE.postcomb,1),1);
            
            for i = 1 : size(LUE.postcomb,1),
                for j = 1:size(LUE.postcomb,2),
                    if LUE.postcomb(i,j) == 0,
                        postcombside{i} = [postcombside{i}; LUE.fac(j)];
                    end;
                end;
            end;
            
            for i = 1 : length(postcombside),
                postcombside{i} = fullfact(postcombside{i});
%                 holder9 = postcombside{i};
%                 postcombside{i} = [1:holder9(1)]';
%                 for j = 2 : length(holder9),
%                     postcombside{i} = [repmat(postcombside{i}, holder9(j),1), repmat([1:holder9(j)]',size(postcombside{i},1),1)];
%                 end;
            end;
            
        else,
            postcombside = [];
            for i = 1:size(LUE.postcomb,1),
                postcombside{i} = 1;
            end;
        end;
        
        for i = 1 : size(LUE.postcomb,1),
            %if sum(LUE.postcomb(i,:)) == 2,
                ap   = find(LUE.postcomb(i,:));

                for m1 = 1 : LUE.fac(ap(1))-1,
                    for m2 = m1 + 1 : LUE.fac(ap(1)),
                        for n1 = 1 : LUE.fac(ap(2))-1,
                            for n2 = n1 + 1 : LUE.fac(ap(2)),
                                for k1 = 1 : size(postcombside{i},1),
                                    %for k2 = 1 : size(postcombside{i},2)
%                                         eval(sprintf('result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.lrt = zeros(size(LUE.effects,1),prod(V(1).dim));',i,m1,m2,n1,n2,k1));
%                                         eval(sprintf('result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.df  = zeros(size(LUE.effects,1),prod(V(1).dim));',i,m1,m2,n1,n2,k1));
%                                         eval(sprintf('result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.p   = zeros(size(LUE.effects,1),prod(V(1).dim));',i,m1,m2,n1,n2,k1));
                                        eval(sprintf('result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.lrt = zeros(1,prod(V(1).dim));',i,m1,m2,n1,n2,k1));
                                        eval(sprintf('result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.df  = zeros(1,prod(V(1).dim));',i,m1,m2,n1,n2,k1));
                                        eval(sprintf('result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.p   = ones( 1,prod(V(1).dim));',i,m1,m2,n1,n2,k1));
                                    %end;
                                end;
                            end;
                        end;
                    end;
                end;
            %end;
        end;
    end;
    
elseif LUE.type == 2,
    
    LUE.ranks        = reshape(nix_rank(reshape(LUE.within,1,prod(size(LUE.within)))),size(LUE.within,1),size(LUE.within,2));
    
    %% emtpy result files anlegen
    result.Img.F     = zeros(1,prod(V(1).dim));
    result.Img.p     = ones(1,prod(V(1).dim));
    result.wt1.F     = zeros(1,prod(V(1).dim));
    result.wt1.p     = ones(1,prod(V(1).dim));
    result.ImgXwt1.F = zeros(1,prod(V(1).dim));
    result.ImgXwt1.p = ones(1,prod(V(1).dim));
    
    if ~isequal(LUE.bt,1),
        result.bt1.F         = zeros(1,prod(V(1).dim));
        result.bt1.p         = ones(1,prod(V(1).dim));
        result.ImgXbt1.F     = zeros(1,prod(V(1).dim));
        result.ImgXbt1.p     = ones(1,prod(V(1).dim));
        result.bt1Xwt1.F     = zeros(1,prod(V(1).dim));
        result.bt1Xwt1.p     = ones(1,prod(V(1).dim));
        result.ImgXbt1Xwt1.F = zeros(1,prod(V(1).dim));
        result.ImgXbt1Xwt1.p = ones(1,prod(V(1).dim));
    end;
    if length(LUE.wt)>1,
        result.wt2.F         = zeros(1,prod(V(1).dim));
        result.wt2.p         = ones(1,prod(V(1).dim));
        result.ImgXwt2.F     = zeros(1,prod(V(1).dim));
        result.ImgXwt2.p     = ones(1,prod(V(1).dim));
        result.wt1Xwt2.F     = zeros(1,prod(V(1).dim));
        result.wt1Xwt2.p     = ones(1,prod(V(1).dim));
        result.ImgXwt1Xwt2.F = zeros(1,prod(V(1).dim));
        result.ImgXwt1Xwt2.p = ones(1,prod(V(1).dim));
    end;
    
    
    sammler     = 0;
    for i = 1:length(V),
        sammler = [sammler; unique(V(i).private.dat(:,:,:))];
    end;
    allgroups   = unique(sammler);
    
    if LUE.praepost,
        imgperm = [];
        for i = 1:length(allgroups) - 1,
            for j = i+1 : length(allgroups),
                imgperm = [imgperm; allgroups(i),allgroups(j)];
            end;
        end;
        
        for i = 1:LUE.btall,
            for y = 1:size(imgperm,1),
                for m = 1:size(LUE.within,2)-1,
                    for n = m + 1 : size(LUE.within,2),
                        eval(sprintf('extresult.G%02d.G%02d.W%02d.W%02d.F = zeros(1,prod(V(1).dim));',i,y,m,n));
                        eval(sprintf('extresult.G%02d.G%02d.W%02d.W%02d.p = ones(1,prod(V(1).dim));',i,y,m,n));
                    end;
                end;
            end;
        end;
        
        LUE.imgperm = imgperm;
    else,
        extresult = [];
    end;
end;

%% Reading & Performing
fprintf('- Reading Data and Performing Analyses:\n'); 

rfilter = Inf(prod(V(1).dim),1);

if LUE.core > 1, try, matlabpool(LUE.core); fprintf('\n\n'); end; end;

fprintf('[%s]',repmat(' ',1,V(1).dim(3)));
if LUE.type == 1,
    
%    disp('Performing Analysis and saving correct data');
    
    
    for slice = 1 : V(1).dim(3), 
        working = zeros(size(LUE.final,1), mat);
        for i = 1:size(LUE.final,1),
            working(i,:) = reshape( V(i).private.dat(:,:,slice), 1, mat);
        end;
        
        tworking = zeros(size(LUE.faccomb,1),size(working,2));
        
        novarh   = zeros(1,prod(size(V(1).dim(1:2))));
        
        for i = 1 : size(tworking,2) % parfor
            if length(unique(working(:,i))) == 1,
                %disp('NOVAR NOCH EINFÜGEN');
                %novar                                = [novar; [(slice - 1) * mat + 1 : (slice - 1) * mat +  prod(V(1).dim(1:2)) ]'];
                %novar((slice - 1) * mat + i)         = 1;
                novarh(i)         = 1;
                
            else,
                for j = 1 : size(LUE.final,1),
                    % zuweisen in tabelle
                    holder3                          = ones(size(LUE.faccomb,1),1);
                    for k = 1 : size(LUE.final,2)-1,
                        holder3                      = holder3 .* (LUE.final{j,k}==LUE.faccomb(:,k));
                    end;
                    holder3                          = holder3 .* (find(working(j,i)==LUE.lesiontypec) == LUE.faccomb(:,end));
                    tworking(find(holder3),i)        = tworking(find(holder3),i) + 1;
                end;
                
                mytable                              = reshape(tworking(:,i),[LUE.fac,LUE.lesiontypes]);
                
                reslog                               = nix_contingency(mytable,LUE.effects);
                
                result.lrt(:,(slice - 1) * mat + i)  = [reslog.lrt];
                result.df(:,(slice - 1) * mat + i)   = [reslog.df];
                result.p(:,(slice - 1) * mat + i)    = [reslog.p];
                
                rweg                                 = min(mytable);
                while any(size(rweg)>1),      rweg   = min(rweg); end;
                rfilter((slice - 1) * mat + i)       = rweg;
                
                
                %Post-hoc
                if LUE.praepost == 1,
                    for j = 1 : size(LUE.postcomb,1),
                        ap   = find(LUE.postcomb(j,:));
                        
                        for m1 = 1 : LUE.fac(ap(1))-1,
                            for m2 = m1 + 1 : LUE.fac(ap(1)),
                                for n1 = 1 : LUE.fac(ap(2))-1,
                                    for n2 = n1 + 1 : LUE.fac(ap(2)),
                                        for k1 = 1 : size(postcombside{j},1),
                                            
                                            astringq = 'pmytable{i} = reshape(mytable(';

                                            runnere  = 1;
                                            runneri  = 1;
                                            for v = 1:length(LUE.fac),
                                                if any(ismember(ap,v)) & (runneri==1),
                                                    astringq = sprintf('%s[%d,%d],',astringq,m1,m2);
                                                    runneri  = runneri + 1;
                                                elseif any(ismember(ap,v)),
                                                    astringq = sprintf('%s[%d,%d],',astringq,n1,n2);
                                                else,
                                                    astringq = sprintf('%s%d,',astringq,postcombside{j}(k1,runnere));
                                                    runnere  = runnere + 1;
                                                end;
                                            end;
                                            
                                            astringq = [astringq,'[1,2]),2,2,2);'];
                                            eval(astringq);
                                            a(i)     = nix_contingency(pmytable{i},[1,1,1]);
                                            
                                            %disp('ACHTUNG HIER EINS WÄHLEN');
                                            %s(i)     = sign(((pmytable{i}(1,1,1)/pmytable{i}(2,1,1))/(pmytable{i}(1,2,1)/pmytable{i}(2,2,2))) / ((pmytable{i}(1,1,2)/pmytable{i}(2,1,2))/(pmytable{i}(1,2,2)/pmytable{i}(2,2,2))));
                                            s(i)     = sign(((pmytable{i}(1,1,1)/(pmytable{i}(2,1,1)+pmytable{i}(1,1,1)))/(pmytable{i}(1,2,1)/(pmytable{i}(1,2,1)+pmytable{i}(2,2,2)))) / ((pmytable{i}(1,1,2)/(pmytable{i}(2,1,2)+pmytable{i}(2,1,2)))/(pmytable{i}(1,2,2)/(pmytable{i}(2,2,2)+pmytable{i}(2,2,2)))));
                                            
                                            if any(isnan(s(i))), 
                                                s(i) = 0; 
                                            end; % !!! evtl. problematisch
                                            
                                            eval(sprintf('result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.lrt((slice - 1) * mat + i) = a(%d).lrt * s(i);',j,m1,m2,n1,n2,k1,i));
                                            eval(sprintf('result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.df((slice - 1) * mat + i)  = a(%d).df;', j,m1,m2,n1,n2,k1,i));
                                            eval(sprintf('result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.p((slice - 1) * mat + i)   = a(%d).p;',  j,m1,m2,n1,n2,k1,i));
                                        end;
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;

            end;
        end;
        
        %LUE.max = max([LUE.max, max(max(max(tworking)))]);
        
        w.working(:,(slice - 1) * mat + 1 : slice * mat) = tworking;
        
        novar((slice - 1) * mat + 1 : slice * mat)       = novarh;
        
        pause(.001);
        eval(sprintf('fprintf(''%s'');',repmat('\b',1,V(1).dim(3)+1)));
        fprintf('%s%s]',repmat('=',1,slice),repmat(' ',1,V(1).dim(3)-slice));
        pause(.001);
    end;
    
    
    
elseif LUE.type == 2,
    for slice = 1 : V(1).dim(3),
        
        working = zeros(size(LUE.final,1), mat);
        for i = 1:size(LUE.final,1),
            working(i,:) = reshape( V(i).private.dat(:,:,slice), 1, mat);
        end;
        w.working(:,(slice - 1) * mat + 1 : slice * mat) = working;
        
        if length(unique(working)) < 2, % checks slice for variance
            novar = [novar; [(slice - 1) * mat + 1 : (slice - 1) * mat +  prod(V(1).dim(1:2)) ]'];
        else,
            res = []; asign = []; idx = [];
            tmpfilter = Inf(size(working,2),1);
            parfor i = 1:size(working,2), % schleife, damit nur Voxel mit Varianz aufgenommen werden und welche, bei denen mindestens 2 Personen zu jeder Voxelgruppe gehören
                ausp = unique(working(:,i));
                if length(ausp) == length(allgroups),
                    r = 1; rfdr = 1;
                    for j = 1:length(ausp),
                        agroup = find(working(:,i) == ausp(j));
                        for k = 1:size(ugroups,1),
                            for o = 1:size(ugroups,2),
%                                 if (slice - 1) * mat + i == 870690,
%                                     1,
%                                 end;
                                if length(find(LUE.groups(agroup,o)==ugroups(k,o))) < tmpfilter(i)%rfilter((slice - 1) * mat + i),
                                    tmpfilter(i) = length(find(LUE.groups(agroup,o)==ugroups(k,o)));
                                end;
                                if length(find(LUE.groups(agroup,o)==ugroups(k,o))) < LUE.fcorr, %*wtall
                                    rfdr = 0;
                                end;
                                if length(find(LUE.groups(agroup,o)==ugroups(k,o))) < 2, %*wtall
                                    r = 0;
                                    %break;
                                end;
                            end;
                        end;
                    end;
                    if r == 1,
                        if rfdr == 1,
                            fdrvar = [fdrvar; (slice - 1) * mat + i];
                        end;
                        
                        genvar = [genvar;(slice - 1) * mat + i];
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        switch 1,
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            case isequal(LUE.bt, 1) & (length(LUE.wt) == 1),
                                res(i).res = nix_f1ldf2([working(:,i), LUE.within], LUE.wt, 1);
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            case isequal(LUE.bt, 1) & (length(LUE.wt) == 2),
                                res(i).res = nix_f1ldf2([working(:,i), LUE.within], LUE.wt(1), LUE.wt(2));
                                
                                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            otherwise,
                                res(i).res = nix_f2ldf1([working(:,i), LUE.groups, LUE.within]);
                        end;
                        
                        idx(i,:) = [(slice - 1) * mat + i, i];
                        
                        if LUE.praepost,
                            for k = 1:LUE.btall,
                                for y = 1:size(imgperm,1),
                                    auswahl  = (LUE.groups == k);
                                    auswahl2 = sign((working(:,i) == imgperm(y,1)) + (working(:,i) == imgperm(y,2)));
                                    auswahl  = find(auswahl .* auswahl2);
                                    gr1      = find(working(auswahl,i) == imgperm(y,1));%allgroups(1));
                                    gr2      = find(working(auswahl,i) == imgperm(y,2));%allgroups(2));
                                    for m = 1:size(LUE.within,2)-1,
                                        for n = m + 1 : size(LUE.within,2),
                                            %res(i).a(k).a(m).a(n).res = nix_f1ldf2([working(auswahl,i), LUE.within(auswahl,[m,n])],2,1);
                                            %asign(i).a(k).a(m).a(n).asign = sign((median(LUE.ranks(auswahl(gr1),m)) - median(LUE.ranks(auswahl(gr2),m))) - (median(LUE.ranks(auswahl(gr1),n)) - median(LUE.ranks(auswahl(gr2),n))));
                                            res(i).a(k).a(y).a(m).a(n).res     = nix_f1ldf2([working(auswahl,i), LUE.within(auswahl,[m,n])],2,1);
                                            asign(i).a(k).a(y).a(m).a(n).asign = sign((median(LUE.ranks(auswahl(gr1),m)) - median(LUE.ranks(auswahl(gr2),m))) - (median(LUE.ranks(auswahl(gr1),n)) - median(LUE.ranks(auswahl(gr2),n))));
                                        end;
                                    end;
                                end;
                            end;
                        end;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    else,
                        nogroupvar = [nogroupvar; (slice - 1) * mat + i];
                    end;
                    
                else,
                    novar = [novar; (slice - 1) * mat + i];
                end;
            end;
            
            rfilter((slice - 1) * mat + 1:slice * mat) = tmpfilter;
            
            if length(find(idx)) > 0,
                aidx = find(idx(:,1));
                for i = 1 : length(aidx),
                    
                    eval(sprintf('result.wt1.F(%d) = res(%d).res.%s.Ft1;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                    eval(sprintf('result.wt1.p(%d) = res(%d).res.%s.pt1;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                    eval(sprintf('result.Img.F(%d) = res(%d).res.%s.Fgr1;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                    eval(sprintf('result.Img.p(%d) = res(%d).res.%s.pgr1;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                    eval(sprintf('result.ImgXwt1.F(%d) = res(%d).res.%s.Fgr1t1;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                    eval(sprintf('result.ImgXwt1.p(%d) = res(%d).res.%s.pgr1t1;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                    
                    switch 1
                        case isequal(LUE.bt, 1)  & (length(LUE.wt) == 2),
                            
                            eval(sprintf('result.wt2.F(%d) = res(%d).res.%s.Ft2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.wt2.p(%d) = res(%d).res.%s.pt2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.wt1Xwt2.F(%d) = res(%d).res.%s.Ft1t2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.wt1Xwt2.p(%d) = res(%d).res.%s.pt1t2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.ImgXwt2.F(%d) = res(%d).res.%s.Fgr1t2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.ImgXwt2.p(%d) = res(%d).res.%s.pgr1t2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.ImgXwt1Xwt2.F(%d) = res(%d).res.%s.Fgr1t1t2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.ImgXwt1Xwt2.p(%d) = res(%d).res.%s.pgr1t1t2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            
                        case ~isequal(LUE.bt, 1) & (length(LUE.wt) == 1),
                            
                            eval(sprintf('result.bt1.F(%d) = res(%d).res.%s.Fgr2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.bt1.p(%d) = res(%d).res.%s.pgr2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.ImgXbt1.F(%d) = res(%d).res.%s.Fgr1gr2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.ImgXbt1.p(%d) = res(%d).res.%s.pgr1gr2;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.bt1Xwt1.F(%d) = res(%d).res.%s.Fgr2t1;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.bt1Xwt1.p(%d) = res(%d).res.%s.pgr2t1;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.ImgXbt1Xwt1.F(%d) = res(%d).res.%s.Fgr1gr2t1;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            eval(sprintf('result.ImgXbt1Xwt1.p(%d) = res(%d).res.%s.pgr1gr2t1;',idx(aidx(i),1),idx(aidx(i),2),LUE.statkind));
                            
                    end;
                    
                    if LUE.praepost,
                        for k = 1:LUE.btall,
                            for y = 1:size(imgperm,1),
                                auswahl  = (LUE.groups == k);
                                auswahl2 = sign((working(:,aidx(i)) == imgperm(y,1)) + (working(:,aidx(i)) == imgperm(y,2)));
                                auswahl  = find(auswahl .* auswahl2);
                                for m = 1:size(LUE.within,2)-1,
                                    for n = m + 1 : size(LUE.within,2),
                                        %eval(sprintf('extresult.G%02d.W%02d.W%02d.F(%d) = res(%d).a(%d).a(%d).a(%d).res.%s.Fgr1t1 * asign(%d).a(%d).a(%d).a(%d).asign;', k, m, n, idx(aidx(i),1), idx(aidx(i),2), k, m, n, LUE.statkind, idx(aidx(i),2), k, m, n));
                                        %eval(sprintf('extresult.G%02d.W%02d.W%02d.p(%d) = res(%d).a(%d).a(%d).a(%d).res.%s.pgr1t1;', k, m, n, idx(aidx(i),1), idx(aidx(i),2), k, m, n, LUE.statkind));
                                        eval(sprintf('extresult.G%02d.G%02d.W%02d.W%02d.F(%d) = res(%d).a(%d).a(%d).a(%d).a(%d).res.%s.Fgr1t1 * asign(%d).a(%d).a(%d).a(%d).a(%d).asign;', k, y, m, n, idx(aidx(i),1), idx(aidx(i),2), k, y, m, n, LUE.statkind, idx(aidx(i),2), k, y, m, n));
                                        eval(sprintf('extresult.G%02d.G%02d.W%02d.W%02d.p(%d) = res(%d).a(%d).a(%d).a(%d).a(%d).res.%s.pgr1t1;', k, y, m, n, idx(aidx(i),1), idx(aidx(i),2), k, y, m, n, LUE.statkind));
                                    end;
                                end;
                            end;
                        end;
                    end;
                    
                end;
            end;
            
        end;
        
        pause(.001);
        eval(sprintf('fprintf(''%s'');',repmat('\b',1,V(1).dim(3)+1)));
        fprintf('%s%s]',repmat('=',1,slice),repmat(' ',1,V(1).dim(3)-slice));
        pause(.001);
    end;
end;
fprintf('\n');

if LUE.core > 1, try, matlabpool('close'); fprintf('\n'); end; end;

%% Write Results into Files

fprintf('\n- Writing Results into nii-File and calculating FDR adjustment\n'); 

V = V(1);

% Group Size
V.dt(1) = 8;
V.fname = fullfile(LUE.resultdir,'MinGroupSize.nii');
rfilter(find(isinf(rfilter))) = 0;
spm_write_vol(V,reshape(rfilter,V(1).dim));

V.dt(1) = 2; % für novar gebraucht

if LUE.type == 1,
    
    if isempty(fdrvar), fprintf('No voxel achieved the threshold of %d\n\n',LUE.fcorr); end;
    
    % No Lesion Variance
    novar2  = zeros(V(1).dim); novar2(find(novar)) = 1;
    V.fname = fullfile(LUE.resultdir,'NoVariance.nii');
    spm_write_vol(V,novar2);
    
    % FDR adjustments
    working   = w.working;
    [c,ia,ic] = unique(working(:,fdrvar)','rows');
    
    % Effects
    V.dt(1) = 16;
    for i = 1 : size(LUE.effects,1),
        fname     = sprintf('%s',sprintf('%s_X_',LUE.namesfinal{find(LUE.effects(i,1:end-1))})); fname = fname(1:end-3);
        if LUE.effects(i,end) == 1, fname = sprintf('%s_X_ImgData',fname); end;
        
        V.fname = fullfile(LUE.resultdir,['F_',fname,'.nii']);
        spm_write_vol(V, reshape(result.lrt(i,:),V(1).dim));
        
        V.fname = fullfile(LUE.resultdir,['df_',fname,'.nii']);
        spm_write_vol(V, reshape(result.df(i,:),V(1).dim));fname
        
        V.fname = fullfile(LUE.resultdir,['p_',fname,'.nii']);
        spm_write_vol(V, reshape(result.p(i,:),V(1).dim));
        
        V.fname = fullfile(LUE.resultdir,['padj_',fname,'.nii']);
        P = ones(size(working,2),1); if ~isempty(fdrvar), [~, ~, fdradj] = nix_fdr(result.p(i,fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
        spm_write_vol(V, reshape(P,V(1).dim));
    end;
    
    if LUE.praepost == 1,
        for j = 1 : size(LUE.postcomb,1),
            ap   = find(LUE.postcomb(j,:));
            for m1 = 1 : LUE.fac(ap(1))-1,
                for m2 = m1 + 1 : LUE.fac(ap(1)),
                    for n1 = 1 : LUE.fac(ap(2))-1,
                        for n2 = n1 + 1 : LUE.fac(ap(2)),
                            for k1 = 1 : size(postcombside{j},1),
                                
                                V.dt(1) = 16;
                                
                                %Naming Files
                                fname   = sprintf('%s',LUE.namesfinal{ap(1)});
                                if LUE.fac(ap(1))>2, fname = sprintf('%s.%s.%s',fname,LUE.levelnames{ap(1)}{m1},LUE.levelnames{ap(1)}{m2}); end;
                                fname   = sprintf('%s_X_%s',fname,LUE.namesfinal{ap(2)});
                                if LUE.fac(ap(2))>2, fname = sprintf('%s.%s.%s',fname,LUE.levelnames{ap(2)}{n1},LUE.levelnames{ap(2)}{n2}); end;
                                fname   = sprintf('%s_X_ImgData',fname);
                                if length(LUE.fac)>2,  
                                    fname = sprintf('%s_|',fname); 
                                    for v = 1:size(postcombside{j},2),
                                        if ~any(ismember(ap,postcombside{j}(k1,v))),
                                            fname = sprintf('%s_%s.%s',fname,LUE.namesfinal{v},LUE.levelnames{v}{postcombside{j}(k1,v)}); 
                                        end;
                                    end;
                                end;
                                
                                %Saving
                                V.fname = fullfile(LUE.resultdir,['Diff_Posthoc_for_',fname,'.nii']);
                                eval(sprintf('P = reshape(result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.lrt,V(1).dim);',j,m1,m2,n1,n2,k1));
                                spm_write_vol(V,P);
                                
                                V.fname = fullfile(LUE.resultdir,['p_Posthoc_for_',fname,'.nii']);
                                eval(sprintf('P = reshape(result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.p,V(1).dim);',  j,m1,m2,n1,n2,k1));
                                spm_write_vol(V,P);
                                
                                V.fname = fullfile(LUE.resultdir,['padj_Posthoc_for_',fname,'.nii']);
                                eval(sprintf('P = reshape(result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.p,V(1).dim);',  j,m1,m2,n1,n2,k1));
                                if ~isempty(fdrvar), [~, ~, fdradj] = nix_fdr(P(fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
                                spm_write_vol(V, reshape(P,V(1).dim));
                                
                                V.dt(1) = 4;
                                V.fname = fullfile(LUE.resultdir,['df_Posthoc_for_',fname,'.nii']);
                                eval(sprintf('P = reshape(result.post%02d.G%02d.G%02d.G%02d.G%02d.G%02d.df,V(1).dim);', j,m1,m2,n1,n2,k1));
                                spm_write_vol(V,P);
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
    
elseif LUE.type == 2,
    genvar = sort(genvar); novar = sort(novar); nogroupvar = sort(nogroupvar); fdrvar = sort(fdrvar);
    w.genvar     = genvar;
    w.novar      = novar;
    w.nogroupvar = nogroupvar;
    w.fdrvar     = fdrvar;
    
    V.fname = fullfile(LUE.resultdir,'NoVariance.nii');
    P = zeros(V.dim); P(novar) = 1; P(nogroupvar) = 2;
    spm_write_vol(V,P);
    
    % adapted FDR
    working   = w.working;
    %fdrvar    = w.fdrvar;
    [c,ia,ic] = unique(working(:,fdrvar)','rows');
    
    % V.fname = fullfile(LUE.resultdir,'NoGroupVar.nii');
    % P = zeros(V.dim); P(nogroupvar) = 1;
    % spm_write_vol(V,P);
    %
    % V.fname = fullfile(LUE.resultdir,'NotPerformable.nii');
    % P = zeros(V.dim); P(nogroupvar) = 1; P(novar) = 1;
    % spm_write_vol(V,P);
    
    % Within Main Effect
    V.dt(1) = 16;
    
    V.fname = fullfile(LUE.resultdir,['F_',LUE.wtnames{1},'.nii']);
    spm_write_vol(V,reshape(result.wt1.F,V.dim(1),V.dim(2),V.dim(3)));
    V.fname = fullfile(LUE.resultdir,['p_',LUE.wtnames{1},'.nii']);
    spm_write_vol(V,reshape(result.wt1.p,V.dim(1),V.dim(2),V.dim(3)));
    V.fname = fullfile(LUE.resultdir,['padj_',LUE.wtnames{1},'.nii']); P = ones(V.dim); 
%     if ~isempty(w.fdrvar), 
%         if LUE.fdrcheck == 1,
%             [~, ~, fdradj] = nix_fdr(result.wt1.p(w.fdrvar),.05,0); 
%             P(w.fdrvar) = fdradj; 
%         elseif LUE.fdrcheck==2,
%             %[~, ~, fdradj] = nix_fdr(result.wt1.p(fdrvar(c)),.05,0);
%             %P(fdrvar(ic)) = fdradj;
%         end;
%     end;
    if ~isempty(w.fdrvar), [~, ~, fdradj] = nix_fdr(result.wt1.p(fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
    spm_write_vol(V,reshape(P,V.dim(1),V.dim(2),V.dim(3)));
    
    % Lesion Main Effect
    V.fname = fullfile(LUE.resultdir,['F_ImgData.nii']);
    spm_write_vol(V,reshape(result.Img.F,V.dim(1),V.dim(2),V.dim(3)));
    V.fname = fullfile(LUE.resultdir,['p_ImgData.nii']);
    spm_write_vol(V,reshape(result.Img.p,V.dim(1),V.dim(2),V.dim(3)));
    V.fname = fullfile(LUE.resultdir,['padj_ImgData.nii']); P = ones(V.dim); %if ~isempty(w.fdrvar), [~, ~, fdradj] = nix_fdr(result.Img.p(w.fdrvar),.05,0);  P(w.fdrvar) = fdradj; end;
    if ~isempty(fdrvar), [~, ~, fdradj] = nix_fdr(result.Img.p(fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
    spm_write_vol(V,reshape(P,V.dim(1),V.dim(2),V.dim(3)));
    
    % Within X Lesion Interaction
    V.fname = fullfile(LUE.resultdir,['F_ImgData_X_',LUE.wtnames{1},'.nii']);
    spm_write_vol(V,reshape(result.ImgXwt1.F,V.dim(1),V.dim(2),V.dim(3)));
    V.fname = fullfile(LUE.resultdir,['p_ImgData_X_',LUE.wtnames{1},'.nii']);
    spm_write_vol(V,reshape(result.ImgXwt1.p,V.dim(1),V.dim(2),V.dim(3)));
    V.fname = fullfile(LUE.resultdir,['padj_ImgData_X_',LUE.wtnames{1},'.nii']); P = ones(V.dim); 
    if ~isempty(fdrvar), [~, ~, fdradj] = nix_fdr(result.ImgXwt1.p(fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
    spm_write_vol(V,reshape(P,V.dim(1),V.dim(2),V.dim(3)));
    
    
    
    switch 1
        case isequal(LUE.bt, 1)  & (length(LUE.wt) == 2),
            
            % 2. Within Main Effect
            V.fname = fullfile(LUE.resultdir,['F_',LUE.wtnames{2},'.nii']);
            spm_write_vol(V,reshape(result.wt2.F,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['p_',LUE.wtnames{2},'.nii']);
            spm_write_vol(V,reshape(result.wt2.p,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['padj_',LUE.wtnames{2},'.nii']); P = ones(V.dim); %if ~isempty(w.fdrvar), [~, ~, fdradj] = nix_fdr(result.wt2.p(w.fdrvar),.05,0);  P(w.fdrvar) = fdradj;end;
            if ~isempty(fdrvar), [~, ~, fdradj] = nix_fdr(result.wt2.p(fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
            spm_write_vol(V,reshape(P,V.dim(1),V.dim(2),V.dim(3)));
            
            % 2. Within X Lesion Interaction
            V.fname = fullfile(LUE.resultdir,['F_ImgData_X_',LUE.wtnames{2},'.nii']);
            spm_write_vol(V,reshape(result.ImgXwt2.F,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['p_ImgData_X_',LUE.wtnames{2},'.nii']);
            spm_write_vol(V,reshape(result.ImgXwt2.p,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['padj_ImgData_X_',LUE.wtnames{2},'.nii']); P = ones(V.dim); %if ~isempty(w.fdrvar), [~, ~, fdradj] = nix_fdr(result.ImgXwt2.p(w.fdrvar),.05,0); P(w.fdrvar) = fdradj;end;
            if ~isempty(fdrvar), [~, ~, fdradj] = nix_fdr(result.ImgXwt2.p(fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
            spm_write_vol(V,reshape(P,V.dim(1),V.dim(2),V.dim(3)));
            
            % 1. Within X 2. Within X Lesion Interaction
            V.fname = fullfile(LUE.resultdir,['F_ImgData_X_',LUE.wtnames{1},'X',LUE.wtnames{2},'.nii']);
            spm_write_vol(V,reshape(result.ImgXwt1Xwt2.F,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['p_ImgData_X_',LUE.wtnames{1},'X',LUE.wtnames{2},'.nii']);
            spm_write_vol(V,reshape(result.ImgXwt1Xwt2.p,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['padj_ImgData_X_',LUE.wtnames{1},'X',LUE.wtnames{2},'.nii']); P = ones(V.dim); %if ~isempty(w.fdrvar), [~, ~, fdradj] = nix_fdr(result.ImgXwt1Xwt2.p(w.fdrvar),.05,0); P(w.fdrvar) = fdradj; end;
            if ~isempty(fdrvar), [~, ~, fdradj] = nix_fdr(result.ImgXwt1Xwt2.p(fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
            spm_write_vol(V,reshape(P,V.dim(1),V.dim(2),V.dim(3)));
            
        case ~isequal(LUE.bt, 1) & (length(LUE.wt) == 1),
            
            % Between Main Effect
            V.fname = fullfile(LUE.resultdir,['F_',LUE.btnames{1},'.nii']);
            spm_write_vol(V,reshape(result.bt1.F,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['p_',LUE.btnames{1},'.nii']);
            spm_write_vol(V,reshape(result.bt1.p,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['padj_',LUE.btnames{1},'.nii']);  P = ones(V.dim); %if ~isempty(w.fdrvar), [~, ~, fdradj] = nix_fdr(result.bt1.p(w.fdrvar),.05,0);P(w.fdrvar) = fdradj; end;
            if ~isempty(fdrvar), [~, ~, fdradj] = nix_fdr(result.bt1.p(fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
            spm_write_vol(V,reshape(P,V.dim(1),V.dim(2),V.dim(3)));
            
            % Between X Lesion Interacton
            V.fname = fullfile(LUE.resultdir,['F_ImgData_X_',LUE.btnames{1},'.nii']);
            spm_write_vol(V,reshape(result.ImgXbt1.F,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['p_ImgData_X_',LUE.btnames{1},'.nii']);
            spm_write_vol(V,reshape(result.ImgXbt1.p,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['padj_ImgData_X_',LUE.btnames{1},'.nii']); P = ones(V.dim); %if ~isempty(w.fdrvar), [~, ~, fdradj] = nix_fdr(result.ImgXbt1.p(w.fdrvar),.05,0);  P(w.fdrvar) = fdradj; end;
            if ~isempty(fdrvar), [~, ~, fdradj] = nix_fdr(result.ImgXbt1.p(fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
            spm_write_vol(V,reshape(P,V.dim(1),V.dim(2),V.dim(3)));
            
            % Between X Within Interacton
            V.fname = fullfile(LUE.resultdir,[sprintf('F_%s_X_%s',LUE.btnames{1},LUE.wtnames{1}),'.nii']);
            spm_write_vol(V,reshape(result.bt1Xwt1.F,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,[sprintf('p_%s_X_%s',LUE.btnames{1},LUE.wtnames{1}),'.nii']);
            spm_write_vol(V,reshape(result.bt1Xwt1.p,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,[sprintf('padj_%s_X_%s',LUE.btnames{1},LUE.wtnames{1}),'.nii']); P = ones(V.dim); %if ~isempty(w.fdrvar), [~, ~, fdradj] = nix_fdr(result.bt1Xwt1.p(w.fdrvar),.05,0);  P(w.fdrvar) = fdradj; end;
            if ~isempty(fdrvar), [~, ~, fdradj] = nix_fdr(result.bt1Xwt1.p(fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
            spm_write_vol(V,reshape(P,V.dim(1),V.dim(2),V.dim(3)));
            
            % Within X Between X Lesion Interaction
            V.fname = fullfile(LUE.resultdir,['F_ImgData_X_',LUE.btnames{1},'X',LUE.wtnames{1},'.nii']);
            spm_write_vol(V,reshape(result.ImgXbt1Xwt1.F,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['p_ImgData_X_',LUE.btnames{1},'X',LUE.wtnames{1},'.nii']);
            spm_write_vol(V,reshape(result.ImgXbt1Xwt1.p,V.dim(1),V.dim(2),V.dim(3)));
            V.fname = fullfile(LUE.resultdir,['padj_ImgData_X_',LUE.btnames{1},'X',LUE.wtnames{1},'.nii']);  P = ones(V.dim); %if ~isempty(w.fdrvar), [~, ~, fdradj] = nix_fdr(result.ImgXbt1Xwt1.p(w.fdrvar),.05,0); P(w.fdrvar) = fdradj; end;
            if ~isempty(fdrvar), [~, ~, fdradj] = nix_fdr(result.ImgXbt1Xwt1.p(fdrvar(ia)),.05,1); P(fdrvar) = fdradj(ic); end;
            spm_write_vol(V,reshape(P,V.dim(1),V.dim(2),V.dim(3)));
            
    end;
    
    
    
    
    if LUE.praepost,
        LUE.postnames = cell(0,0);
        for k = 1:LUE.btall,
            auswahl  = (LUE.groups == k);
            for y = 1:size(imgperm,1),
                for m = 1:size(LUE.within,2)-1,
                    for n = m + 1 : size(LUE.within,2),
                        % FileName
                        if ~isequal(LUE.btall,1), fname = [LUE.btnames{1},LUE.levelnames{1}{k},'_']; else, fname = []; end;%num2str(k),'_']; else, fname = []; end;
                        if size(imgperm,1) > 1, fname = [fname, 'Img',num2str(round(imgperm(y,1))),'_Img',num2str(round(imgperm(y,2))),'_']; end;
                        fname = [fname, LUE.wtnamesxls{m}, '_', LUE.wtnamesxls{n}];
                        LUE.postnames(end+1,:) = {fname,m,n,median(LUE.groups(auswahl)),imgperm(y,1),imgperm(y,2)};
                        
                        % Create Files
                        V.fname = fullfile(LUE.resultdir,['Diff_Posthoc_for_',fname,'.nii']);
                        eval(sprintf('spm_write_vol(V,reshape(extresult.G%02d.G%02d.W%02d.W%02d.F,%d,%d,%d));',k,y,m,n,V.dim(1),V.dim(2),V.dim(3)));
                        V.fname = fullfile(LUE.resultdir,['p_Posthoc_for_',fname,'.nii']);
                        eval(sprintf('spm_write_vol(V,reshape(extresult.G%02d.G%02d.W%02d.W%02d.p,%d,%d,%d));',k,y,m,n,V.dim(1),V.dim(2),V.dim(3)));
                        
                           % FDR Korrektur
                        [c,ia,ic] = unique(working(find(auswahl),fdrvar)','rows');
                        V.fname = fullfile(LUE.resultdir,['padj_Posthoc_for_',fname,'.nii']);                        
                        %P = ones(V.dim); if ~isempty(w.fdrvar), eval(sprintf('[~, ~, fdradj] = nix_fdr(extresult.G%02d.G%02d.W%02d.W%02d.p(w.fdrvar),.05,0);',k,y,m,n)); P(w.fdrvar) = fdradj; end;
                        P = ones(V.dim); if ~isempty(fdrvar), eval(sprintf('[~, ~, fdradj] = nix_fdr(extresult.G%02d.G%02d.W%02d.W%02d.p(fdrvar(ia)),.05,1);',k,y,m,n)); P(fdrvar) = fdradj(ic); end;                                               
                        spm_write_vol(V,P);
                    end;
                end;
            end;
        end;
    end;
end;
disp('Zeile 276 noch wählen welches');
w.LUE = LUE;
%w.rfilter = rfilter;
