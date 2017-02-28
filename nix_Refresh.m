function nix_Refresh(varargin)

warning('off');
pause(.001); % so that everything from buttons and inputs can get in order

global hans LUE inn

clc
LUE.saveposthoc = [];
%inn       = matfile(fullfile(hans.effectpath,LUE.niicontfile));
%inn       = load(fullfile(hans.effectpath,LUE.niicontfile));

if ~any(ismember([6,18,26],str2num(get(hans.clustui(3),'String')))),
    check = 0;
    while check == 0,
        answer = inputdlg(sprintf('Please set a valid connectivity value.\nValid values are 6(surface), 18(edge) or 26(corner).\nFor further information see the manual or in Matlab: help spm_bwlabel.'),'Connectivity Value',1,{'26'});
        try, if any(ismember([6,18,26],str2num(answer{1}))), check = 1; end; end;
    end;
    set(hans.clustui(3),'String',answer{1});
end;

sel            = get(hans.effectpopup,'Value');
oldnfilter     = hans.nfilter;
hans.nfilter   = str2num(get(hans.clustui(4),'string'));
oldpuval       = hans.puval;
hans.puval     = str2num(get(hans.clustui(1),'String'));

%% Funkt Bild aufnehmen
if ~hans.alreadyloaded(sel,1)% ((nfilter==0)&(~hans.alreadyloaded(sel,1)))|(sel==1), % Standard
    
    afile = hans.efiles{sel};
    
    if sel == 1,
        hans.resultinf(sel,1)  = spm_vol(fullfile(hans.effectpath,[afile,'.nii']));
    else,
        if hans.efilestype(sel,1)==1,
            hans.resultinf(sel,1)  = spm_vol(fullfile(hans.effectpath,['F_',afile,'.nii']));
        else,
            hans.resultinf(sel,1)  = spm_vol(fullfile(hans.effectpath,['Diff_',afile,'.nii']));
        end;
        
        hans.plink(sel,1) = spm_vol(fullfile(hans.effectpath,['p_',afile,'.nii']));
        if LUE.type == 2,
            hans.padjlink(sel,1)   = spm_vol(fullfile(hans.effectpath,['padj_',afile,'.nii']));
        end;
        
        hans.alreadyloaded(sel,1)  = 1;
        anew                       = 1;
    end;
    
else,
    anew        = 1;
    
end;

% ClusterInfos befüllen
if sel > 1,
    if (~isequal(hans.clust(sel).conn,str2num(get(hans.clustui(3),'String')))) | (anew==1) | (hans.nfilter~=oldnfilter) | (oldpuval ~= hans.puval),
        hans.clust(sel).conn = str2num(get(hans.clustui(3),'String'));
        
        %[aclust,anum] = spm_bwlabel(hans.resultinf(sel).private.dat(:,:,:) .* (hans.plink(sel).private.dat(:,:,:)<str2num(get(hans.clustui(1),'String'))) .* (hans.rfilter >= hans.nfilter),hans.clust(sel).conn);
        %[aclust,anum]              = spm_bwlabel(hans.resultinf(sel).private.dat(:,:,:) .* (hans.rfilter >= hans.nfilter),hans.clust(sel).conn);
        [aclust,anum]              = spm_bwlabel(hans.resultinf(sel).private.dat(:,:,:) .* (hans.rfilter >= hans.nfilter) .* (hans.plink(sel).private.dat(:,:,:) <= str2num(get(hans.clustui(1),'String'))),hans.clust(sel).conn);
        %[aclust,anum] = spm_bwlabel(hans.plink(sel).private.dat(:,:,:) .* (hans.plink(sel).private.dat(:,:,:)<str2num(get(hans.clustui(1),'String'))),hans.clust(sel).conn);
        hans.clust(sel).clust      = aclust;
        %hans.clust(sel).out = zeros(anum,4);
        hans.clust(sel).out        = [];
        hans.clustparts(sel).parts = [];
        for i = 1 : anum,
            ach2                               = find(aclust==i);
            if length(ach2) >= 5,
                locmax                         = find(abs(hans.resultinf(sel).private.dat(ach2))==max(abs(hans.resultinf(sel).private.dat(ach2)))); locmax = locmax(1);
                hans.clust(sel).out(end+1,1:3) = hans.aXYZ(:,ach2(locmax));
                hans.clust(sel).out(end,4:5)   = [length(ach2), i];
                
                schnitt = sort(hans.atlas(ach2));
                weg     = (schnitt<0); schnitt(weg) = []; % necessary, because spm12 in windows produces strange imcalc results
                schnitt = [0; schnitt];
                amount  = hist(schnitt, schnitt(end)) / (length(schnitt)-1);
                [~,b]   = sort(amount,'descend');
                ma      = max(find(amount(b)>=.05));
                b       = b(1:ma);
                weg = find(b==1); b(weg) = [];
                hans.clustparts(sel).parts{end+1} = '';
                for bbb = 1 : length(b),
                    ac  = find([hans.atlasref{:,1}] == b(bbb));
                    hans.clustparts(sel).parts{end} = sprintf('%s | %s: %1.1f%%',hans.clustparts(sel).parts{end}, hans.atlasref{ac,2}, 100*amount(b(bbb)));
                end;
                hans.clustparts(sel).parts{end} = hans.clustparts(sel).parts{end}(4:end);
            end;
        end;
        if ~isempty(hans.clust(sel).out),
            [~,b] = sort(hans.clust(sel).out(:,4),'descend'); %sortieren
            %hans.clust(sel).out = [hans.clust(sel).out(b(end:-1:1),:),b(end:-1:1)]; %sortieren
            hans.clust(sel).out = [hans.clust(sel).out(b,:)];%,b];
            hans.clustparts(sel).parts = hans.clustparts(sel).parts(b);
        else,
            hans.clustparts(sel).parts = [];
        end;
        hans.last = -1;
    end;
else,
    hans.clust(sel).out = [];
    hans.last = -1;
end;

%% ClusterInfos befüllen
if ~isequal(hans.last,sel),
    hans.last = sel;
    astring = cell(size(hans.clust(sel).out,1),1);
    for i = 1:size(hans.clust(sel).out,1),
        astring{i} = sprintf('% 7.1f % 7.1f% 7.1f % 7.1d  %s',hans.clust(sel).out(i,1:4),hans.clustparts(sel).parts{i});
    end;
    astring{end+1} = 'None';
    set(hans.clusterpopup,'String',astring);
    if sel == 1, set(hans.clusterpopup,'Value',1); end;
end;

%% Strukt Bild aufnehmen und anpassen
if isempty(hans.underlay),
    hans.underlayinf = spm_vol(hans.underlayfile);
    [underlay, uXYZ] = spm_read_vols(hans.underlayinf);
    %hans.uXYZ = round(uXYZ*100)/100;
    %nochfolgend Zurechtschneiden, damit funktionell und strukturell
    %gleiche Ausschnitte
    if (max(uXYZ(1,:))*1.05<max(hans.aXYZ(1,:)))|(min(uXYZ(1,:))*1.05>min(hans.aXYZ(1,:))), error('The Functional Picture exceeds the Underlay Picture. This was not designated. Please choose another Underlay.'); end;
    if (max(uXYZ(2,:))*1.05<max(hans.aXYZ(2,:)))|(min(uXYZ(2,:))*1.05>min(hans.aXYZ(2,:))), error('The Functional Picture exceeds the Underlay Picture. This was not designated. Please choose another Underlay.'); end;
    if (max(uXYZ(3,:))*1.05<max(hans.aXYZ(3,:)))|(min(uXYZ(3,:))*1.05>min(hans.aXYZ(3,:))), error('The Functional Picture exceeds the Underlay Picture. This was not designated. Please choose another Underlay.'); end;
    
    raus = find(round(uXYZ(1,:))>max(round(hans.aXYZ(1,:)))); uXYZ(:,raus) = []; underlay(raus) = [];
    raus = find(round(uXYZ(1,:))<min(round(hans.aXYZ(1,:)))); uXYZ(:,raus) = []; underlay(raus) = [];
    hans.underlayinf.dim(1) = (max(uXYZ(1,:)) + abs(min(uXYZ(1,:)))) / abs(hans.underlayinf.mat(1,1)) + 1;
    
    raus = find(round(uXYZ(2,:))>max(round(hans.aXYZ(2,:)))); uXYZ(:,raus) = []; underlay(raus) = [];
    raus = find(round(uXYZ(2,:))<min(round(hans.aXYZ(2,:)))); uXYZ(:,raus) = []; underlay(raus) = [];
    hans.underlayinf.dim(2) = (max(uXYZ(2,:)) + abs(min(uXYZ(2,:)))) / abs(hans.underlayinf.mat(2,2)) + 1;
    
    raus = find(round(uXYZ(3,:))>max(round(hans.aXYZ(3,:)))); uXYZ(:,raus) = []; underlay(raus) = [];
    raus = find(round(uXYZ(3,:))<min(round(hans.aXYZ(3,:)))); uXYZ(:,raus) = []; underlay(raus) = [];
    hans.underlayinf.dim(3) = (max(uXYZ(3,:)) + abs(min(uXYZ(3,:)))) / abs(hans.underlayinf.mat(3,3)) + 1;
    
    %reshape
    abst = max(diff(sort(uXYZ')));
    hans.underlayinf.dim(1) = length(min(uXYZ(1,:)):abst(1):max(uXYZ(1,:)));
    hans.underlayinf.dim(2) = length(min(uXYZ(2,:)):abst(2):max(uXYZ(2,:)));
    hans.underlayinf.dim(3) = length(min(uXYZ(3,:)):abst(3):max(uXYZ(3,:)));
    underlay = reshape(underlay,hans.underlayinf.dim(1),hans.underlayinf.dim(2),hans.underlayinf.dim(3));
    
    hans.underlay = underlay;
    hans.uXYZ = uXYZ;
    
    acur = diff(uXYZ');
    for i = 1:3
        af              = find(acur(:,i)~=0);
        hans.uXYZdir(i) = mode(acur(af,i));
    end;
end;

%% atlas File einladen
if isempty(hans.atlas),
    V             = spm_vol(hans.atlasfile);
    hans.atlas    = round(spm_read_vols(V));
    
    hans.atlasref = [];
    fid           = fopen(fullfile(LUE.resultdir,'atlas.txt'),'r');
    string        = fgets(fid);
    while ~isequal(string,-1),
        try,
            weg = strfind(string,sprintf('\t'));
            hans.atlasref{end+1,1} = str2num(string(1:weg(1)-1));
            hans.atlasref{end,  2} = deblank(string(weg(1)+1:end));
            %             as    = textscan(string,'%d\t%s');
            %             if ~isempty(as{1}),
            %                 hans.atlasref{end+1,1} = as{1};
            %                 hans.atlasref{end,2}   = as{2}{1};
            %             end;
        end;
        string    = fgets(fid);
    end;
    fclose(fid);
end;

%% aktuelles x,y,z bestimmen
if length(varargin)<1|5,
    x = str2num(get(hans.edit(1),'String'));
    y = str2num(get(hans.edit(2),'String'));
    z = str2num(get(hans.edit(3),'String'));
end;

if length(varargin)>0,
    if varargin{1}==5,
        aentf = hans.clustcoord{sel} - repmat([x,y,z],size(hans.clustcoord{sel},1),1);
        aentf = sum(aentf'.^2)';
        aentf = find(aentf == min(aentf));
        set(hans.clusterpopup,'Value',aen>tf);
        varargin{1} = 4;
    end;
    
    if varargin{1}==1,
        a = get(gca,'currentPoint');
        x = str2num(get(hans.edit(1),'String'));
        y = a(1,1); ybs = [1,hans.resultinf(sel).dim(2)]; yrs = [min(hans.aXYZ(2,:)), max(hans.aXYZ(2,:))]; y = (y-ybs(1))/diff(ybs) * diff(yrs) + yrs(1);
        z = a(1,2); zbs = [1,hans.resultinf(sel).dim(3)]; zrs = [max(hans.aXYZ(3,:)),min(hans.aXYZ(3,:))]; z = (z-zbs(1))/diff(zbs) * diff(zrs) + zrs(1);
    elseif varargin{1}==2,
        a = get(gca,'currentPoint');
        x = a(1,1); xbs = [1,hans.resultinf(sel).dim(1)]; xrs = [min(hans.aXYZ(1,:)), max(hans.aXYZ(1,:))]; x = (x-xbs(1))/diff(xbs) * diff(xrs) + xrs(1);
        y = str2num(get(hans.edit(2),'String'));
        z = a(1,2); zbs = [1,hans.resultinf(sel).dim(3)]; zrs = [max(hans.aXYZ(3,:)), min(hans.aXYZ(3,:))]; z = (z-zbs(1))/diff(zbs) * diff(zrs) + zrs(1);
    elseif varargin{1}==3,
        a = get(gca,'currentPoint');
        x = a(1,1); xbs = [1,hans.resultinf(sel).dim(1)]; xrs = [min(hans.aXYZ(1,:)), max(hans.aXYZ(1,:))]; x = (x-xbs(1))/diff(xbs) * diff(xrs) + xrs(1);
        y = a(1,2); ybs = [1,hans.resultinf(sel).dim(2)]; yrs = [max(hans.aXYZ(2,:)), min(hans.aXYZ(2,:))]; y = (y-ybs(1))/diff(ybs) * diff(yrs) + yrs(1);
        z = str2num(get(hans.edit(3),'String'));
    elseif varargin{1}==4,
        selcluster = get(hans.clusterpopup,'Value');
        if ~(selcluster>size(hans.clust(sel).out,1)),
            x = hans.clust(sel).out(selcluster,1);
            y = hans.clust(sel).out(selcluster,2);
            z = hans.clust(sel).out(selcluster,3);
        end;
    end;
end;

%% x,y,z anpassen, vorbereitungen

% adists = hans.aXYZ-repmat([x;y;z],1,size(hans.aXYZ,2)); adists = sum(adists.^2);
% aact = find(adists==min(adists)); aact = aact(1);
% h = round(hans.aXYZ(:,aact)'*10)/10; x=h(1); y=h(2); z=h(3);
[h, aact] = nix_dist(x,y,z);
x=h(1); y=h(2); z=h(3);
LUE.aact  = aact;
LUE.coord = [x,y,z];

udists = hans.uXYZ-repmat([x;y;z],1,size(hans.uXYZ,2)); udists = sum(udists.^2);
uact = find(round(udists*100)/100==min(round(udists*100)/100)); uact = uact(1);

cross_farbe = [0,.2,1];
cont_color = [.2,.2,.2; .0,.6,.2];
ajet_abstufungen = 100;
if hans.efilestype(sel)==1,
    ajet       = jet(ajet_abstufungen*2);
    ajet       = ajet(size(ajet,1)/2+1:end,:);
    functrange = [0, str2num(get(hans.thres(2),'String'))];
else,
    ajet       = jet(ajet_abstufungen);
    functrange = [-1*str2num(get(hans.thres(2),'String')), str2num(get(hans.thres(2),'String'))];
end;
amaprange  = [0, str2num(get(hans.thres(1),'String'))];
amapint    = [0,str2num(get(hans.thres(3),'String'))];

ov_verkl      = max(hans.underlayinf.dim);
re_verkl      = max(hans.resultinf(sel).dim);

%% malen
for i = 1:3,
    % Strukt malen
    anot = [1:3]; anot(i) = [];
    auswahl = find(round(hans.uXYZ(i,uact)*100)/100 == (round(hans.uXYZ(i,:)*100)/100));
    ubild = rot90(reshape(hans.underlay(auswahl),hans.underlayinf.dim(anot(1)),hans.underlayinf.dim(anot(2))),1);
    %if i~=1, ubild = fliplr(ubild); end;
    %if hans.underlayinf.mat(i,i) < 0, ubild = fliplr(ubild); end;
    %if hans.uXYZdir(i) > 0, ubild = fliplr(ubild); end;
    if hans.aXYZdir(anot(1)) < 0, ubild = fliplr(ubild); end;
    if hans.aXYZdir(anot(2)) < 0, ubild = flipud(ubild); end;
    
    ubild = ubild/max(max(ubild)); ubild(:,:,2) = ubild(:,:,1); ubild(:,:,3) = ubild(:,:,1);
    ubild(find(ubild<0)) = 0; ubild(find(ubild>1)) = 1;
    ui = image(round(ubild*100)/100,'Parent',hans.und(i));
    axis(hans.und(i),'off'); %'equal');
    if hans.underlayinf.dim(anot(1))<ov_verkl, m = get(hans.und(i),'XLIM'); set(hans.und(i),'XLIM', m + (ov_verkl-diff(m))/2* [-1,1]); end; %Achsenanpassung
    if hans.underlayinf.dim(anot(2))<ov_verkl, m = get(hans.und(i),'YLIM'); set(hans.und(i),'YLIM', m + (ov_verkl-diff(m))/2* [-1,1]); end;
    
    %Funkt malen
    auswahl = find((hans.aXYZ(i,aact)==hans.aXYZ(i,:)));
    if isequal(hans.efiles{sel},'NoVariance'), % damit novariance immer angezeigt wird
        hans.nfilter = 0; 
    end; 
    
    %abild = rot90(reshape(hans.resultinf(sel).private.dat(auswahl),hans.resultinf(sel).dim(anot(1)),hans.resultinf(sel).dim(anot(2))),1);
    abild = rot90(reshape(hans.resultinf(sel).private.dat(auswahl) .* (hans.rfilter(auswahl)'>=hans.nfilter),hans.resultinf(sel).dim(anot(1)),hans.resultinf(sel).dim(anot(2))),1);
    
    %if i~=1, abild = fliplr(abild); end;
    %if hans.resultinf(sel).mat(i,i) > 0, abild = fliplr(abild); end;
    if hans.aXYZdir(anot(1)) < 0, abild = fliplr(abild); end;
    if hans.aXYZdir(anot(2)) < 0, abild = flipud(abild); end;
    
    % Alpha-Map bestimmen
    alpha_map = (abs(abild)-amaprange(1))/((amaprange(2)-amaprange(1)));
    alpha_map(find(alpha_map>1)) = 1; alpha_map(find(alpha_map<0)) = 0; % falls außerhalb der Grenzen
    alpha_map = (alpha_map + amapint(1)) * (amapint(2) - amapint(1));
    
    abild = (abild-functrange(1))/(functrange(2)-functrange(1));
    abild(find(abild>1)) = 1; abild(find(abild<.01)) = .01; % falls außerhalb der Grenzen
    abild = round(abild*ajet_abstufungen);
    
    %%% !!! fuer nans
    anan = find(isnan(abild)); abild(anan) = 1;
    abild = ajet(abild,:);
    abild(anan,:) = 0;
    
    abild = reshape(abild,hans.resultinf(sel).dim(anot(2)),hans.resultinf(sel).dim(anot(1)),3);
    ai = image(abild,'Parent',hans.brah(i));
    if hans.resultinf(sel).dim(anot(1))<re_verkl, m = get(hans.brah(i),'XLIM'); set(hans.brah(i),'XLIM', m + (re_verkl-diff(m))/2* [-1,1]); end; %Achsenanpassung
    if hans.resultinf(sel).dim(anot(2))<re_verkl, m = get(hans.brah(i),'YLIM'); set(hans.brah(i),'YLIM', m + (re_verkl-diff(m))/2* [-1,1]); end;
    set(ai,'AlphaData',alpha_map); %'ButtonDownFcn',sprintf('nix_Refresh(%d)',i)
    axis(hans.brah(i),'off','equal');
    
    if sel > 1,
        % Contouren für p und padj machen
        contp        = rot90(reshape(((hans.plink(sel).private.dat(auswahl)) + (hans.rfilter(auswahl)' < hans.nfilter)),   hans.resultinf(sel).dim(anot(1)),hans.resultinf(sel).dim(anot(2))),1);
        if LUE.type == 2,
            contpadj = rot90(reshape(hans.padjlink(sel).private.dat(auswahl) + (hans.rfilter(auswahl)' < hans.nfilter),hans.resultinf(sel).dim(anot(1)),hans.resultinf(sel).dim(anot(2))),1);
        end;
        
        if hans.aXYZdir(anot(1)) < 0, contp = fliplr(contp); if LUE.type == 2, contpadj = fliplr(contpadj); end; end;
        if hans.aXYZdir(anot(2)) < 0, contp = flipud(contp); if LUE.type == 2, contpadj = flipud(contpadj); end; end;
        hold(hans.brah(i),'on');
        contp    = contp    < str2num(get(hans.clustui(1),'String'));
        contour(contp   ,1,'Color',cont_color(1,:),'Parent',hans.brah(i));
        if LUE.type == 2,
            contpadj = contpadj < str2num(get(hans.clustui(2),'String'));
            contour(contpadj,1,'Color',cont_color(2,:),'Parent',hans.brah(i));
        end;
        hold(hans.brah(i),'off');
        if hans.resultinf(sel).dim(anot(1))<re_verkl, m = get(hans.brah(i),'XLIM'); set(hans.brah(i),'XLIM', m + (re_verkl-diff(m))/2* [-1,1]); end;
        if hans.resultinf(sel).dim(anot(2))<re_verkl, m = get(hans.brah(i),'YLIM'); set(hans.brah(i),'YLIM', m + (re_verkl-diff(m))/2* [-1,1]); end;
    end;
    
    % Cross malen
    %cross = NaN(size(abild,1),size(abild,2),3);
    cross = zeros(hans.resultinf(sel).dim); cbild = [];
    aus = find((hans.aXYZ(i,aact)==hans.aXYZ(i,:))&(hans.aXYZ(anot(1),aact)==hans.aXYZ(anot(1),:))); cross(aus) = 1;
    aus = find((hans.aXYZ(i,aact)==hans.aXYZ(i,:))&(hans.aXYZ(anot(2),aact)==hans.aXYZ(anot(2),:))); cross(aus) = 1;
    aus = find((hans.aXYZ(anot(1),aact)==hans.aXYZ(anot(1),:))&(hans.aXYZ(anot(2),aact)==hans.aXYZ(anot(2),:))); cross(aus) = 0; % macht den einen Mittelvoxel wieder rein
    cross = rot90(reshape(cross(auswahl),hans.resultinf(sel).dim(anot(1)),hans.resultinf(sel).dim(anot(2))),1);
    %if i~=1, cross = fliplr(cross); end;
    %if hans.resultinf(sel).mat(i,i) > 0, cross = fliplr(cross); end;
    %if hans.aXYZdir(i) > 0, cross = fliplr(cross); end;
    if hans.aXYZdir(anot(1)) < 0, cross = fliplr(cross); end;
    if hans.aXYZdir(anot(2)) < 0, cross = flipud(cross); end;
    cbild(:,:,1) = cross*cross_farbe(1);cbild(:,:,2) = cross*cross_farbe(2);cbild(:,:,3) = cross*cross_farbe(3);
    ci = image(cbild,'Parent',hans.cross(i));
    set(ci,'ButtonDownFcn',sprintf('nix_Refresh(%d)',i),'AlphaData',cross);
    axis(hans.cross(i),'off','equal');
    if hans.resultinf(sel).dim(anot(1))<re_verkl, m = get(hans.cross(i),'XLIM'); set(hans.cross(i),'XLIM', m + (re_verkl-diff(m))/2* [-1,1]); end;
    if hans.resultinf(sel).dim(anot(2))<re_verkl, m = get(hans.cross(i),'YLIM'); set(hans.cross(i),'YLIM', m + (re_verkl-diff(m))/2* [-1,1]); end;
end;

% Colourbar mit eigener Unterlage malen
% Unterlage
unt = ones(100,20);
breite = 5;
for k = 1:size(unt,2),
    for j = 1:breite,
        unt(j + mod(ceil((k/breite)),2)*breite : breite*2:end,k) = .3;
    end;
    
end;
unt(:,:,1) = unt; unt(:,:,2) = unt(:,:,1); unt(:,:,3) = unt(:,:,1);
unterlage = image(unt,'Parent',hans.colourbar);
hold(hans.colourbar,'on');
% Colorbar
cb = reshape(ajet,size(ajet,1),1,3);  cb = repmat(cb,1,size(unt,2));
cb(:,:,1) = cb(end:-1:1,:,1); cb(:,:,2) = cb(end:-1:1,:,2); cb(:,:,3) = cb(end:-1:1,:,3);
oberlage = image(cb,'Parent',hans.colourbar);

set(hans.colourbar,'XTick',[],'XColor','white');
set(hans.colourbar,'YTick',[1,25,50,75,99],'Ycolor','white');
set(hans.colourbar,'YTickLabel',{functrange(2) : diff(fliplr(functrange))/4 : functrange(1)});
cb_alpha = [functrange(2) : diff(fliplr(functrange))/99 : functrange(1)]';
cb_alpha = (abs(cb_alpha)-amaprange(1))/((amaprange(2)-amaprange(1)));%cb_alpha = (abs(cb_alpha)-imageprops.transp(1))/((imageprops.transp(2)-imageprops.transp(1)));
cb_alpha(find(cb_alpha>1)) = 1; cb_alpha(find(cb_alpha<0)) = 0;
cb_alpha = (cb_alpha + amapint(1)) * (amapint(2) - amapint(1));
cb_alpha = repmat(cb_alpha,1,size(unt,2));
set(oberlage,'AlphaData',cb_alpha);
hold(hans.colourbar,'off');

set(hans.edit(1),'String',num2str(h(1)));
set(hans.edit(2),'String',num2str(h(2)));
set(hans.edit(3),'String',num2str(h(3)));

%% Cluster-Popup auf aktuellen Voxel aktualisieren
cla(hans.raw);
clc;
if sel > 1,
    if hans.clust(sel).clust(aact) == 0,
        set(hans.clusterpopup,'Value',size(hans.clust(sel).out,1)+1);
    else,
        aval = find(hans.clust(sel).out(:,end) == hans.clust(sel).clust(aact));
        set(hans.clusterpopup,'Value',aval);
    end;
    
    %% Rohdaten plotten als Diagramm
    if hans.resultinf(1).private.dat(aact) == 1
        set(hans.raw,'XLim',[0,1],'YLIM',[0,1]);
        text(.5,.5,'No Variance in this Voxel','Parent',hans.raw,'HorizontalAlignment','center','FontSIze',26);
    elseif hans.resultinf(1).private.dat(aact) == 2
        set(hans.raw,'XLim',[0,1],'YLIM',[0,1]);
        text(.5,.5,'Not enough Group Variance in this Voxel','Parent',hans.raw,'HorizontalAlignment','center','FontSIze',26);
        
        %     if hans.novar.private.dat(aact)==1, % Gar keine Varianz vorgegeben
        %         set(hans.raw,'XLim',[0,1],'YLim',[0,1]);
        %         text(.5,.5,'No Variance in this Voxel','Parent',hans.raw,'HorizontalAlignment','center','FontSIze',26);
        %         %    if hans.efilestype(sel)==2, text(.5, 1.065,sprintf('Diff = %1.3f',aresultdata(aact)),'Parent',hans.raw,'HorizontalAlignment','center','FontSize',16); end;
        %     elseif hans.nogroupvar.private.dat(aact)==1, % Nicht genug Gruppenvarianz
        %         set(hans.raw,'XLim',[0,1],'YLim',[0,1]);
        %         text(.5,.5,'No Group Variance in this Voxel','Parent',hans.raw,'HorizontalAlignment','center','FontSIze',26);
        % if hans.efilestype(sel)==2, text(.5, 1.065,sprintf('Diff = %1.3f',aresultdata(aact)),'Parent',hans.raw,'HorizontalAlignment','center','FontSize',16); end;
    else, % Gerechnet result plotting eingeleitet
        clc;
        aeffect = get(hans.effectpopup,'String');
        if hans.efilestype(sel)==1, % um max interaction anzuzeigen bei PostHoc
            aeffect = aeffect{sel};
        else,
            aeffect = aeffect{hans.maxinteraction};
        end;
        
        
        %% Ausgabe
        if LUE.type == 1,
            
            nix_result_posthoc(1,aeffect,aact,sel); % boxplots in figure
            
        elseif LUE.type == 2,
            
            % this voxel
            strvox = sprintf('Voxel %1.1f %1.1f %1.1f',x,y,z);
            LUE.saveposthoc = [LUE.saveposthoc, sprintf('%s\n',repmat('-',length(strvox)+2+16*2,1))];
            LUE.saveposthoc = [LUE.saveposthoc, sprintf('---------------- %s ----------------\n',strvox)];
            LUE.saveposthoc = [LUE.saveposthoc, sprintf('%s\n\n\n',repmat('-',length(strvox)+2+16*2,1))];
            
            adata         = [];
            if any(strfind(aeffect,'ImgData')),
                adata = [adata,inn.working(:,aact)];
            end;
            nix_result_posthoc(1,aeffect,aact,sel,adata); % boxplots in figure
            nix_result_posthoc(2,aeffect,aact,sel,adata); % descriptives + Posthoc-output
            
            % weist die typische Verteilung der Personen zu
            strvox = sprintf('General for Cluster');
            LUE.saveposthoc = [LUE.saveposthoc, sprintf('%s\n',repmat('-',length(strvox)+2+16*2,1))];
            LUE.saveposthoc = [LUE.saveposthoc, sprintf('---------------- %s ----------------\n',strvox)];
            LUE.saveposthoc = [LUE.saveposthoc, sprintf('%s\n\n\n',repmat('-',length(strvox)+2+16*2,1))];
            
            actcluster = hans.clust(sel).clust(aact);
            if ~isequal(actcluster,0),
                weg        = find(hans.clust(sel).clust == actcluster);
                adata      = mode(inn.working(:,weg)*-1,2)*-1;
                
                try,
                    nix_result_posthoc(2,aeffect,aact,sel,adata); % descriptives + Posthoc-output
                catch,
                    LUE.saveposthoc = [LUE.saveposthoc, sprintf('Not possible. The homogeneity of overlaps within this cluster is too small.\n')];
                end;
            else,
                LUE.saveposthoc = [LUE.saveposthoc, sprintf('Voxel outside a cluster was chosen.\n')];
            end;
            
            fprintf('%s',LUE.saveposthoc);
        end;
    end;
else,
    set(hans.raw,'XLim',[0,1],'YLIM',[0,1]);
    if hans.resultinf(sel).private.dat(aact) == 1
        text(.5,.5,'No Variance in this Voxel','Parent',hans.raw,'HorizontalAlignment','center','FontSIze',26);
    elseif hans.resultinf(sel).private.dat(aact) == 2
        text(.5,.5,'Not enough Group Variance in this Voxel','Parent',hans.raw,'HorizontalAlignment','center','FontSIze',26);
    else
        text(.5,.5,'Analysis was performed','Parent',hans.raw,'HorizontalAlignment','center','FontSIze',26);
    end;
end;
hans.sel = sel;