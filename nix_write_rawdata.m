function nix_write_rawdata,
% writes Infos of the present voxel for further plotting or statistical
% analyses in long format in text file (for SPSS, R, Excel ...); tab 
% delimited; text file is created in the result-folder

global LUE hans inn

aact    = LUE.aact;
sel     = get(hans.effectpopup,'Value');
aeffect = hans.efiles{sel};

fid  = fopen(fullfile(LUE.resultdir,sprintf('RawData__%1.1f_%1.1f_%1.1f__%s.txt',LUE.coord,aeffect)),'w+');

if LUE.type == 1, % cont table
    
    % Header
    fwrite(fid,sprintf('%s\t',LUE.namesfinal{:}));
    fwrite(fid,sprintf('lesion_voxel'));
    fwrite(fid,sprintf('\n'));
    
    % Content
    for i = 1:size(LUE.final,1),
        for j = 1:size(LUE.final,2),
            if j < size(LUE.final,2),
                fwrite(fid,sprintf('%s\t',LUE.levelnames{j}{LUE.final{i,j}}));
            else,
                fwrite(fid,sprintf('%s\t',LUE.final{i,j}));
            end;
        end;
        
        try, 
            V = spm_vol(LUE.final{i,j}); 
        catch,
            errordlg(sprintf('An original file was not found:\n%s\n\nWriting aborted',LUE.final{i,j}),'Error');
            return;
        end;
        fwrite(fid,sprintf('%d\n',round(V.private.dat(aact)))); %LUE.levelnames{j}{LUE.final{i,j}}));
    end;
    
elseif LUE.type == 2, % rep measure stats
    
    actcluster = hans.clust(sel).clust(aact);
    if ~isequal(actcluster,0),
        weg        = find(hans.clust(sel).clust == actcluster);
        adata      = mode(inn.working(:,weg)*-1,2)*-1;
    else,
        adata      = nan(size(inn.working,1),1);
    end;
    
    if LUE.btall > 1,
        astart = 1;
        fwrite(fid,sprintf('%s\t',LUE.namesfinal{1}));
        acor   = 1;
    else,
        astart = 2;
        acor   = 0;
    end;
    
    % Header
    for i = 2 : length(LUE.namesfinal) - 1,
        fwrite(fid,sprintf('%s_raw\t',LUE.namesfinal{i}));
        fwrite(fid,sprintf('%s_ranks\t',LUE.namesfinal{i}));
    end;
    fwrite(fid,sprintf('File\t'));
    fwrite(fid,sprintf('lesion_voxel\tlesion_cluster'));
    fwrite(fid,sprintf('\n'));
    
    % Content
    for i = 1:size(LUE.final,1),
        for j = astart:size(LUE.final,2),
            if j == 1,
               fwrite(fid,sprintf('%s\t',LUE.levelnames{j}{LUE.final{i,j}}));
               
            elseif j < size(LUE.final,2),
                fwrite(fid,sprintf('%1.3f\t',LUE.final{i,j}));
                fwrite(fid,sprintf('%1.1f\t',LUE.ranks(i,j-astart+1-acor)));
                
            else,
                fwrite(fid,sprintf('%s\t',LUE.final{i,j}));
                
            end;
        end;
        fwrite(fid,sprintf('%d\t%d\n',round(inn.working(i,aact)),round(adata(i))));
    end;
    
end;

% Additional Stats
fwrite(fid,sprintf('#Effect: %s | Coordinates: %1.1f %1.1f %1.1f | Statistics: %s',aeffect,LUE.coord,LUE.astat));

% close and finish
fclose(fid);