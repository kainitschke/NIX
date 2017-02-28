function nix_write_cluster

global LUE hans

aact    = LUE.aact;
sel     = get(hans.effectpopup,'Value');
aeffect = hans.efiles{sel};
pval    = str2num(get(hans.clustui(1),'String'));

to_write = get(hans.clusterpopup,'String');

fid  = fopen(fullfile(LUE.resultdir,sprintf('ClusterData__%1.1f_%1.1f_%1.1f__%s__%1.5f.txt',LUE.coord,aeffect,pval)),'w+');

for i = 1:length(to_write),
    fwrite(fid,sprintf('%s\n',to_write{i}));
end;

fwrite(fid,sprintf('#Effect: %s | p-threshold for clusters: %1.5f',aeffect,pval));

fclose(fid);