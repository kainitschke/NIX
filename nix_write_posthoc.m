function nix_write_posthoc,
% writes the Output of the Matlab-Window into text file

global LUE hans

aact    = LUE.aact;
sel     = get(hans.effectpopup,'Value');
aeffect = hans.efiles{sel};

fid = fopen(fullfile(LUE.resultdir,sprintf('PostHocTest__%1.1f_%1.1f_%1.1f__%s.txt',LUE.coord,aeffect)),'w+');

fwrite(fid,sprintf('###\t%s\n',aeffect));
fwrite(fid,sprintf('###\t%1.1f\t%1.1f\t%1.1f\n\n\n',LUE.coord));

fwrite(fid, LUE.saveposthoc);

fclose(fid);