function nix_write_mult_coord,

global LUE hans

fid = fopen(fullfile(LUE.resultdir,sprintf('%s_MultipleCoordsTests.txt',datestr(now,'yyyymmdd_HHMMSS'))),'a+');
fwrite(fid,hans.SaveMultCond);
fclose(fid);