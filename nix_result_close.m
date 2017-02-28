function nix_result_close,

global LUE hans inn

try, delete(hans.fh); end;
try, delete(hans.voxelan1.fh); end;

clearvars -global LUE hans inn