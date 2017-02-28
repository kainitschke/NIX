function nix_manual,

apath = which('nix.m');
apath = apath(1:end-6);
afile = dir(fullfile(apath,'NIX_Manual*.pdf'));
open(fullfile(apath,afile(end).name));