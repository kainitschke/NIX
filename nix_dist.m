function [h, aact] = nix_dist(x,y,z)
% find nearest existing voxel of input
% kai.nitschke@uniklinik-freiburg.de

global hans

adists = hans.aXYZ-repmat([x;y;z],1,size(hans.aXYZ,2)); adists = sum(adists.^2);
aact = find(adists==min(adists)); aact = aact(1);
h = round(hans.aXYZ(:,aact)'*10)/10; x=h(1); y=h(2); z=h(3);