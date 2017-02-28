function result = nii_loglin(contab,margin)

contab = zeros(3,4,5);
contab(:,:,1) = [35,26,36,29;21,12,26,27;9,26,23,38];
contab(:,:,2) = [43,26,89,26;78,42,32,9;34,39,42,51];
contab(:,:,3) = [51,72,62,21;6,52,71,32;35,36,48,21];
contab(:,:,4) = [12,26,53,42;96,75,81,56;32,65,75,63];
contab(:,:,5) = [21,18,24,33;22,56,36,91;32,65,63,14];
margin = [1,2;2,3];%[1:3];

dtab = size(contab);
nvar = length(dtab);
ncon = length(margin);
conf = zeros(nvar,ncon);
nmar = 0;

for i = 1 : length(margin)
    tmp          = margin(:,i)';
    conf(tmp, i) = tmp;
    nmar         = nmar + prod(dtab(tmp));
end;

ntab = prod(size(contab));