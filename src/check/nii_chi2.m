function [chi2,df,p] = nii_chi2(contab),
% Function to check contigency tables for k x 2

Nab = sum(contab);
N = sum(Nab);

halter = 0;
for i = 1 : size(contab,1),
    halter = halter + contab(i,1).^2 / sum(contab(i,:));
end;

chi2 = N^2 / (Nab(1)*Nab(2)) * (halter - Nab(1)^2/N);

df = size(contab,1) - 1;

p = 1 - chi2cdf(chi2,df);