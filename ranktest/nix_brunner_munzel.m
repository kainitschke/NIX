function [stats, df, p] = nix_brunner_munzel(a,b)
% a and b vectors

if size(a,1) > size(a,2), a = a'; b = b'; end;

na = length(a);
nb = length(b);

ra = nix_rank(a);
rb = nix_rank(b);
r  = nix_rank([a,b]);

ma = mean(r(1:na));
mb = mean(r(na+1:end));

pst = (mb - (nb + 1)/2) / na;
va = sum((r(1:na)       - ra - ma + (na + 1) / 2) .^ 2) / (na - 1);
vb = sum((r(na+1 : end) - rb - mb + (nb + 1) / 2) .^ 2) / (nb - 1);

stats = na * nb * (mb - ma) / (na + nb) / sqrt(na * va + nb * vb);
df = ((na * va + nb * vb) ^ 2) / (((na * va) ^ 2) / (na - 1) + ((nb * vb) ^ 2) / (nb - 1));
p = 2 * min( cdf('t',abs(stats),df), 1 - cdf('t',abs(stats),df));

