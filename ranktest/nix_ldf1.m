function result = nix_ldf1(data)

t = size(data,2);
s = size(data,1); n = s;
N = s * t;

allrank = nix_rank(reshape(data,N,1)');
allrank_res = reshape(allrank,n,t);
    
rankmean = mean(allrank_res);
RTE = (rankmean - 0.5)/N;

C = eye(t,t) - ones(t,t) / t;
V = (n/(N^2)) * nix_covr('ldf1',allrank_res, rankmean, t);
CVC = C*V*C;

p = C * RTE';
%% Wald Statistic
result.wald.Ft1 = n * p' * pinv(CVC) * p;
result.wald.dft1 = sum(diag(CVC * pinv(CVC)));
result.wald.pt1 = 1 - cdf('chi2',result.wald.Ft1,result.wald.dft1, Inf);
%% Hotelling Statistics
result.hotelling.At = result.wald.Ft1 * (n - t + 1)/((t - 1) * (n - 1));
result.hotelling.df1t = t - 1;
result.hotelling.df2t = n - t + 1;
result.hotelling.pt = 1 - fcdf(result.hotelling.At,result.hotelling.df1t,result.hotelling.df2t);
%% ANOVA Statistics
T = C * pinv(C * C) * C;
result.anova.Ft1  = n * p' * T * p/sum(diag(T * V));
result.anova.dft1 = (sum(diag(T * V)))^2/sum(diag(T * V * T * V));
result.anova.pt1  = 1 - fcdf(result.anova.Ft1,result.anova.dft1,Inf);