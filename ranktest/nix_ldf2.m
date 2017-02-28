function result = nix_ldf2(data,t1,t2)
% data = MxN matrix with M Subjects and N = t1*t2. T1 is nested in t2
% Example for nesting:
% T2 1 1 1 2 2 2
% T1 1 2 3 1 2 3

s = size(data,1); n = s;
N = s * t1 * t2;

allrank = nix_rank(reshape(data,N,1)');
allrank_res = reshape(allrank,n,t1*t2); 
allrank_res_ch = []; for i = 1:t1, allrank_res_ch = [allrank_res_ch, allrank_res(:,([0:t2-1])*t1+i)]; end;

R1 = reshape(allrank_res,n*t1,t2);
R2 = []; for i = 1:t1, R2 = [R2, reshape(allrank_res(:,([0:t2-1])*t1+i),n*t2,1)]; end;

rankmeans = mean(allrank_res);
rankmeans_ch = mean(allrank_res_ch);
R1means = mean(R1);
R2means = mean(R2);
Rallmeans = [rankmeans_ch, R2means, R1means];

RTE = (Rallmeans(1:t1*t2) - 0.5)/N;
V = nix_covr('ldf2', t1 * t2, n, ones(t1*t2)*n, allrank_res_ch);

Ct1               = eye(t1) - ones(t1,t1) / t1; Ct1 = kron(Ct1,ones(1,t2)/t2);
Ct2               = eye(t2) - ones(t2,t2) / t2; Ct2 = kron(ones(1,t1)/t1,Ct2);
Ct1t2              = kron((eye(t1,t1) - ones(t1,t1)/t1), eye(t2) - ones(t2,t2) / t2);

%% Wald Statistic
result.wald.Wt1   = n * (Ct1 * RTE')'  * pinv(Ct1  * V' * Ct1')  * (Ct1 * RTE');
result.wald.Wt2   = n * (Ct2 * RTE')'  * pinv(Ct2  * V' * Ct2')  * (Ct2 * RTE');
result.wald.Wt1t2  = n * (Ct1t2 * RTE')' * pinv(Ct1t2 * V' * Ct1t2') * (Ct1t2 * RTE');
result.wald.dft1  = rank(Ct1);
result.wald.dft2  = rank(Ct2);
result.wald.dft1t2 = rank(Ct1t2);
result.wald.pt1   =  1 - cdf('chi2',result.wald.Wt1, result.wald.dft1, Inf);
result.wald.pt2   =  1 - cdf('chi2',result.wald.Wt2, result.wald.dft2, Inf);
result.wald.pt1t2 =  1 - cdf('chi2',result.wald.Wt1t2, result.wald.dft1t2, Inf);

%% ANOVA Statistics
Bt1   = Ct1'  * pinv(Ct1  * Ct1')  * Ct1;
Bt2   = Ct2'  * pinv(Ct2  * Ct2')  * Ct2;
Bt1t2 = Ct1t2' * pinv(Ct1t2 * Ct1t2') * Ct1t2;
h1    = Bt1 * V;
h2    = Bt2 * V;
ht1t2 = Bt1t2 * V;
result.anova.Ft1    = n/sum(diag(h1)) * (RTE * Bt1 * RTE');
result.anova.Ft2    = n/sum(diag(h2)) * (RTE * Bt2 * RTE');
result.anova.Ft1t2  = n/sum(diag(ht1t2)) * (RTE * Bt1t2 * RTE');
result.anova.dft1   = sum(diag(Bt1 * V))^2 / sum(diag(Bt1 * V * Bt1' * V'));
result.anova.dft2   = sum(diag(Bt2 * V))^2 / sum(diag(Bt2 * V * Bt2' * V'));
result.anova.dft1t2 = sum(diag(Bt1t2 * V))^2 / sum(diag(Bt1t2 * V * Bt1t2' * V'));
result.anova.pt1    = 1 - fcdf(result.anova.Ft1, result.anova.dft1, Inf); % !!! INF WOULD BE APPROPRIATE
result.anova.pt2    = 1 - fcdf(result.anova.Ft2, result.anova.dft2, Inf); % !!! INF WOULD BE APPROPRIATE
result.anova.pt1t2  = 1 - fcdf(result.anova.Ft1t2, result.anova.dft1t2, Inf); % !!! INF WOULD BE APPROPRIATE