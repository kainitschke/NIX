function result = nix_f2ldf1(data)
% data = Mx(2+N) matrix with M Subjects and N = WithinVariable. 
% The first and second columns contain the groups

groups2 = data(:,2); [~,b] = sort(groups2); data = data(b,:);
groups1 = data(:,1); [~,b] = sort(groups1); data = data(b,:);
groups1 = data(:,1); groups2 = data(:,2); data = data(:,3:end);
t = size(data,2);
n = size(data,1); 
N = n * t;

tg1 = unique(groups1); gr1 = length(tg1); for i = 1:length(tg1), tgn1(i) = sum((groups1==tg1(i))); gindex1{i} = find(groups1==tg1(i)); end; %tgn1 = tgn1(end:-1:1); 
tg2 = unique(groups2); gr2 = length(tg2); for i = 1:length(tg2), tgn2(i) = sum((groups2==tg2(i))); gindex2{i} = find(groups2==tg2(i)); end; %tgn2 = tgn2(end:-1:1);
tg12 = unique([groups1,groups2],'rows'); gr12 = size(tg12,1); for i = 1:gr12, tgn12(i) = sum((groups1==tg12(i,1)) .* (groups2==tg12(i,2))); gindex12{i} = find((groups1==tg12(i,1)) .* (groups2==tg12(i,2))); end;

allrank = nix_rank(reshape(data,N,1)');
allrank_res = reshape(allrank,n,t); 

Ra = {};
for i = 1:gr1, Ra{1}{1}{i} = reshape(allrank_res(gindex1{i},:),1,tgn1(i)*t); 
for j = 1:gr2, Ra{1}{2}{i}{j} = allrank_res(gindex12{(i-1)*gr2+j},:); end; end;
%for i = 1:gr12, Ra{1}{2}{i} = allrank_res(gindex12{i},:); end;
for i = 1:gr1, Ra{1}{3}{i} = allrank_res(gindex1{i},:); end;
for i = 1:gr2, Ra{2}{1}{i} = reshape(allrank_res(gindex2{i},:),1,tgn2(i)*t); end;
%for i = 1:gr12, Ra{2}{2}{i} = allrank_res(gindex12{i},:); end;
for i = 1:gr2, Ra{2}{2}{i} = allrank_res(gindex2{i},:); end;

holder1 = cellfun(@mean, Ra{1}{1}, 'UniformOutput', false);
holder2 = cellfun(@mean, Ra{2}{1}, 'UniformOutput', false);
holder3 = cellfun(@mean, [Ra{1}{2}{:}], 'UniformOutput', false);
holder4 = cellfun(@mean, holder3, 'UniformOutput', false);
holder5 = cellfun(@mean, Ra{1}{3}, 'UniformOutput', false);
holder6 = cellfun(@mean, Ra{2}{2}, 'UniformOutput', false);
Rmeans = [holder1{:}, holder2{:}, mean(allrank_res), ...
    holder4{:}, ...
    holder5{:}, ...
    holder6{:}, ...
    holder3{:}];
Ros    = [cellfun('size', Ra{1}{1}, 2),  cellfun('size', Ra{2}{1}, 2), repmat(size(allrank_res,1),1,size(allrank_res,2)), ...
    cellfun(@prod, cellfun(@size, [Ra{1}{2}{:}], 'UniformOutput', false)), ...
    reshape(repmat(cell2mat(cellfun(@length, [Ra{1}{3}], 'UniformOutput', false)),t,1),1,t*gr1), ...
    reshape(repmat(cell2mat(cellfun(@length, [Ra{2}{2}], 'UniformOutput', false)),t,1),1,t*gr2), ...
    reshape(repmat(cell2mat(cellfun(@length, [Ra{1}{2}{:}], 'UniformOutput', false)),t,1),1,gr12*t)]; 

RTE = (Rmeans - .5) / N;
pvec = RTE(end-gr1*gr2*t + 1: end);

V = nix_covr('f2ldf1', gr1, gr2, t, Ra{1}{2}, holder3, n);

Pgr1 = eye(gr1) - ones(gr1,gr1)/gr1;
Pgr2 = eye(gr2) - ones(gr2,gr2)/gr2;
Pt   = eye(t) - ones(t,t)/t;
gr11 = ones(1,gr1) / gr1;
gr21 = ones(1,gr2) / gr2;
t1    = ones(1,t) / t;
Cgr1 = kron(Pgr1, kron(gr21, t1));
Cgr2 = kron(gr11, kron(Pgr2, t1));
Ct = kron(gr11, kron(gr21, Pt));
Cgr1gr2 = kron(Pgr1, kron(Pgr2, t1));
Cgr1t = kron(Pgr1, kron(gr21, Pt));
Cgr2t = kron(gr11, kron(Pgr2, Pt));
Cgr1gr2t = kron(Pgr1, kron(Pgr2, Pt));

%% Wald statistics
result.wald.Fgr1     = n * (Cgr1 * pvec')'     * pinv(Cgr1     * V * Cgr1')     * (Cgr1 * pvec');
result.wald.Fgr2     = n * (Cgr2 * pvec')'     * pinv(Cgr2     * V * Cgr2')     * (Cgr2 * pvec');
result.wald.Ft1       = n * (Ct * pvec')'       * pinv(Ct       * V * Ct')       * (Ct   * pvec');
result.wald.Fgr1gr2  = n * (Cgr1gr2 * pvec')'  * pinv(Cgr1gr2  * V * Cgr1gr2')  * (Cgr1gr2 * pvec');
result.wald.Fgr1t1    = n * (Cgr1t * pvec')'    * pinv(Cgr1t    * V * Cgr1t')    * (Cgr1t * pvec');
result.wald.Fgr2t1    = n * (Cgr2t * pvec')'    * pinv(Cgr2t    * V * Cgr2t')    * (Cgr2t * pvec');
result.wald.Fgr1gr2t1 = n * (Cgr1gr2t * pvec')' * pinv(Cgr1gr2t * V * Cgr1gr2t') * (Cgr1gr2t * pvec');

result.wald.dfgr1     = rank(Cgr1*V*(Cgr1)');
result.wald.dfgr2     = rank(Cgr2*V*(Cgr2)');
result.wald.dft1       = rank(Ct*V*(Ct)');
result.wald.dfgr1gr2  = rank(Cgr1gr2*V*(Cgr1gr2)');
result.wald.dfgr1t1    = rank(Cgr1t*V*(Cgr1t)');
result.wald.dfgr2t1    = rank(Cgr2t*V*(Cgr2t)');
result.wald.dfgr1gr2t1 = rank(Cgr1gr2t*V*(Cgr1gr2t)');

result.wald.pgr1     = 1 - cdf('chi2',result.wald.Fgr1, result.wald.dfgr1, Inf);
result.wald.pgr2     = 1 - cdf('chi2',result.wald.Fgr2, result.wald.dfgr2, Inf);
result.wald.pt1       = 1 - cdf('chi2',result.wald.Ft1, result.wald.dft1, Inf);
result.wald.pgr1gr2  = 1 - cdf('chi2',result.wald.Fgr1gr2, result.wald.dfgr1gr2, Inf);
result.wald.pgr1t1    = 1 - cdf('chi2',result.wald.Fgr1t1, result.wald.dfgr1t1, Inf);
result.wald.pgr2t1    = 1 - cdf('chi2',result.wald.Fgr2t1, result.wald.dfgr2t1, Inf);
result.wald.pgr1gr2t1 = 1 - cdf('chi2',result.wald.Fgr1gr2t1, result.wald.dfgr1gr2t1, Inf);

%% Anova statistics
Btgr1     = Cgr1'     * pinv(Cgr1     * Cgr1')     * Cgr1;
Btgr2     = Cgr2'     * pinv(Cgr2     * Cgr2')     * Cgr2;
Btt       = Ct'       * pinv(Ct       * Ct')       * Ct;
Btgr1gr2  = Cgr1gr2'  * pinv(Cgr1gr2  * Cgr1gr2')  * Cgr1gr2;
Btgr1t    = Cgr1t'    * pinv(Cgr1t    * Cgr1t')    * Cgr1t;
Btgr2t    = Cgr2t'    * pinv(Cgr2t    * Cgr2t')    * Cgr2t;
Btgr1gr2t = Cgr1gr2t' * pinv(Cgr1gr2t * Cgr1gr2t') * Cgr1gr2t;
TVgr1     = Btgr1 * V;
TVgr2     = Btgr2 * V;
TVt       = Btt * V;
TVgr1gr2  = Btgr1gr2 * V;
TVgr1t    = Btgr1t * V;
TVgr2t    = Btgr2t * V;
TVgr1gr2t = Btgr1gr2t * V;

result.anova.Fgr1     = (n/sum(diag(TVgr1))) * (pvec * Btgr1 * pvec');
result.anova.Fgr2     = (n/sum(diag(TVgr2))) * (pvec * Btgr2 * pvec');
result.anova.Ft1       = (n/sum(diag(TVt))) * (pvec * Btt * pvec');
result.anova.Fgr1gr2  = (n/sum(diag(TVgr1gr2))) * (pvec * Btgr1gr2 * pvec');
result.anova.Fgr1t1    = (n/sum(diag(TVgr1t))) * (pvec * Btgr1t * pvec');
result.anova.Fgr2t1    = (n/sum(diag(TVgr2t))) * (pvec * Btgr2t * pvec');
result.anova.Fgr1gr2t1 = (n/sum(diag(TVgr1gr2t))) * (pvec * Btgr1gr2t * pvec');

result.anova.dfgr1     = (sum(diag(Btgr1 * V))^2)     / (sum(diag(Btgr1 * V * Btgr1 * V)));
result.anova.dfgr2     = (sum(diag(Btgr2 * V))^2)     / (sum(diag(Btgr2 * V * Btgr2 * V)));
result.anova.dft1       = (sum(diag(Btt * V))^2)       / (sum(diag(Btt * V * Btt * V)));
result.anova.dfgr1gr2  = (sum(diag(Btgr1gr2 * V))^2)  / (sum(diag(Btgr1gr2 * V * Btgr1gr2 * V)));
result.anova.dfgr1t1    = (sum(diag(Btgr1t * V))^2)    / (sum(diag(Btgr1t * V * Btgr1t * V)));
result.anova.dfgr2t1    = (sum(diag(Btgr2t * V))^2)    / (sum(diag(Btgr2t * V * Btgr2t * V)));
result.anova.dfgr1gr2t1 = (sum(diag(Btgr1gr2t * V))^2) / (sum(diag(Btgr1gr2t * V * Btgr1gr2t * V)));

result.anova.pgr1     = 1 - fcdf(result.anova.Fgr1,     result.anova.dfgr1,     Inf);
result.anova.pgr2     = 1 - fcdf(result.anova.Fgr2,     result.anova.dfgr2,     Inf);
result.anova.pt1       = 1 - fcdf(result.anova.Ft1,       result.anova.dft1,       Inf);
result.anova.pgr1gr2  = 1 - fcdf(result.anova.Fgr1gr2,  result.anova.dfgr1gr2,  Inf);
result.anova.pgr1t1    = 1 - fcdf(result.anova.Fgr1t1,    result.anova.dfgr1t1,    Inf);
result.anova.pgr2t1    = 1 - fcdf(result.anova.Fgr2t1,    result.anova.dfgr2t1,    Inf);
result.anova.pgr1gr2t1 = 1 - fcdf(result.anova.Fgr1gr2t1, result.anova.dfgr1gr2t1, Inf);
