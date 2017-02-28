function result = nix_f1ldf2(data,t1,t2)
% result = nix_f1ldf2(data,t1,t2)
% t1 is the size of the first within factor
% t2 is the size of the second within factor
% data = Mx(1+N) matrix with M Subjects and N = t1*t2. 
% The first column contains the groups
% T1 is nested in t2
% Example for nesting:
% T2 1 1 1 2 2 2
% T1 1 2 3 1 2 3

groups = data(:,1); [a,b] = sort(groups); groups = groups(b); data = data(b,2:end);
s = size(data,1); n = s;
N = s * t1 * t2;
tg = unique(groups); gr = length(tg); for i = 1:length(tg), tgn(i) = sum((groups==tg(i))); gindex{i} = find(groups==tg(i)); end;
%na.a.vec = tgn

allrank = nix_rank(reshape(data,N,1)');
allrank_res = reshape(allrank,n,t1*t2); 
allrank_res_ch = []; for i = 1:t1, allrank_res_ch = [allrank_res_ch, allrank_res(:,([0:t2-1])*t1+i)]; end;
%R1 = Rs |R2 = Rr
R1 = []; for i = 1:t1, R1 = [R1, reshape(allrank_res(:,([0:t2-1])*t1+i),n*t2,1)]; end;
R2 = reshape(allrank_res,n*t1,t2);


for i = 1:length(tg),
    R1g{i} = []; R2g{i} = []; Rallg{i} = []; %mRallg = []; %mR1g = []; mR2g = [];
    Rallg{i} = allrank_res_ch(gindex{i},:);
    for j = 1: size(Rallg{i},2),
        R1g{i}(tgn(i)*mod(j-1,t2)+1: tgn(i)*mod(j-1,t2)+1 + tgn(i) - 1, floor((j-1)/t2)+1) = Rallg{i}(:,j);
        R2g{i}(tgn(i)*floor((j-1)/t2)+1: tgn(i)*floor((j-1)/t2)+1 + tgn(i) - 1, mod(j-1,t2)+1) = Rallg{i}(:,j);

    end;
    Rmeans(i) = mean(reshape(Rallg{i},prod(size(Rallg{i})),1));
end;
h1 = cellfun(@mean, R1g, 'UniformOutput', false);
h2 = cellfun(@mean, R2g, 'UniformOutput', false);
ha = cellfun(@mean, Rallg, 'UniformOutput', false);
Rmeans = [Rmeans, mean(R1), mean(R2), [h1{:}],  mean(allrank_res_ch), [h2{:}], [ha{:}]];

RTE = (Rmeans - .5) / N;

Pgr = amat(gr); % PA
Dgr = ones(1,gr) / gr;  %A1
Pt1 = amat(t1); %PC
Dt1 = ones(1,t1) / t1; %C1
Pt2 = amat(t2); %PT
Dt2 = ones(1,t2) / t2; %T1

Cgr = kron(Pgr, kron(Dt1, Dt2)); %CA
Ct1 = kron(Dgr, kron(Pt1, Dt2)); %CC
Ct2 = kron(Dgr, kron(Dt1, Pt2)); %CT

Cgrt1 = kron(Pgr, kron(Pt1, Dt2));
Ct1t2 = kron(Dgr, kron(Pt1, Pt2));
Cgrt2 = kron(Pgr, kron(Dt1, Pt2));
Cgrt1t2 = kron(Pgr, kron(Pt1, Pt2));
 
V = nix_covr('f1ldf2', n, t1, t2, allrank_res_ch, Rmeans((gr+t1+t2+t1*gr+t1*t2+t2*gr+1):(gr+t1+t2+t1*gr+t1*t2+t2*gr+gr*t1*t2)), tgn, gr);
pvec = RTE((gr+t1+t2+t1*gr+t1*t2+t2*gr+1):(gr+t1+t2+t1*gr+t1*t2+t2*gr+gr*t1*t2));

%% Wald statistics
result.wald.Fgr1     = n * (Cgr * pvec')' * pinv(Cgr * V * (Cgr)') * (Cgr * pvec');
result.wald.Ft1     = n * (Ct1 * pvec')' * pinv(Ct1 * V * (Ct1)') * (Ct1 * pvec');
result.wald.Ft2     = n * (Ct2 * pvec')' * pinv(Ct2 * V * (Ct2)') * (Ct2 * pvec');
result.wald.Fgr1t1   = n * (Cgrt1 * pvec')' * pinv(Cgrt1 * V * (Cgrt1)') * (Cgrt1 * pvec');
result.wald.Ft1t2   = n * (Ct1t2 * pvec')' * pinv(Ct1t2 * V * (Ct1t2)') * (Ct1t2 * pvec');
result.wald.Fgr1t2   = n * (Cgrt2 * pvec')' * pinv(Cgrt2 * V * (Cgrt2)') * (Cgrt2 * pvec');
result.wald.Fgr1t1t2 = n * (Cgrt1t2 * pvec')' * pinv(Cgrt1t2 * V * (Cgrt1t2)') * (Cgrt1t2 * pvec');

result.wald.dfgr1     = rank(Cgr*V*(Cgr)');
result.wald.dft1     = rank(Ct1*V*(Ct1)');
result.wald.dft2     = rank(Ct2*V*(Ct2)');
result.wald.dfgr1t1   = rank(Cgrt1*V*(Cgrt1)');
result.wald.dft1t2   = rank(Ct1t2*V*(Ct1t2)');
result.wald.dfgr1t2   = rank(Cgrt2*V*(Cgrt2)');
result.wald.dfgr1t1t2 = rank(Cgrt1t2*V*(Cgrt1t2)');

result.wald.pgr1     = 1 - cdf('chi2',result.wald.Fgr1, result.wald.dfgr1, Inf);
result.wald.pt1     = 1 - cdf('chi2',result.wald.Ft1, result.wald.dft1, Inf);
result.wald.pt2     = 1 - cdf('chi2',result.wald.Ft2, result.wald.dft2, Inf);
result.wald.pgr1t1   = 1 - cdf('chi2',result.wald.Fgr1t1, result.wald.dfgr1t1, Inf);
result.wald.pgr1t2   = 1 - cdf('chi2',result.wald.Fgr1t2, result.wald.dfgr1t2, Inf);
result.wald.pt1t2   = 1 - cdf('chi2',result.wald.Ft1t2, result.wald.dft1t2, Inf);
result.wald.pgr1t1t2 = 1 - cdf('chi2',result.wald.Fgr1t1t2, result.wald.dfgr1t1t2, Inf);

%% Anova statistics
RTEAov = RTE(gr + t1 + t2 + t1 * gr + t1 * t2 + t2 * gr + 1 : gr + t1 + t2 + t1 * gr + t1 * t2 + t2 * gr + gr * t1 * t2);

Btgr     = Cgr' * (pinv(Cgr * Cgr')) * Cgr;
Btt1     = Ct1' * (pinv(Ct1 * Ct1')) * Ct1;
Btt2     = Ct2' * (pinv(Ct2 * Ct2')) * Ct2;
Btgrt1   = Cgrt1' * (pinv(Cgrt1 * Cgrt1')) * Cgrt1;
Btt1t2   = Ct1t2' * (pinv(Ct1t2 * Ct1t2')) * Ct1t2;
Btgrt2   = Cgrt2' * (pinv(Cgrt2 * Cgrt2')) * Cgrt2;
Btgrt1t2 = Cgrt1t2 * (pinv(Cgrt1t2 * Cgrt1t2')) * Cgrt1t2;

TVgr = Btgr * V;
TVt1 = Btt1 * V;
TVt2 = Btt2 * V;
TVt1t2 = Btt1t2 * V;
TVgrt1 = Btgrt1 * V;
TVgrt2 = Btgrt2 * V;
TVgrt1t2 = Btgrt1t2 * V;

result.anova.Fgr1     = (n/sum(diag(TVgr))) * (RTEAov * Btgr * RTEAov');
result.anova.Ft1     = (n/sum(diag(TVt1))) * (RTEAov * Btt1 * RTEAov');
result.anova.Ft2     = (n/sum(diag(TVt2))) * (RTEAov * Btt2 * RTEAov');
result.anova.Ft1t2   = (n/sum(diag(TVt1t2))) * (RTEAov * Btt1t2 * RTEAov');
result.anova.Fgr1t1   = (n/sum(diag(TVgrt1))) * (RTEAov * Btgrt1 * RTEAov');
result.anova.Fgr1t2   = (n/sum(diag(TVgrt2))) * (RTEAov * Btgrt2 * RTEAov');
result.anova.Fgr1t1t2 = (n/sum(diag(TVgrt1t2))) * (RTEAov * Btgrt1t2 * RTEAov');

result.anova.dfgr1     = (sum(diag(Btgr * V))^2)/(sum(diag(Btgr * V * Btgr * V)));
result.anova.dft1     = (sum(diag(Btt1 * V))^2)/(sum(diag(Btt1 * V * Btt1 * V)));
result.anova.dft2     = (sum(diag(Btt2 * V))^2)/(sum(diag(Btt2 * V * Btt2 * V)));
result.anova.dft1t2   = (sum(diag(Btt1t2 * V))^2)/(sum(diag(Btt1t2 * V * Btt1t2 * V)));
result.anova.dfgr1t1   = (sum(diag(Btgrt1 * V))^2)/(sum(diag(Btgrt1 * V * Btgrt1 * V)));
result.anova.dfgr1t2   = (sum(diag(Btgrt2 * V))^2)/(sum(diag(Btgrt2 * V * Btgrt2 * V)));
result.anova.dfgr1t1t2 = (sum(diag(Btgrt1t2 * V))^2)/(sum(diag(Btgrt1t2 * V * Btgrt1t2 * V)));

result.anova.pgr1     = 1 - fcdf(result.anova.Fgr1, result.anova.dfgr1, Inf);
result.anova.pt1     = 1 - fcdf(result.anova.Ft1, result.anova.dft1, Inf);  
result.anova.pt2     = 1 - fcdf(result.anova.Ft2, result.anova.dft2, Inf);  
result.anova.pt1t2   = 1 - fcdf(result.anova.Ft1t2, result.anova.dft1t2, Inf);
result.anova.pgr1t1   = 1 - fcdf(result.anova.Fgr1t1, result.anova.dfgr1t1, Inf);
result.anova.pgr1t2   = 1 - fcdf(result.anova.Fgr1t2, result.anova.dfgr1t2, Inf);
result.anova.pgr1t1t2 = 1 - fcdf(result.anova.Fgr1t1t2, result.anova.dfgr1t1t2, Inf);
                                            
function result = amat(x),
result = eye(x,x) - ones(x,x) / x;