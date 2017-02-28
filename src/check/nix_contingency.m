function [stats,p] = nix_contingency(data,exact,contcorr)
% 2x2x3
% data = reshape([1:4],2,2)
% data(:,:,2) = reshape([5:8],2,2)
% data(:,:,3) = reshape([9:12],2,2)

% 3x3x3
% data = reshape([1:9],3,3)
% data(:,:,2) = reshape([10:18],3,3)
% data(:,:,3) = reshape([19:27],3,3)

% 2x2xk möglich
[I,J,K] = size(data);

if (I==2) & (J==2),
    if exact == 0,
        for i = 1:K, s1(:,i) = sum(data(:,:,i),2); end;
        for i = 1:K, s2(:,i) = sum(data(:,:,i),1); end;
        for i = 1:K, n(i) = sum(sum(data(:,:,i))); end;
        delta = sum(reshape([data(1,1,:)],1,K) - s1(1,:) .* s2(1,:) ./ n);
        if (contcorr) & (abs(delta) >= .5),
            yates = .5;
        else,
            yates = 0;
        end;
        stats = (abs(delta) - yates).^2 / sum( prod([s1;s2]) ./ (n.^2 .* (n-1)))
        p = 1-fcdf(stats,1,Inf);
        % Common Odds ration fehlt
        
    else,
        for i = 1:K, mn(:,i) = sum(data(:,:,i),1); end;
        m = mn(1,:); n = mn(2,:);
        
    end;
end;
1