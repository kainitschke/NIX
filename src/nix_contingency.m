function result = nix_contingency(contab,effects)
% effects are a M x N matrix that only contains ones and zeros. The number
% of columns have to be equal to the number of dimensions of contab. Every
% row of "effects comprises the effect that should be tested among the
% dimensions. 
%
% Example: contab is a 2x2x2 matrix
% effects = [1 0 0; 1 0 1]
% --> this will test for the effect of the 1st dimension and the
% interaction effects of the 1st and the 3rd dimension

%mytable = zeros(2,4,2,3); mytable(1,:,1,1) = [ 62,37,32,70]; mytable(2,:,1,1) = [10,50,90,12]; mytable(1,:,2,1) = [100,79,95,88];
%mytable(2,:,2,1) = [94,62,11,77]; mytable(1,:,1,2) = [69,61,57,68];mytable(2,:,1,2) = [72,14,35,34];mytable(1,:,2,2) = [33,37,47,89];
%mytable(2,:,2,2) = [21,44,26, 5];mytable(1,:,1,3) = [45,91,58,97];mytable(2,:,1,3) = [90,16,45,78];mytable(1,:,2,3) = [96, 6,78,22];mytable(2,:,2,3) = [24,23,96,88];
%contab = mytable;
%disp('ACHTUNG alles davor raus');
% mytable(:,:,1) = [26,10;38,8;4,1;1,0];
% mytable(:,:,2) = [8,8;43,24;14,17;3,7];

dtab = size(contab);
nvar = length(dtab);

for i = 1:size(effects,1), 
    
    margin        = [];
    einsen        = find(effects(i,:)==1);
    nullen        = find(effects(i,:)==0);
    for a = 1 : length(einsen),
        margin{a} = einsen(a);
        for b = 1 : length(nullen),
            margin(a) = {[margin{a}, nullen(b)]};
        end;
    end;

    ncon          = length(margin);
    conf          = zeros(nvar,ncon);
    nmar          = 0;
    
    for k = 1 : length(margin)
        tmp                      = margin{k};
        conf([1:length(tmp)], k) = tmp;
        nmar                     = nmar + prod(dtab(tmp));
    end;
    
    ntab          = numel(contab);
    
    [z,dev,nlast] = loglin(conf,contab,ones(1,prod(size(contab))),nmar,.1,20);
    
    string        = 'fit=reshape(z,';
    for j = 1 : length(dtab),
        string    = sprintf('%s%d,',string,dtab(j));
    end;
    string        = sprintf('%s);',string(1:end-1));
    eval(string);
    
    observed      = contab(:)';
    expected      = fit(:)';
    
    pearson       = sum((observed - expected).^2 ./ expected);
    
    observed      = contab((contab .* fit) > 0)';
    expected      = fit((contab .* fit) > 0)';
    
    lrt           = 2 * sum(observed .* log(observed ./ expected));
    ncon          = length(margin);
    conf          = zeros(nvar,ncon);
    nmar          = 0;
    
    for k = 1 : length(margin)
        tmp                      = margin{k};
        conf([1:length(tmp)], k) = tmp;
        nmar                     = nmar + prod(dtab(tmp));
    end;
    
    ntab          = numel(contab);
    
    [z,dev,nlast] = loglin(conf,contab,ones(1,prod(size(contab))),nmar,.1,20);
    
    string        = 'fit=reshape(z,';
    for j = 1 : length(dtab),
        string    = sprintf('%s%d,',string,dtab(j));
    end;
    string        = sprintf('%s);',string(1:end-1));
    eval(string);
    
    observed      = contab(:)';
    expected      = fit(:)';
    
    pearson       = sum((observed - expected).^2 ./ expected);
    
    observed      = contab((contab .* fit) > 0)';
    expected      = fit((contab .* fit) > 0)';
    
    lrt           = 2 * sum(observed .* log(observed ./ expected));
    
    df            = zeros(1,2^nvar);
    
    for k = 1 : length(margin),
        terms = ll_subset(margin{k});
        for j = 1 : length(terms),
            df(sum(2 .^ (terms{j} - 1))) = prod(dtab(terms{j}) - 1);
        end;
    end;
    
    result(i).lrt = lrt;
    %result(i).pea = pearson;
    result(i).df  = ntab - sum(df) - 1;
    result(i).p   = 1-chi2cdf(lrt, result(i).df);
    df            = zeros(1,2^nvar);
    
    for k = 1 : length(margin),
        terms = ll_subset(margin{k});
        for j = 1 : length(terms),
            df(sum(2 .^ (terms{j} - 1))) = prod(dtab(terms{j}) - 1);
        end;
    end;
    
    result(i).lrt = lrt;
    %result(i).pea = pearson;
    result(i).df  = ntab - sum(df) - 1;
    result(i).p   = 1-chi2cdf(lrt, result(i).df);
    
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = ll_subset(x)

y = cell(1,1);

for i = 1 : length(x),
    for j = 1:length(y),
        y{end+1} = [y{j},x(i)];
    end;
end;
if isempty(y{1}), y = y(2:end); end;