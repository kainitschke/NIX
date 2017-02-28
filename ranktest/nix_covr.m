function V = nix_covr(base, varargin) 
if isequal(base,'ldf1'),
    %% LD.F1
    ranks = varargin{1}; rmean = varargin{2}; t = varargin{3};
    V = zeros(t,t);
    mm = ones(1,size(ranks,2)) * size(ranks,1);
    for i = 1:t
        V(i, i) = sum(((ranks(:, i)' - rmean(i)).^2)/(mm(i) * (mm(i) - 1)));
        if (i < t)
            %%%j = i + 1
            for j = (i + 1):t
                MM = size(ranks,1);
                K = (MM - 1) * (MM - 1) + MM - 1;
                V(i, j) = sum(repmat((1/K),1,MM) .* ((ranks(:,i) - rmean(i))' .* (ranks(:,j) - rmean(j))'));
                V(j, i) = V(i, j);
            end;
        end;
    end;

elseif isequal(base,'ldf2'),
    %% LD.F2
  facs = varargin{1}; N = varargin{2}; colsum = varargin{3}; allranks = varargin{4};
  rankmeans = mean(allranks);
  
    V = zeros(facs, facs);
    for s1 = 1:facs 
      for s2 = 1:facs 
        if s1 == s2 
          holder   = (allranks(:,s1) - rankmeans(s1)) .^ 2;
          V(s1,s2) = V(s1, s2) + N * sum(holder)./((N*facs)^2 .* N .* (N - 1));
        else,
          holder    = (allranks(:,s1) - rankmeans(s1)) .* (allranks(:,s2) - rankmeans(s2));
          ss        = (N - 1) .* (N - 1) + N - 1;
          V(s1, s2) = V(s1, s2) + N .* sum(holder)./((N*facs)^2 * ss);
        end;
      end;
    end;
    
elseif isequal(base,'f1ldf2'),
    %% F1.LD.F2
    n = varargin{1}; t1 = varargin{2}; t2 = varargin{3}; data = varargin{4}; Rmeans = varargin{5}; tgn = varargin{6}; gr = varargin{7};
    K = zeros(gr*t1*t2,gr*t1*t2);
    
    for i = 1:gr,
        V = zeros(t1*t2,t1*t2);
        for s1 = 1 : t1*t2,
            for s2 = 1 : t1*t2,
                if (s1 == s2), %potientieller Problemherd
                    holder   = (data(sum([1,tgn(1:i-1)]):sum(tgn(1:i)),s1) - Rmeans((i-1)*t1*t2+s1)).^2;
                    V(s1,s2) = V(s1,s2) + tgn(i) * sum(holder) /((n*t1*t2)^2 * tgn(i) * (tgn(i) - 1));
                else,
                    holder   = (data(sum([1,tgn(1:i-1)]):sum(tgn(1:i)),s1) - Rmeans((i-1)*t1*t2+s1)) .* (data(sum([1,tgn(1:i-1)]):sum(tgn(1:i)),s2) - Rmeans((i-1)*t1*t2+s2));
                    ss       = (tgn(i) - 1) * (tgn(i) - 1) + tgn(i)  - 1;
                    V(s1,s2) = V(s1,s2) + tgn(i) * (sum(holder))/((n*t1*t2)^2 * ss);
                end;
            end;
        end;
        
        K((((i-1)*(t1*t2)+1):((i-1)*(t1*t2)+t1*t2)),(((i-1)*(t1*t2)+1):((i-1)*(t1*t2)+t1*t2))) = (n/tgn(i)) * V;
        
    end;
    V = K;
    
elseif isequal(base,'f2ldf1'),
    %% F2.LD.F1
    gr1 = varargin{1}; gr2 = varargin{2}; t = varargin{3}; Ra = varargin{4}; Rmeans = varargin{5}; n = varargin{6}; 
    K = zeros(t*gr1*gr2, t*gr1*gr2);
    
    for i = 1:gr1,
        for j = 1:gr2,
            V = zeros(t,t);
            for s1 = 1:t,
                for s2 = 1:t,
                    if (s1==s2),
                        holder   = (Ra{i}{j}(:,s1) - mean(Ra{i}{j}(:,s1))) .^ 2;
                        V(s1,s2) = V(s1,s2) + size(Ra{i}{j},1) * sum(holder) / ((n*t)^2 * size(Ra{i}{j},1) * (size(Ra{i}{j},1) - 1));
                    else,
                        holder   = (Ra{i}{j}(:,s1) - mean(Ra{i}{j}(:,s1))) .* (Ra{i}{j}(:,s2) - mean(Ra{i}{j}(:,s2)));
                        ss       = (size(Ra{i}{j},1) - 1) * (size(Ra{i}{j},1) - 1) + size(Ra{i}{j},1) - 1;
                        V(s1,s2) = V(s1,s2) + size(Ra{i}{j},1) * sum(holder) / ((n*t)^2 * ss);
                    end;
                end;
            end;
            
            K(((((i - 1) * gr2 + (j - 1)) * t + 1):(((i - 1) * gr2 + (j - 1)) * t + t)), ((((i - 1) * gr2 + (j - 1)) * t + 1):(((i - 1) * gr2 + (j - 1)) * t + t))) = (n / size(Ra{i}{j},1)) * V;
        end;
    end;
    V = K;
    
end;
