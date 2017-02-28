function result = nix_bt2(av,uv1,uv2)

% addpath(genpath('\\AFS\.fbi.ukl.uni-freiburg.de\apps\fbi-matlab\NIX'))

if length(unique(uv1))>2 | length(unique(uv2))>2,
    error('Not for more than two levels');
end;

avcor_int = av;
avcor_main = av;
uv1fac = unique(uv1);
uv2fac = unique(uv2);

for i = 1:length(uv1fac), if ~isequal(uv1fac(i),i), weg = find(uv1 == uv1fac(i)); uv1(weg) = i; end; end;
for j = 1:length(uv2fac), if ~isequal(uv2fac(j),j), weg = find(uv2 == uv2fac(j)); uv2(weg) = j; end; end;

%% Interaction
for i = 1:length(uv1fac), 
    uv1i  = uv1 == uv1fac(i);
    muv1i = mean(av(find(uv1i)));
    for j = 1:length(uv2fac),
        uv2i  = uv2 == uv2fac(j);
        muv2i = mean(av(find(uv2i)));
        weg   = find(uv1i .* uv2i);
        avcor_int(weg) = av(weg) - muv1i - muv2i;
    end;
end;

rank_avcor_int = nix_rank(avcor_int);
[~,an] = anovan(rank_avcor_int,{uv1,uv2},'model','interaction','display','off');

%% Main Effect
sammler = cell(length(uv1fac),length(uv2fac));
for k = 1:length(av),
    sammler{uv1(k),uv2(k)} = [sammler{uv1(k),uv2(k)}; k]; % av(k)
end;

for i = 1:length(uv1fac), 
    for j = 1:length(uv2fac), 
        sammler2(i,j) = mean(av(sammler{i,j}));
    end;
end;

for k = 1:length(av),
    holder = ones(length(uv1fac),length(uv2fac));
    holder(uv1(k),:) = 0; holder(:,uv2(k)) = 0;
    holder(uv1(k),uv2(k)) = 1;
    avcor_main(k) = av(k) - mean(sammler2(find(holder)));
end;

rank_avcor_main = nix_rank(avcor_main);
[~,an2] = anovan(rank_avcor_main,{uv1,uv2},'model','interaction','display','off');

%% Delete
fprintf('\nMain Effects\n');
fprintf('F: %1.3f\n',an2{2,6});
fprintf('p: %1.3f\n\n',an2{2,7});
fprintf('F: %1.3f\n',an2{3,6});
fprintf('p: %1.3f\n\n',an2{3,7});

fprintf('Interaction:\n');
fprintf('F: %1.3f\n',an{4,6});
fprintf('p: %1.3f\n',an{4,7});
