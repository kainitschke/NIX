function [columnnames,result] = nix_txteinlesen(txtfile)
global LUE

if isequal(txtfile(end-2:end),'xls') | isequal(txtfile(end-3:end),'xlsx'),
    [N,T,R] = xlsread(txtfile); clear N T;
end;

eval(sprintf('ch1 = nix_menu2(''Selection Variable:'',%s''None-->All rows are chosen'');',sprintf('''%s'',',R{1,:})));

if ch1<=size(R,2),
    groups = unique([R{2:end,ch1}]);
    list_groups = sprintf('%d,',groups); list_groups(end) = [];
    eval(sprintf('ch2 = nix_menu2(''Which Group:'',%s);',list_groups));
    ch2 = find(ismember([R{2:end,ch1}],groups(ch2))) + 1;
else,
    ch2 = [2:size(R,1)];
    ch1 = [];
end;


ch3 = []; K = [];
for i = 1:LUE.wtall, 
    p = R(1,:); p(ch1) = []; p(K) = []; %p(ch3) = [];
    pstr = sprintf('''%s'',',p{:}); pstr(end) = [];
    
    %eval(sprintf('K = menu2(''Within-Variable #%d:'',%s);',i,pstr));
    astr = [];
    for j = 1:size(LUE.wtgroups,2),
        astr = sprintf('%s, %s:%d',astr,LUE.wtnames{j},LUE.wtgroups(i,j));
    end;
    %eval(sprintf('K = menu2(''%s #%d:'',%s);',LUE.wtnames{LUE.wtgroups(i,1)},LUE.wtgroups(i,2),pstr));
    eval(sprintf('K = nix_menu2(''%s'',%s);',astr(3:end),pstr));
    ch3 = [ch3; find(ismember([R(1,:)],p(K(end))))];
    if ~isempty(ch1),
        K = [ch3(ch1>ch3);ch3(ch1<ch3)-1];
    else,
        K = ch3;
    end;
end;

result = cell2mat((R(ch2,ch3)));
columnnames = R(1,ch3);