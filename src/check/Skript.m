
y = round(randn(1,40)*20);
A = ceil(rand(1,40)*2);
B = ceil(rand(1,40)*2);

y = y+(A*15);

fprintf('y=c(');fprintf('%d,',[y,y]); fprintf('\b)\n');
fprintf('A=c(');fprintf('%d,',[A,A]); fprintf('\b)\n');
fprintf('B=c(');fprintf('%d,',[B,B]); fprintf('\b)\n');
fprintf('time=c(');fprintf('%d,',[repmat(1,length(y),1),repmat(2,length(y),1)]); fprintf('\b)\n');
fprintf('subj=c(');fprintf('%d,',[1:length(y),1:length(y)]); fprintf('\b)\n');
fprintf('nparLD(y ~ A*B*time, subject=subj)\n');

nix_bt2(y,A,B)