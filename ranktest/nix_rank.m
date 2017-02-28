function result = nix_rank(amatrix)
% Funktion f�r Ranktest
% amatrix ist ein 1xM-Vektor
 %rminus(amatrix) % Minimum Ranks
 %rplus(amatrix) % Maximum Ranks
 for i = 1:size(amatrix,1),
    %result(i,:) = 1/2 * (rminus(amatrix(i,:))+rplus(amatrix(i,:))); % Midranks
    result(i,:) = tiedrank(amatrix(i,:));
 end;

% function result = rminus(amatrix)
% % Minimum Ranks
% ranks = (repmat(amatrix,size(amatrix,2),size(amatrix,1))) - (repmat(amatrix',size(amatrix,1),size(amatrix,2)));
% ranks = cminus(ranks);
% result = 1 + sum(ranks);
% 
% function result = rplus(amatrix)
% % Maximum Ranks
% ranks = (repmat(amatrix,size(amatrix,2),size(amatrix,1))) - (repmat(amatrix',size(amatrix,1),size(amatrix,2)));
% ranks = cplus(ranks);
% result = sum(ranks);
% 
% function result = cgen(amatrix)
% result = 1/2 * (cplus(amatrix) + cminus(amatrix));
% 
% function result = cminus(amatrix)
% result = (amatrix > 0);
% 
% function result = cplus(amatrix)
% result = (amatrix >= 0);