function [detvtx,sided_pval,pth,m0] = nix_fdr_stepup(pval,sgn,maskvtx,rate,tail)
% [detvtx,sided_pval,pth,m0] = nix_fdr_stepup2(pval,maskvtx,rate,tail)
%
% Two-stage FDR approach to achieve tighter control of the FDR. This
% procedure is more powerful than the original FDR procedure implemented in
% nix_fdr_stepup.
%
% Input
% pval: P-values.
% sgn: Sign of vertex-wise contrasts.
% maskvtx: Mask's vertices (1-based). Default [] (all vertices included).
% rate: Expected FDR.
% tail: Must be -1 for left-sided or 0 for two-sided or 1 for right-sided
% hypothesis testing.
%
% Output
% detvtx: Detected vertices (1-based).
% spval: Unsigned sided p-values.
% pth: FDR threshold.
% m0: Estimated number of null vertices.
%
% $Revision: 1.1.1.1 $  $Date: 2012/02/02 01:25:52 $
% Original Author: Jorge Luis Bernal Rusiel
% CVS Revision Info:
%    $Author: jbernal$
%    $Date: 2012/02/02 21:25:58 $
%    $Revision: 1.1 $
% References: Benjamini, Y., Krieger, A.M., Yekutieli, D. (2006). Adaptive
% linear step-up procedures that control the false discovery rate.
% Biometrika, 93, 491-507.
%
if nargin < 2
    error('Too few inputs');
elseif nargin < 5
    tail = 0;
    if nargin < 4
        rate = 0.05;
        if nargin < 3
            maskvtx = [];
        end
    end;
end;
nv0 = length(pval);
if isempty(maskvtx)
    maskvtx = 1:nv0;
end;
p = pval(maskvtx).*sgn(maskvtx);
p(pval(maskvtx)==1) = 1;
nv = length(p);

%% First stage (m0 estimation)
q0 = rate/(1+rate);
spval = sidedPval(p,tail);
pth0 = lme_mass_FDR(spval,q0);
if tail==0
    % two-sided thresh
    detv0 = maskvtx(spval <= pth0);
elseif  tail==-1
    % left-sided thresh
    vtx = maskvtx(p<0);
    detv0 = vtx(spval(p<0) <= pth0);
elseif tail==1
    % right-sided thresh
    vtx = maskvtx(p>0);
    detv0 = vtx(spval(p>0) <= pth0);
end;
ndetv0 = length(detv0);
m0 = nv-ndetv0;
%% Second stage
if (ndetv0 ~= 0) && (ndetv0 ~= nv)
    % one sided-thresh
    pth = lme_mass_FDR(spval,q0*nv/m0);
    if tail==0
        % two-sided thresh
        detvtx = maskvtx(spval <= pth);
    elseif  tail==-1
        % left-sided thresh
        detvtx = vtx(spval(p<0) <= pth);
    elseif tail==1
        % right-sided thresh
        detvtx = vtx(spval(p>0) <= pth);
    end;
else
    detvtx = detv0;
    pth = pth0;
end;
sided_pval = ones(1,nv0);
sided_pval(maskvtx) = spval;



%% AUXILIAR function
function [spval] = sidedPval(pval,tail)
spval = abs(pval);
if tail==-1
    spval(pval<0) = spval(pval<0)*0.5;
    spval(pval>0) = 1-0.5*spval(pval>0);
elseif tail==1
    spval(pval>0) = spval(pval>0)*0.5;
    spval(pval<0) = 1-0.5*spval(pval<0);
elseif tail~=0
    error('Tail must be -1 for left-sided or 0 for two-sided or 1 for right-sided hypothesis testing');
end




function pthresh = lme_mass_FDR(p,fdr)
% pthresh = lme_mass_FDR(p,fdr)
% This function has slightly modified the freesurfer's fast_fdrthresh
% function (jbernal modification 2010)
%
% p = list of p values between -1 and 1
% fdr = false discovery rate, between 0 and 1
%
% Based on Tom's FDR.m from 
%   http://www.sph.umich.edu/~nichols/FDR/FDR.m
% The threshold returned from this function is based on an 
% assumption of "independence or positive dependence",
% which should be "reasonable for imaging data".
%
% $Id: fast_fdrthresh.m,v 1.1 2004/10/30 00:36:45 greve Exp $
%

if(nargin ~= 2)
  fprintf('pthresh = lme_mass_FDR(p,fdr)\n');
  return;
end
% pthresh = [];
p = sort(abs(p(:)));
Nv = length(p(:));
nn = [1:Nv]';
imax = max(find(p <= fdr*nn/Nv));
if(~isempty(imax))
  %fprintf('imax = %d\n',imax);
  pthresh = p(imax);
else
    %This is just to not return the min(p) in this case
    pthresh = min(p)/10;
end
return;