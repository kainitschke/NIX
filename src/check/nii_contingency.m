function [chi2, df, p] = nii_contingency(contab),

if length(size(contab)) == 2,
    %if size(contab,2) > size(contab,1),
    %    contab = contab';
    %end;
    [chi2.MAIN, df.MAIN, p.MAIN] = nii_chi2(contab);
   
elseif length(size(contab)) == 3,
    
    [chi2.MAIN, df.MAIN, p.MAIN] = nii_chi2(reshape(contab,size(contab,1),size(contab,2)*size(contab,3))');
    
    [chi2.A, df.A, p.A] = nii_chi2(reshape(sum(contab,2),size(contab,1),size(contab,3)));
    [chi2.B, df.B, p.B] = nii_chi2(reshape(sum(contab,3),size(contab,1),size(contab,2)));
    
    %for i = 1:size(contab,3),
    %    halter( = 
    %end;
    halter = 0;
    for a = 1 : size(contab,1),
        for b = 1 : size(contab,2),
            for c = 1 : size(contab,3),
               halter = halter + (contab(a,b,c) 
            end;
        end;
    end;
    [chi2.AxB, df.AxB, p.AxB] = nii_chi2(reshape(sum(contab,2),size(contab,1),size(contab,3)));
    
else,
    
    error('not implemented');
    
end;