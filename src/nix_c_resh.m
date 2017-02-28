function result = nix_c_resh(tab,amountfac,richtung)

if     richtung == 1,
    
    result     = reshape(tab, prod(amountfac(1:ceil(length(amountfac)/2))), prod(amountfac(ceil(length(amountfac)/2)+1:end)));
    
elseif richtung == 2,
    
    string     = 'result = reshape(tab,';
    for i = 1 : length(amountfac),
        string = sprintf('%s%d,',string,amountfac(i));
    end;
    string     = sprintf('%s);',string(1:end-1));
    eval(string);
    %prod(LUE.amountfactors(1:ceil(end/2))),prod(LUE.amountfactors(ceil(end/2)+1:end)));
    
else,
    error('No direction specified');
end;