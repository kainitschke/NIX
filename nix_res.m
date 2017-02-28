function res = nix_res(pvertausgleich,pverh),
% search for largest monitor as output for the NIX-Viewer

set(0,'units','pixels');
mm = get(0, 'MonitorPositions');

if size(mm,1) > 1, % check for largest monitor
    sm      = 1; % standard monitor
    msizes  = [mm(:,3)-mm(:,1)+1,  mm(:,4)-mm(:,2)+1];
    weg     = find(max(msizes(:,1)) == msizes(:,1)); weg = weg(1);
    if max(msizes(:,2)) == msizes(weg,2),
        sm  = weg;
    else,
        for i = 1 : size(msizes,1),
            if msizes(i,1) * pverh < msizes(i,2),
                nm(i) = msizes(i,1) * pverh;
            else,
                nm(i) = msizes(i,2);
            end;
        end;
        weg = find(nm == max(nm));
        sm = weg(1);
    end;
    res = mm(sm,:);
else,
    res     = get(0,'screenSize'); 
end;

% if length(res)>2, 
%     res     = res(3:4); 
% end; 
%res(2)      = res(2) - pvertausgleich;
%res(4)      = res(4) - pvertausgleich;
