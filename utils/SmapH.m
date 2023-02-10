function res = SmapH(maps, x)
    [sx,sy,nv] = size(maps);
%     res = zeros(sx,sy);

   res = zeros(sx,sy,nv);
%     h = conj(maps);
%    for n=1:nv
%     res(:,:,n) = sum(conj( maps(:,:,n)).*x,3);
%    end
res = sum(conj( maps ).*x,3);
return;