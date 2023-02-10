function res = Smap(maps, x)
   [sx,sy,nv] = size(maps);
   res = zeros(sx,sy,nv);
%      res = zeros(sx,sy);
   for n=1:nv
       res(:,:,n) =  maps(:,:,n).*x;
   end
% res =  maps(:,:,1).*x;
return;
