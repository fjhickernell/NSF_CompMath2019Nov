function [kerval,kdiageval,errKNull,theta] = ...
   GaussKernel(t,x,theta,transYes)
   [nt,d] = size(t);
   [nx,dx] = size(x);
      dth = length(theta);
   if ~(dx == d) || ~((mod(dth, d) == 0) || (dth == 1))
      error(['dim of t = ' int2str(d) ...
      ', dim of x = ' int2str(dx) ...
      ', and dim of theta = ' int2str(dth) ...
      ' do not all match'])
   end
   if (dth == d)
      if transYes
         theta = exp(theta);
      end
      tmx = (reshape(t,[nt,1,d]) - reshape(x,[1,nx,d])) ...
         .* reshape(theta,[1 1 dth]);
      normtmx2 = reshape(sum(tmx.^2,3),[nt,nx]);
      kerval = exp(-normtmx2);
      if nargout > 1
         kdiageval = ones(nx,1);
         errKNull = max(kdiageval);
      end
   end
end

