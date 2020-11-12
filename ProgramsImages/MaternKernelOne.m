function [kerval,kdiageval,errKNull] = ...
   MaternKernelOne(t,x)
   [nt,d] = size(t);
   [nx,dx] = size(x);
   tmx = (reshape(t,[nt,1,d]) - reshape(x,[1,nx,d]));
   normtmx = reshape(sqrt(sum(tmx.^2,3)),[nt,nx]);
   kerval = exp(-normtmx);
   if nargout > 1
      kdiageval = ones(size(x,1),1);
      errKNull = max(kdiageval);
   end
end

