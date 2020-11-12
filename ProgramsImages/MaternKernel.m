function [kerval,kdiageval,errKNull,theta] = ...
   MaternKernel(t,x,theta,transYes)
   if nargin < 4, transYes = false; end
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
      normtmx = reshape(sqrt(sum(tmx.^2,3)),[nt,nx]);
      kerval = (1 + normtmx) .*  exp(-normtmx);
      if nargout > 1
         kdiageval = ones(nx,1);
         errKNull = max(kdiageval);
      end
   else
      dd = dth/2;
      if transYes
         theta = [exp(theta(1:dd)) theta(dd+1:dth)];
      end
      thetaa = theta(1:dd);
      thetab = theta(dd+1:dth);
      tmx = (reshape(t,[nt,1,d]) - reshape(x,[1,nx,d])) ...
         .* reshape(thetaa,[1 1 dd]);
      normtmx = reshape(sqrt(sum(tmx.^2,3)),[nt,nx]);
      tmpl = reshape(sum((reshape(t,[nt,1,d]) + reshape(x,[1,nx,d]))...
         .* reshape(thetab,[1 1 dd]),3),[nt,nx]);
      kerval = exp(tmpl).*(1 + normtmx) .*  exp(-normtmx); 
      if nargout > 1
         kdiageval = exp(sum(x.*(2*thetab),2));
         errKNull = max(kdiageval);
      end
   end
end

