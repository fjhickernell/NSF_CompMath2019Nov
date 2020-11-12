function [errKXx,errKX,whKX,errPXx,errPX] = powerfun(Kmat, Kdateval, Kdiageval, Peval, PTKinv)
errKXx = sqrt(abs(Kdiageval - sum(Kdateval .* (Kmat\Kdateval),1)'));
[errKX,whKX] = max(errKXx);
if nargin >= 5
   PTKinvKdateval = PTKinv*Kdateval;
   errPXx = sqrt(abs(sum((Peval-PTKinvKdateval').^2,2)));
   errPX = max(errPXx);
end

   
