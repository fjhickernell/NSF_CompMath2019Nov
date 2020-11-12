function x = seqFixedDes(nrange, d, shift, transform)
%SEQFIXEDDES gives an extensible lattice design that is shifted
defShift = 1/3;
if ~exist('d','var')
   if ~exist('shift','var')
      d = 1;
      shift = defShift;
   end
else
   if ~exist('shift','var')
      shift = defShift*ones(1,d);
   end
end
if ~exist('nrange','var')
   nrange = 1;
end
x = lattice_gen(nrange(1),nrange(end),d);
x = mod(x + shift,1);

if nargin < 4
   transform = '';
end

switch transform
   case 'Chebyshev'
      x = (1+sin(pi*(-1/2 + x)))/2;
end
   


