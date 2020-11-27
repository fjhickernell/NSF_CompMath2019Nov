%% Plot of IID and Sobol Points
function PlotPoints(bwcolor)
gail.InitializeDisplay %clean up 
format long

if ~nargin
   bwcolor = 'color';
end

if strcmp(bwcolor,'bw')
   black = [0 0 0];
   MATLABBlue = black;
   MATLABOrange = black;
   MATLABPurple = black;
   MATLABGreen = black;
   %MATLABCyan = black;
end

%% Plot IID Random Points
rng(2947)
n = 64;
d = 2;
tick = 0:0.25:1;
xIID = rand(n,d);
figure
plot(xIID(:,1),xIID(:,2),'.','color',MATLABBlue)
xlabel('\(x_{i1}\)')
ylabel('\(x_{i2}\)')
title('Uniform IID Data Sites')
axis square
set(gca,'xtick',tick,'ytick',tick)
print -depsc IIDPoints.eps

%% Plot Unscramled Sobol Points
xUSobol = net(sobolset(d),n);
figure
plot(xUSobol(:,1),xUSobol(:,2),'.','color',MATLABOrange)
xlabel('\(x_{i1}\)')
ylabel('\(x_{i2}\)')
title('Unscrambled Sobol'' Nodes')
axis square
set(gca,'xtick',tick,'ytick',tick)
print -depsc USobolPoints.eps

%% Plot Scramled Sobol Points
xSSobol = net(scramble(sobolset(d),'MatousekAffineOwen'),n);
figure
plot(xSSobol(:,1),xSSobol(:,2),'.','color',MATLABPurple)
xlabel('\(x_{i1}\)')
ylabel('\(x_{i2}\)')
title('Uniform Sobol'' Data Sites')
axis square
set(gca,'xtick',tick,'ytick',tick)
print -depsc SSobolPoints.eps

figure
xChebSSobol = (1+sin(pi*(-1/2 + xSSobol)))/2;
figure
plot(xChebSSobol(:,1),xChebSSobol(:,2),'.','color',MATLABPurple)
xlabel('\(x_{i1}\)')
ylabel('\(x_{i2}\)')
title('Chebyshev Sobol'' Data Sites')
axis square
set(gca,'xtick',tick,'ytick',tick)
print -depsc ChebSSobolPoints.eps


%% Plot Lattice Points
rng(47)
xlattice = gail.lattice_gen(1,n,d);
shift = rand(1,d);
sxlat = mod(bsxfun(@plus,xlattice,shift),1); 
figure
plot(sxlat(:,1),sxlat(:,2),'.','color',MATLABGreen)
xlabel('\(x_{i1}\)')
ylabel('\(x_{i2}\)')
title('Uniform Lattice Data Sites')
axis square
set(gca,'xtick',tick,'ytick',tick)
print -depsc ShiftedLatticePoints.eps

figure
xChebLat = (1+sin(pi*(-1/2 + xlattice)))/2;
figure
plot(xChebLat(:,1),xChebLat(:,2),'.','color',MATLABPurple)
xlabel('\(x_{i1}\)')
ylabel('\(x_{i2}\)')
title('Chebyshev Lattice Data Sites')
axis square
set(gca,'xtick',tick,'ytick',tick)
print -depsc ChebLatticePoints.eps


%% Kernel Space Filling Designs
d = 4;
xkernel = zeros(n,d);
xkernel(1,:) = 0.5*ones(1,d);
theta = -ones(1,d);
kernel = @(t,x) MaternKernelOne(t,x);
meval = 16;
xeval = seqFixedDes([1 2^meval], d);
for i = 2:n
   [Kmat, Kdateval, Kdiageval] = KMP(xkernel(1:i-1,:), xeval, kernel);
   [~, ~, whKX] = powerfun(Kmat, Kdateval, Kdiageval);
   xkernel(i,:) = xeval(whKX,:);
end
figure
plot(xkernel(:,1),xkernel(:,2),'.','color',MATLABGreen)
xlabel('\(x_{i1}\)')
ylabel('\(x_{i2}\)')
title('Kernel Based Data Sites')
axis square
set(gca,'xtick',tick,'ytick',tick)
print -depsc KernelSpaceFillingPoints.eps



