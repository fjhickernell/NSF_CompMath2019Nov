%% Adaptive sampling and approximation example
clearvars
InitializeDisplay
f = @(x) exp(-10*x).*sin(8*x);

%% Plot function
xData = [0:0.1:0.6 0.8:0.1:1]';
fData = f(xData);
xPlot = (0:0.002:1)';
fPlot = f(xPlot);
figure(1);
plot(xPlot,fPlot,xData,fData,'.')
xlabel('\(x\)')
ylabel('\(f(x)\)')
axis([0 1 -0.2 0.4])
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandData.eps')

%% Compute and plot approximation
s = 1; %scale parameter
theta = 1; %shape parameter
dist = @(x,y) abs(x - y');
kernel = @(z,s,theta) s*(1 + theta*z).*exp(-theta*z);
KDataData = kernel(dist(xData,xData),s,theta);
coeff = KDataData\fData;
KPlotData = kernel(dist(xPlot,xData),s,theta);
fAppPlot = KPlotData*coeff;
figure(2)
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot);
xlabel('\(x\)')
legend(h,{'\(f(x)\)','\(f(x_i)\)','APP\((f,n)(x)\)'})
legend('boxoff')
axis([0 1 -0.2 0.4])
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppx.eps')

%% Compute and plot approximation also using small design
xSmallData = [0 0.4 0.6 1]';
fSmallData = f(xSmallData);
KSmallDataData = kernel(dist(xSmallData,xSmallData),s,theta);
coeffSmall = KSmallDataData\fSmallData;
KPlotSmallData = kernel(dist(xPlot,xSmallData),s,theta);
fAppPlotSmall = KPlotSmallData*coeffSmall;
figure;
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot,xPlot,fAppPlotSmall);
xlabel('\(x\)')
legend(h,{'\(f(x)\)','\(f(x_i)\)','APP\((f,10)(x)\)', 'APP\((f,4)(x)\)'})
legend('boxoff')
axis([0 1 -0.2 0.4])
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppxSmall.eps')


%% Next data point based on prediction error
A = 1.3;
normf = sqrt(coeff'*fData);
RMSPE = real(sqrt(kernel(0,s,theta) - ...
   sum(KPlotData.*(KDataData\KPlotData')',2))) .* normf;
[~,whBad] = max(RMSPE);
xBad = xPlot(whBad);
fBad = f(xBad);
figure
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot, ...
   xPlot,fAppPlot + A*[-1,1].*RMSPE);
hold on
h = [h; scatter(xBad,fBad,200,MATLABPurple,'filled','d')];
set(h(4:5),'color',MATLABGreen)
xlabel('\(x\)')
lgd = legend(h([1:4 6]),{'\(f(x)\)','\(f(x_i)\)','APP\((f,10)(x)\)', ...
   'APP\((f,10)(x) \pm \)ERR\((f,10,x)\)', ...
   '\(\bigl(x_{\textrm{bad}},f(x_{\textrm{bad}})\bigr)\)'});
lgd.NumColumns = 2;
legend('boxoff')
axis([0 1 -0.2 0.4])
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppxAndRMSPE.eps')

whMiss = find((fPlot > fAppPlot + A*RMSPE + 1000*eps) | ...
   (fPlot < fAppPlot - A*RMSPE - 1000*eps))
Miss = [xPlot(whMiss) fPlot(whMiss) fAppPlot(whMiss) + A*[-1 1].*RMSPE(whMiss)]



