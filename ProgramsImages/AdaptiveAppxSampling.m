%% Adaptive sampling and approximation example
clearvars
InitializeDisplay
f = @(x) exp(-10*x).*sin(8*x);

%% Plot function
xData = [0:0.1:0.6 0.8:0.1:1]';
fData = f(xData);
xPlot = (0:0.005:1)';
fPlot = f(xPlot);
figure(1);
plot(xPlot,fPlot,xData,fData,'.')
xlabel('\(x\)')
ylabel('\(f(x)\)')
axis([0 1 -0.2 0.4])
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
legend(h([1 3]),{'\(f(x)\)','APP\((f)(x)\)'})
legend('boxoff')
axis([0 1 -0.2 0.4])
print('-depsc','fandDataAndAppx.eps')

%% Next data point based on prediction error
normf = sqrt(coeff'*fData);
RMSPE = real(sqrt(kernel(0,s,theta) - ...
   sum(KPlotData.*(KDataData\KPlotData')',2))) .* normf;
[~,whBad] = max(RMSPE);
xBad = xPlot(whBad);
fBad = f(xBad);
figure(3)
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot, ...
   xPlot,fAppPlot + [-1,1].*RMSPE);
hold on
h = [h; scatter(xBad,fBad,150,MATLABPurple,'filled','d')];
set(h(4:5),'color',MATLABGreen)
xlabel('\(x\)')
lgd = legend(h([1 3 4 6]),{'\(f(x)\)','APP\((f)(x)\)', ...
   'APP\((f)(x) \pm \)RMSPE\((x)\)', ...
   '\((x_{\textrm{bad}},f(x_{\textrm{bad}})\)'})
legend('boxoff')
axis([0 1 -0.2 0.4])
print('-depsc','fandDataAndAppxAndRMSPE.eps')



