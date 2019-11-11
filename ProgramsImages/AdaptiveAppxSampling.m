%% Adaptive sampling and approximation example
clearvars
InitializeDisplay
f0 = @(x) exp(-6*x).*sin(8*x);
f = @(x) f0(x+0.01);
axisBox = [0 1 -0.2 0.5];

%% Plot function
xData = [0:0.1:0.6 0.8:0.1:1]';
fData = f(xData);
xPlot = (0:0.002:1)';
fPlot = f(xPlot);
figure(1)
plot(xPlot,fPlot,xData,fData,'.')
xlabel('\(x\)')
ylabel('\(f(x)\)')
axis(axisBox)
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
figure
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot);
xlabel('\(x\)')
[lgd,icons] = legend(h,{'\(f(x)\)','\(f(x_i)\)','APP\((f,n)(x)\)'});
legend('boxoff')
axis(axisBox)
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
figure
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot,xPlot,fAppPlotSmall);
xlabel('\(x\)')
legend(h,{'\(f(x)\)','\(f(x_i)\)','APP\((f,10)(x)\)', 'APP\((f,4)(x)\)'})
legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppxSmall.eps')


%% Next data point based on prediction error, no inference
A = 1;
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
[lgd,icons] = legend(h([1:4 6]),{'\(f(x)\)','\(f(x_i)\)','APP\((f,10)(x)\)', ...
   'APP\((f,10)(x) \pm \)ERR\((f,10,x)\)', ...
   '\(\bigl(x_{\textrm{bad}},f(x_{\textrm{bad}})\bigr)\)'});
lgd.NumColumns = 2;
legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppxAndRMSPE.eps')

whMiss = find((fPlot > fAppPlot + A*RMSPE + 1000*eps) | ...
   (fPlot < fAppPlot - A*RMSPE - 1000*eps))
Miss = [xPlot(whMiss) fPlot(whMiss) fAppPlot(whMiss) + A*[-1 1].*RMSPE(whMiss)]

%% Find next data point for optimization
[~,whBadMin] = min(fAppPlot - A.*RMSPE);
xBadMin = xPlot(whBadMin);
fBadMin = f(xBadMin);
[fDataMin,whDataMin] = min(fData);
xDataMin = xData(whDataMin);
figure
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fAppPlot, ...
   xPlot,fAppPlot + A*[-1,1].*RMSPE);
hold on
h = [h; scatter(xBad,fBad,200,MATLABPurple,'filled','d'); ...
   scatter(xDataMin,fDataMin,800,MATLABGreen,'filled','p'); ...
   scatter(xBadMin,fBadMin,200,MATLABCyan,'filled','s')];
set(h(4:5),'color',MATLABOrange)
xlabel('\(x\)')
[lgd,icons] = legend(h([1:4 6:8]),{'\(f(x)\)','\(f(x_i)\)',...
   '\(\widehat{\textrm{APP}}(10)(x)\)', ...
   '\(\widehat{\textrm{APP}}(10)(x) \pm \)ERR\((10,x)\)', ...
   '\(\bigl(x_{11},f(x_{11})\bigr)\) for \(\widehat{\textrm{APP}}\)', ...
   '\(\bigl(x_{\textrm{min}},f(x_{\textrm{min}})\bigr)\)', ... 
   '\(\bigl(x_{11},f(x_{11})\bigr)\) for \(\widehat{\textrm{MIN}}\)'});
lgd.NumColumns = 2;
legend('boxoff')
icons(16).Children.MarkerSize = 15;
icons(17).Children.MarkerSize = 15;
icons(18).Children.MarkerSize = 15;
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppxAndRMSPEAndMin.eps')

%% Infer theta using empirical Bayes
Ktheta = @(logth) kernel(dist(xData,xData),s,exp(logth));
objective = @(K,y) mean(log(eig(K))) + log(y'*(K\y));
figure
logthetaRange = log(0.01):0.005:log(100);
semilogx(exp(logthetaRange), ...
   arrayfun(@(lgth) objective(Ktheta(lgth),fData), logthetaRange))
logthopt = fminbnd(@(logth) objective(Ktheta(logth),fData),-5,5);
thetaopt = exp(logthopt)

%% Compute and plot approximation with optimal theta
theta = thetaopt;
KOptDataData = kernel(dist(xData,xData),s,theta);
coeffOpt = KOptDataData\fData;
KOptPlotData = kernel(dist(xPlot,xData),s,theta);
fOptAppPlot = KOptPlotData*coeffOpt;
figure
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptAppPlot);
xlabel('\(x\)')
legend(h,{'\(f(x)\)','\(f(x_i)\)','APP\((f,n)(x)\)'})
legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndOptAppx.eps')

%% Next data point based on prediction error for optimal theta
AOpt = 1;
normf = sqrt(coeffOpt'*fData);
kernelDiag = @(s,theta) s;
RMSPEOpt = real(sqrt(kernelDiag(s,theta) - ...
   sum(KOptPlotData.*(KOptDataData\KOptPlotData')',2))) .* normf;
[~,whBadOpt] = max(RMSPEOpt);
xBadOpt = xPlot(whBadOpt);
fBadOpt = f(xBadOpt);
figure
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptAppPlot, ...
   xPlot,fOptAppPlot + AOpt*[-1,1].*RMSPEOpt);
hold on
h = [h; scatter(xBadOpt,fBadOpt,200,MATLABPurple,'filled','d')];
set(h(4:5),'color',MATLABGreen)
xlabel('\(x\)')
lgd = legend(h([1:4 6]),{'\(f(x)\)','\(f(x_i)\)','APP\((f,10)(x)\)', ...
   'APP\((f,10)(x) \pm \)ERR\((f,10,x)\)', ...
   '\(\bigl(x_{\textrm{bad}},f(x_{\textrm{bad}})\bigr)\)'});
lgd.NumColumns = 2;
legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppxAndRMSPEOpt.eps')

whMissOpt = find((fPlot > fOptAppPlot + AOpt*RMSPEOpt + 1000*eps) | ...
   (fPlot < fOptAppPlot - AOpt*RMSPEOpt - 1000*eps))
Miss = [xPlot(whMissOpt) fPlot(whMissOpt) fOptAppPlot(whMissOpt) + ...
   AOpt*[-1 1].*RMSPEOpt(whMissOpt)]



%% Infer y-varying kernel using empirical Bayes
kernel = @(z,s,theta,a,x,y) s*exp(a*(x+y')).*(1 + theta*z).*exp(-theta*z);
Ktheta = @(logth,a) kernel(dist(xData,xData),s,exp(logth),a,xData,xData);
objective = @(K,y) mean(log(max(eig(K),100*eps))) + log(y'*(K\y));
% theta = 1;
% aRange = -15:0.01:0;
% plot(aRange, ...
%    arrayfun(@(a) objective(Ktheta(log(theta),a),fData), aRange))
% aopt = fminbnd(@(a) objective(Ktheta(log(theta),a),fData),-15,0)
[bopt,objmin] = fminsearch(@(b) objective(Ktheta(b(1),b(2)),fData),[2,-10])
thetaopt = exp(bopt(1))
aopt = bopt(2)

%% Compute and plot approximation with optimal y-varying
theta = thetaopt;
a = aopt;
KOptyDataData = kernel(dist(xData,xData),s,theta,a,xData,xData);
coeffOpty = KOptyDataData\fData;
KOptyPlotData = kernel(dist(xPlot,xData),s,theta,a,xPlot,xData);
fOptyAppPlot = KOptyPlotData*coeffOpty;
figure
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptyAppPlot);
xlabel('\(x\)')
legend(h,{'\(f(x)\)','\(f(x_i)\)','APP\((f,n)(x)\)'})
legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndOptyAppx.eps')

%% Next data point based on prediction error for y-varying
n = length(xData);
%AOpty = -tinv(0.005,n)/sqrt(length(xData))
AOpty = 1;
normf = sqrt(coeffOpty'*fData);
kernelDiag = @(s,theta,a,x) s*exp(a*(2*x));
RMSPEOpty = real(sqrt(kernelDiag(s,theta,a,xPlot) - ...
   sum(KOptyPlotData.*(KOptyDataData\KOptyPlotData')',2))) .* normf;
[~,whBadOpty] = max(RMSPEOpty);
xBadOpty = xPlot(whBadOpty);
fBadOpty = f(xBadOpty);
figure
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptyAppPlot, ...
   xPlot,fOptyAppPlot + AOpty*[-1,1].*RMSPEOpty);
hold on
h = [h; scatter(xBadOpty,fBadOpty,200,MATLABPurple,'filled','d')];
set(h(4:5),'color',MATLABGreen)
xlabel('\(x\)')
lgd = legend(h([1:4 6]),{'\(f(x)\)','\(f(x_i)\)','APP\((f,10)(x)\)', ...
   'APP\((f,10)(x) \pm \)ERR\((f,10,x)\)', ...
   '\(\bigl(x_{\textrm{bad}},f(x_{\textrm{bad}})\bigr)\)'});
lgd.NumColumns = 2;
legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppxAndRMSPEOpty.eps')

whMissOpty = find((fPlot > fOptyAppPlot + AOpty*RMSPEOpty + 1000*eps) | ...
   (fPlot < fOptyAppPlot - AOpty*RMSPEOpty - 1000*eps))
Miss = [xPlot(whMissOpty) fPlot(whMissOpty) fOptyAppPlot(whMissOpty) + ...
   AOpty*[-1 1].*RMSPEOpty(whMissOpty)]

%% Find next data point for optimization based on prediction error for y-varying
[~,whBadOptyMin] = min(fAppPlot - A.*RMSPEOpty);
xBadOptyMin = xPlot(whBadOptyMin);
fBadOptyMin = f(xBadOptyMin);
figure
h = plot(xPlot,fPlot,xData,fData,'.',xPlot,fOptyAppPlot, ...
   xPlot,fOptyAppPlot + AOpty*[-1,1].*RMSPEOpty);
hold on
h = [h; scatter(xBadOpty,fBadOpty,200,MATLABPurple,'filled','d'); ...
   scatter(xBadOptyMin,fBadOptyMin,200,MATLABCyan,'filled','s')];
set(h(4:5),'color',MATLABGreen)
xlabel('\(x\)')
lgd = legend(h([1:4 6:7]),{'\(f(x)\)','\(f(x_i)\)','APP\((10)(x)\)', ...
   'APP\((10)(x) \pm \)ERR\((10,x)\)', ...
   '\(\bigl(x_{11},f(x_{11})\bigr)\) for APP', ...
   '\(\bigl(x_{11},f(x_{11})\bigr)\) for \(\widehat{\textrm{MIN}}\)'});
lgd.NumColumns = 2;
legend('boxoff')
axis(axisBox)
set(gca,'PlotBoxAspectRatio',[1.5 1 1]);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 1.4*pos(3:4)])
print('-depsc','fandDataAndAppxAndRMSPEOpty.eps')

