function  plotErrMeasurePlots(returnTS,calibratedParams,errMeasureParams)

axisFont = 8;
labelFont = 8;
titleFont = 10;

dates = returnTS.getDates;
JSum = errMeasureParams.Jsum;
zYSum = errMeasureParams.zYsum;
zVSum = errMeasureParams.zVsum;
varianceSum = errMeasureParams.varianceSum;
sampleSize = errMeasureParams.iterCount;

expZY = zYSum/sampleSize;
expZV = zVSum/sampleSize;
expJ = round(JSum/sampleSize);
expVariance = varianceSum/sampleSize;

jumpInVol = expJ .* expZV;
jumpInRet = expJ .* expZY;
%% plot jumps in vol and return
figure()

t = tiledlayout(2,1);
title(t,'Jump Plot','fontsize',titleFont,'fontname','Calibri','fontweight','bold')

t1 = nexttile;
plot(dates ,jumpInRet,'-','linewidth',0.75);
datetick('x', 'mmm yyyy', 'keepticks')
xlabel(t1,'Date','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
title(t1,'Jumps in Returns')
%legend('\Return Jumps','fontsize',labelFont,'fontname','Calibri','fontweight','bold','Location','southwest')


t2 = nexttile;
plot(dates ,jumpInVol,'-','linewidth',0.75);
title(t2,'Jumps in Volatility')
datetick('x', 'mmm yyyy', 'keepticks')
xlabel(t2,'Dates','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
%legend('Vol Jumps','fontsize',labelFont,'fontname','Calibri','fontweight','bold','Location','southwest')

%% plot vol
figure()

t = tiledlayout(1,1);
title(t,'Volatility Plot','fontsize',titleFont,'fontname','Calibri','fontweight','bold')

t1 = nexttile;
plot(dates ,expVariance,'-','linewidth',0.75);
datetick('x', 'mmm yyyy', 'keepticks')
xlabel(t1,'Date','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
title(t1,'Volatility')

% t2 = nexttile;
% plot(dates ,expZV,'-','linewidth',0.75);
% datetick('x', 'mmm yyyy', 'keepticks')
% xlabel(t1,'Date','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
% set(gca,'fontsize',axisFont);
% title(t2,'expZV')
% 
% t3 = nexttile;
% plot(dates ,expZY,'-','linewidth',0.75);
% datetick('x', 'mmm yyyy', 'keepticks')
% xlabel(t1,'Date','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
% set(gca,'fontsize',axisFont);
% title(t3,'expZY')



%% plot residuals + garch residuals

% compute SVCJ res
sigma=expVariance.^.5;
returns = returnTS.getValues;
muY = calibratedParams.muY;
resReturns = (returns(2:end)-muY-jumpInRet(2:end))./sigma(1:end-1);
MSE_SVCJ = mse(expVariance,returns.^2);


% compute GARCH res
Mdl = garch(1,1);
EstMdl = estimate(Mdl,returns);
 
vG = infer(EstMdl,returns);
MSE_G = mse(vG,returns.^2);
Res_G = returns./sqrt(vG);

figure()

t = tiledlayout(1,3);
title(t,'Residuals','fontsize',titleFont,'fontname','Calibri','fontweight','bold')

t0 = nexttile;
normplot(resReturns);
set(gca,'fontsize',axisFont);
set(gca, 'xlim', [-4 4]);
title(t0,'Log Returns')

t1 = nexttile;
qqplot(resReturns);
set(gca,'fontsize',axisFont);
set(gca, 'xlim', [-4 4]);
title(t1,['SVCJ (MSE = ',num2str(MSE_SVCJ),' )'])


t2 = nexttile;
qqplot(Res_G);
set(gca,'fontsize',axisFont);
set(gca, 'xlim', [-4 4]);
title(t2,['GARCH (MSE = ',num2str(MSE_G),' )'])
end

