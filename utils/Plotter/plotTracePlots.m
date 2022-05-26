function plotTracePlots(paramSet)


iterCount = paramSet.iterCount;
l = linspace(1,iterCount,iterCount)';

% set axis number font size 
axisFont = 8;
labelFont = 8;
titleFont = 10;
%% plot muY
% dev: set tiled layout
figure()

t = tiledlayout(1,1);
title(t,'Trace Plot','fontsize',titleFont,'fontname','Calibri','fontweight','bold')
t1 = nexttile;

param = paramSet.muY;
plot(l ,[param.evolve,param.meanEvolve]);

xlabel(t1,'Iterations','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
title(t1,'\mu_{y}')
legend('\mu_{y}','\mu_{y} moving average','fontsize',labelFont,'fontname','Calibri','fontweight','bold')

%% plot a,b,k, rho
figure()

t = tiledlayout(2,2);
title(t,'Trace Plot','fontsize',titleFont,'fontname','Calibri','fontweight','bold')

t1 = nexttile;
param = paramSet.alpha;
plot(l ,[param.meanEvolve,param.evolve]);
xlabel(t1,'Iterations','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
title(t1,'\alpha')
legend('\alpha','\alpha moving average','fontsize',labelFont,'fontname','Calibri','fontweight','bold','Location','northeast')


t2 = nexttile;
param = paramSet.beta;
plot(l ,[param.meanEvolve,param.evolve]);
title(t2,'\beta')
xlabel(t2,'Iterations','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
legend('\beta','\beta moving average','fontsize',labelFont,'fontname','Calibri','fontweight','bold','Location','northwest')

t3 = nexttile;
param = paramSet.kappa;
plot(l ,[param.meanEvolve,param.evolve]);
title(t3,'\kappa')
xlabel(t3,'Iterations','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
legend('\kappa','\kappa moving average','fontsize',labelFont,'fontname','Calibri','fontweight','bold','Location','northeast')

t4 = nexttile;
param = paramSet.lambda;
plot(l ,[param.meanEvolve,param.evolve]);
title(t4,'\lambda')
xlabel(t4,'Iterations','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
legend('\lambda','\lambda moving average','fontsize',labelFont,'fontname','Calibri','fontweight','bold','Location','northeast')

t.Padding = 'compact';
t.TileSpacing = 'compact';

%% plot muVJ and sigma VJ
figure()

t = tiledlayout(2,2);
title(t,'Trace Plot','fontsize',titleFont,'fontname','Calibri','fontweight','bold')

t1 = nexttile;
param = paramSet.muV_J;
plot(l ,[param.meanEvolve,param.evolve]);
xlabel(t1,'Iterations','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
title(t1,'\mu^{j}_{v}')
legend('\mu^{j}_{v}','\mu^{j}_{v} moving average','fontsize',labelFont,'fontname','Calibri','fontweight','bold','Location','northeast')


t2 = nexttile;
param = paramSet.sigmaY_J;
plot(l ,[param.meanEvolve,param.evolve]);
xlabel(t2,'Iterations','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
title(t2,'\sigma^{j}_{y}')
legend('\sigma^{j}_{y}','\sigma^{j}_{y} moving average','fontsize',labelFont,'fontname','Calibri','fontweight','bold','Location','northeast')

t3 = nexttile;
param = paramSet.muY_J;
plot(l ,[param.meanEvolve,param.evolve]);
xlabel(t3,'Iterations','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
title(t3,'\mu^{j}_{y}')
legend('\mu^{j}_{y}','\mu^{j}_{y} moving average','fontsize',labelFont,'fontname','Calibri','fontweight','bold','Location','southeast')


t3 = nexttile;
param = paramSet.rhoJ;
plot(l ,[param.meanEvolve,param.evolve]);
title(t3,'\rho_{j}')
xlabel(t3,'Iterations','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
legend('\rho_{j}','\rho_{j} moving average','fontsize',labelFont,'fontname','Calibri','fontweight','bold')

t.Padding = 'compact';
t.TileSpacing = 'compact';

%% Plot sigma V and rho
figure()

t = tiledlayout(2,1);
title(t,'Trace Plot','fontsize',titleFont,'fontname','Calibri','fontweight','bold')

t1 = nexttile;
param = paramSet.sigmaV;
plot(l ,[param.meanEvolve,param.evolve]);
title(t1,'\sigma_{v}')
xlabel(t1,'Iterations','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
legend('\sigma_{v}','\sigma_{v} moving average','fontsize',labelFont,'fontname','Calibri','fontweight','bold','Location','southwest')


t2 = nexttile;
param = paramSet.rho;
plot(l ,[param.meanEvolve,param.evolve]);
title(t2,'\rho')
xlabel(t2,'Iterations','fontsize',labelFont,'fontname','Calibri','fontweight','bold')
set(gca,'fontsize',axisFont);
legend('\rho','\rho moving average','fontsize',labelFont,'fontname','Calibri','fontweight','bold','Location','southwest')

t.Padding = 'compact';
t.TileSpacing = 'compact';

end

