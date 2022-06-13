function plotModelIVs(volTable,coinName,valuationDate,hestonParams, batesParams)
uniqueMaturities = sort(datetime(unique(volTable.maturity)));
% figure('units','normalized','outerposition',[0 0 1 1])
fh = figure();
fh.WindowState = 'maximized';

t = tiledlayout(3,ceil((length(uniqueMaturities))/3));
title(t,[coinName,' K v/s IV Plot for Valuation Date: ',datestr(valuationDate,'ddmmmyyyy')],'fontsize',12,'fontname','Calibri','fontweight','bold')
subtitle(t,['Heston - ',...
    '  v0: ',num2str(hestonParams.v0),....
    '   \theta: ',num2str(hestonParams.theta),....
    '   \kappa: ',num2str(hestonParams.kappa),....
    '   \rho: ',num2str(hestonParams.rho),....
    '   \sigma: ',num2str(hestonParams.sigma),newline,...
    'Bates - ',...
    '  v0: ',num2str(batesParams.v0),....
    '   \theta: ',num2str(batesParams.vT),....
    '   \kappa: ',num2str(batesParams.kappa),....
    '   \rho: ',num2str(batesParams.rho),....
    '   \sigma: ',num2str(batesParams.sigma),...
    '   \lambda: ',num2str(batesParams.lambda),....
    '   \mu_J: ',num2str(batesParams.muJ),....
    '   \sigma_J: ',num2str(batesParams.sigmaJ)],....
    'fontsize',10,'fontname','Calibri','fontweight','bold');

for i = 1:length(uniqueMaturities)

[idx,index] = ismember(volTable.maturity,repmat(uniqueMaturities(i),length(volTable.maturity),1),'rows','legacy');
tempTable = volTable(find(idx),:);
tempTable = sortrows(tempTable,"strike");

maturityInDays = abs(days365(datenum(valuationDate),datenum(uniqueMaturities(i))));

nexttile
hold on 
plotTitle = ['Maturity Date: ',char(uniqueMaturities(i)),' (',num2str(maturityInDays),'d)'];
get(gca,'fontname') ; % shows you what you are using.
set(gca,'fontname','Calibri')  % Set it to times
set(gca,'fontsize',8);
% ytickformat("usd");
ytickformat('%g%%');
xtickformat("usd");
xtickformat('$%,.0f');
xlabel('K','fontweight','bold');
ylabel('Implied Volatility','fontweight','bold');
ylim([round(min(tempTable.marketVol)*100 - 5),round(max(tempTable.marketVol)*100 + 5)])

hestonRMSE = sqrt(mean((tempTable.marketVol - tempTable.hestonVol).^2));
batesRMSE = sqrt(mean((tempTable.marketVol - tempTable.batesVol).^2));

plot(tempTable.strike,round(tempTable.hestonVol*100,2),'-o','MarkerSize',2,'LineWidth',1.2,'Color',[0.0 0.0 0.6],'MarkerFaceColor',[0.0 0.0 0.6])

plot(tempTable.strike,round(tempTable.batesVol*100,2),'-o','MarkerSize',2,'LineWidth',1.2,'Color',[0.6 0.0 0.0],'MarkerFaceColor',[0.6 0.0 0.0])

plot(tempTable.strike,round(tempTable.marketVol*100,2),'s','MarkerSize',4,'LineWidth',1.2,'Color',[0 0.6 0.0],'MarkerFaceColor',[0 0.6 0.0])

% patch([tempTable.strike fliplr(tempTable.strike)],...
%     [round(tempTable.batesVol*100,2) fliplr(round(tempTable.hestonVol*100,2))], [0.0 0.6 0.1])

legend(['Heston RMSE = ',num2str(round(hestonRMSE*100,2)),' %'],...
    ['Bates RMSE = ',num2str(round(batesRMSE*100,2)),' %'],'Market ','fontsize',7,'fontname','Calibri','fontweight','bold','Location','NorthWest')
title(plotTitle)

hold off
 end
end

