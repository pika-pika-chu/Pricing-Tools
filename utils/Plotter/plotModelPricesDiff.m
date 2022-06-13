function plotModelPricesDiff(priceTable)

figure()
hold on 
plotTitle = ['Model Prices Deltas w/ Market Price'];
get(gca,'fontname') ; % shows you what you are using.
set(gca,'fontname','Calibri')  % Set it to times
set(gca,'fontsize',8);
% ytickformat("usd");
ytickformat('%g%%');
xtickformat("usd");
xtickformat('$%,.0f');
xlabel('K');
ylabel('Call Price Diff');
plot(priceTable.strike,[(priceTable.marketPrice -   priceTable.hestonPrice)./priceTable.marketPrice*100,...
    (priceTable.marketPrice -   priceTable.batesPrice)./priceTable.marketPrice*100,...
    (priceTable.marketPrice -   priceTable.blackPrice)./priceTable.marketPrice*100],'-o','LineWidth',1.5)
legend('Heston \delta','Bates \delta','Black \delta','fontsize',8,'fontname','Calibri','fontweight','bold')
title(plotTitle)

hold off



end

