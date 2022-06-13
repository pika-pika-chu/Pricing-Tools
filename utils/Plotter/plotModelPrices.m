function plotModelPrices(priceTable,futsPriceTable)
uniqueMaturities = (unique(priceTable.maturity));

for i = 1:length(uniqueMaturities)
[idx,index] = ismember(priceTable.maturity,repmat(uniqueMaturities(i),length(priceTable.maturity),1),'rows','legacy');
tempTable = priceTable(find(idx),:);
tempTable = sortrows(tempTable,"strike");

% get atm
 [~,index] = ismember(uniqueMaturities(i),futsPriceTable.date);
if ~isempty(index)
     sInitial = futsPriceTable.sInitial(index);

end

figure()
hold on 
plotTitle = ['Maturity ',char(uniqueMaturities(i))];
get(gca,'fontname') ; % shows you what you are using.
set(gca,'fontname','Calibri')  % Set it to times
set(gca,'fontsize',8);
ytickformat("usd");
ytickformat('$%,.0f');
xtickformat("usd");
xtickformat('$%,.0f');
xlabel('K');
ylabel('Call Price');


plot(tempTable.strike,tempTable.hestonPrice,'-o','MarkerSize',2,'LineWidth',1.2,'Color',[0.0 0.0 0.6],'MarkerFaceColor',[0.0 0.0 0.6])
plot(tempTable.strike, tempTable.batesPrice,'-o','MarkerSize',2,'LineWidth',1.2,'Color',[0.6 0.0 0.0],'MarkerFaceColor',[0.6 0.0 0.0])
plot(tempTable.strike,[tempTable.marketQuote],'s','MarkerSize',5,'LineWidth',1.2,'Color',[0 0.6 0.0],'MarkerFaceColor',[0 0.6 0.0])
legend('Heston','Bates','Market ','fontsize',8,'fontname','Calibri','fontweight','bold')
title(plotTitle)

hold off

end

end

