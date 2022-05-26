function  plotOptionPriceGrid(priceTable,sInitial_t0)
figure()
strikes = unique(priceTable.strike);
nosStrikes = length(strikes);
maturities = unique(priceTable.maturity)';
maturitiesInDays = abs(daysact( datestr(sInitial_t0),datestr(maturities)));
nosCol = length(maturities);
prices = priceTable.price;

func =  @(x) strcat("K = ",x);
rowNames = func(string(strikes));

func =  @(x) strcat(x," Day");
colNames = func(string(maturitiesInDays));
format longG 
prices = reshape(prices,nosStrikes,nosCol);
prices = arrayfun(@(x) strcat('$ ',sprintf('%.2f\n',x)),prices,'un',0);
T = table(prices,'RowNames',rowNames);
uitable('Data',T{:,:},'ColumnName',colNames,...
    'RowName',T.Properties.RowNames,'Units', 'Normalized', 'Position',[0, 0, 1, 1]);

end

