function list = createEuropeanOptionList(fileName,sInitial_t0)


fileLocation = ['utils\InstrumentCreater\TestInstruments\',fileName];
tt = readtable(fileLocation,"FileType","text","Delimiter",",");

for i = 1:height(tt)
    id = tt.id(i);
    name = tt.name{i};
    ticker = tt.ticker{i};
    type = tt.type{i};
    strike = tt.strike(i);

    maturity = datestr(datetime(sInitial_t0) + caldays(tt.maturity(i)),'yyyy-mm-dd');
    euOption = EuropeanOption();
    euOption.setId(id);
    euOption.setTicker(ticker);
    euOption.setType(type);
    euOption.setName(name);
    euOption.setMaturity(maturity);
    euOption.setStrike(strike);
    if contains(fileName,'Heston')
        price = tt.marketPrice(i);
        euOption.setPrice(price);
    end

    list(i) = euOption;
end
end

