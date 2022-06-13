clc
clear all

valuationDate = '11-Jun-2022';
coinName = 'ETH';
bidAsk = 'ask';
folderLocation = ['res\',coinName,'\',datestr(valuationDate,'ddmmmyyyy')];
marketInstruments = getDeribitOptionList(folderLocation,valuationDate,bidAsk);


for i = 1:length(marketInstruments)
    instrument = marketInstruments(i);
    id(i,1) = instrument.getId;
    strike(i,1) = instrument.getStrike;
    maturity{i,1} = instrument.getMaturity;
    marketVol(i,1) = instrument.getImpliedVol;
    marketQuote(i,1) = instrument.getPrice;
end

%% BSM
modelName = 'BSM';
blackEuropeanOptionPricer = EuropeanOptionPricer(valuationDate,marketInstruments,modelName);
blackEuropeanOptionPricer.calculatePrice;
blackEuropeanOptionPricer.calculateImpliedVol;
blackPricelist = blackEuropeanOptionPricer.getInstrumentList;

for i = 1:length(marketInstruments)
    instrument = marketInstruments(i);
    blackPrice(i,1) = instrument.getPrice;
    blackVol(i,1) = instrument.getImpliedVol;
    instrument.setPrice(marketQuote(i,1));
    instrument.setImpliedVol(marketVol(i,1));
end

%% Heston
modelName = 'Heston'; 
hestonEuropeanOptionPricer = EuropeanOptionPricer(valuationDate,marketInstruments,modelName);
hestonEuropeanOptionPricer.calculatePrice;
hestonPricelist = hestonEuropeanOptionPricer.getInstrumentList;
modelName = 'BSM';
optionPricer = EuropeanOptionPricer(valuationDate,hestonPricelist,modelName);
optionPricer.calculateImpliedVol;
hestonPricelist = optionPricer.getInstrumentList;

for i = 1:length(marketInstruments)
    instrument = marketInstruments(i);
    hestonPrice(i,1) = instrument.getPrice;
    hestonVol(i,1) = instrument.getImpliedVol;
    instrument.setPrice(marketQuote(i,1));
    instrument.setImpliedVol(marketVol(i,1));
end


%% Bates
modelName = 'Bates';
batesEuropeanOptionPricer = EuropeanOptionPricer(valuationDate,marketInstruments,modelName);
batesEuropeanOptionPricer.calculatePrice;
batesPricelist = batesEuropeanOptionPricer.getInstrumentList;
modelName = 'BSM';
optionPricer = EuropeanOptionPricer(valuationDate,batesPricelist,modelName);
optionPricer.calculateImpliedVol;
batesPricelist = optionPricer.getInstrumentList;

for i = 1:length(marketInstruments)
    instrument = marketInstruments(i);
    batesPrice(i,1) = instrument.getPrice;
    batesVol(i,1) = instrument.getImpliedVol;
end

priceTable = table(id,strike,maturity,marketQuote,blackPrice,hestonPrice,batesPrice);
volTable = table(id,strike,maturity,marketVol,blackVol,hestonVol,batesVol);


futurePriceLoc = [folderLocation,'\FuturesPrices.txt'];
futsPriceTable = readtable(futurePriceLoc,"FileType","text","Delimiter",",");


% plotModelPrices(priceTable,futsPriceTable)
hestonParams = hestonEuropeanOptionPricer.getModel.getCalibratedParam;
batesParams = batesEuropeanOptionPricer.getModel.getCalibratedParam;

plotModelIVs(volTable,coinName,valuationDate,hestonParams,batesParams)





