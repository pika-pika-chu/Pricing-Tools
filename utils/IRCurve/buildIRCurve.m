function irCurve = buildIRCurve(settle,currency)
% dev: sources dates and type from .txt file to build IR curve
fileName = [settle,'-',char(currency),'.txt'];
fileLocation = ['utils\IRCurve\IRCurveData\',fileName];
tt = readtable(fileLocation,"FileType","text","Delimiter",",");

if isempty(tt)
Logger.getInstance.log(LogType.FATAL,...
                        ['File not found for IR curve. ', fileName]);

end

rateType = tt.type{1};

if strcmp(rateType,'zero')
rateType = IRCurveType.ZERO;

elseif strcmp(rateType,'forward')
rateType = IRCurveType.FORWARD;

elseif strcmp(rateType,'discount')
rateType = IRCurveType.DISCOUNT;
end

dates = tt.date;
values = tt.value;
ticker = [settle,'-',char(currency)];
irCurve = IRCurveObj(ticker,settle,currency,dates,values,rateType);

end

