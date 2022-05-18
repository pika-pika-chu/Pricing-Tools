function price = isPrice(TSObj)
price = false;
if (TSObj.getType == TimeseriesType.PRICE)
    price = true;
end
end

