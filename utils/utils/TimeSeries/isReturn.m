function ret = isReturn(TSObj)
ret = false;
if (TSObj.getType == TimeseriesType.RETURN)
    ret = true;
end
end

