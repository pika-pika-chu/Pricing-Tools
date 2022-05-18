function fatal = isFatal(msgObj)
fatal = false;
if (msgObj.logType == LogType.FATAL)
    fatal = true;
end
end

