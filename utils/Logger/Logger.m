classdef (Sealed) Logger < handle
    % @dev:Singlenton logger class

    %% Restricted methods
    methods (Access = private)
        function obj = Logger
        end
    end

    %% Properies
    properties (Access = private)
        msg = struct('fileName',[],...
            'funcName',[],...
            'lineNum',[],...
            'logType',[],...
            'logString',[]);
        msgList = char([]);
    end

    %% Constructor
    methods (Static)
        function singleObj = getInstance
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                localObj = Logger;
            end
            singleObj = localObj;
        end
    end

    %% Methods
    methods (Access = public)
        %log
        function obj = log(obj,logType,logString)
            msgObj = createMsg(obj,logType,logString);
            msgString = printLog(obj,msgObj);
            obj.msgList{end+1,1} = msgString;

            % @dev: check if last log is FATAL
            % if yes, kill!!!
            if isFatal(msgObj)
                error('FATAL err encountered, process terminated')
            end
        end

        %saveOutput
        function saveOutput(obj,outputPath,varargin)
            if (isstring(outputPath))
                outputPath = char(outputPath);
            end

            if (~isempty(varargin))
                var1 = varargin{1};
                var2 = varargin{2};
                if (strcmp(var1,'timestamp') && var2)
                outputPath = strsplit(outputPath,'.');
                outputPath = [outputPath{1},' ',datestr(now,'HHMMSS'),'.',outputPath{2}];
                end
            end

            writecell(obj.msgList,outputPath);
        end
    end

    %% Methods
    methods (Access = private)
        %printLog
        function msgString = printLog(obj,msgObj)
            msgString = [datestr(now,'HH:MM:SS'),...
                ' ',char(msgObj.logType),...
                ' in file: ',msgObj.fileName,...
                ' in function: ',msgObj.funcName,...
                ' at line: ',num2str(msgObj.lineNum),...
                ': ',msgObj.logString];
            disp(msgString)
        end

        %createMsg
        function msgObj = createMsg(obj,logType,logString)
            if (~isenum(logType))
                logType = LogType.NONE;
            end

            if (isstring(logString))
                logString = char(logString);
            end

            stack = dbstack();
            lastEntry = height(stack);
            if (lastEntry == 0)
                msgObj.fileName = 'None';
                msgObj.funcName = 'None';
                msgObj.lineNum = 'None';
                msgObj.logType = logType;
                msgObj.logString = logString;
            else 
                ST = stack(3);
                msgObj.fileName = ST.file;
                msgObj.funcName = ST.name;
                msgObj.lineNum = ST.line;
                msgObj.logType = logType;
                msgObj.logString = logString;
            end
        end
    end

    %% Getters and Setters
    methods
        %getMsg
        function msgList = getMsgList(obj)
            msgList = obj.msgList;
        end
    end
end
