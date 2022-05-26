classdef TimeSeries < handle
    % @dev: class to build standardised timeseries
    %% Private Properties
    properties (Access = private)
        name; % dev: char
        dates; % dev: yyyymmdd
        values; % dev: double vec
        type; % dev: TimeseriesType.ENUM
    end

    %% Constructor
    methods
        function this = TimeSeries(name,dates,values,type)
            this.name = name;
            if ~isenum(type) % dev: check input type is timeseries ENUM
                Logger.getInstance().log(LogType.WARN,'Parsed timeseries type not an ENUM of TimeseriesType, set to UN_ASSIGNED');
                this.type = TimeseriesType.UN_ASSIGNED;
            else
                this.type = type;
            end

            if (length(dates) == length(values)) % dev: check input lengths are identical
                this.dates = dates;
                this.values = values;
                Logger.getInstance.log(LogType.INFO,['Input parsed and ',char(this.type), ' timeseries built for: ',name]);
            else
                Logger.getInstance.log(LogType.FATAL,'Input dates and values do no match in length, cannot build timeseries');
            end
        end

    end

    %% Public methods
    methods (Access = public)

        % getReturnTimeSeries
        function ReturnTSObj = getReturnTimeseries(this, method,varargin)
            % dev: check price timeseries is used to build return timeseries
            if ~isPrice(this)
                Logger.getInstance.log(LogType.FATAL,'Cannot build RETURN type timeseries from non-price timeseries');
            end
            % check if multiplicative factor is parsed
            if (~isempty(varargin))
                var1 = varargin{1};
                var2 = varargin{2};
                if (strcmp(var1,'annulisationFactor') && var2)
                    multFact = var2;
                end
            else
                multFact = 250; % dev: default value
            end
            ReturnTSObj = buildReturnTimeseries(this,method,multFact);
        end

        %getTimeseriesBetweenDates
        function TSObj = getTimeseriesBetweenDates(this,startDate,endDate)

            % dev: get dates [T_0, T_n] and
            % check for start > T_0 and end < T_n
            if (datenum(this.dates(1)) > datenum(str2double(startDate)) || ...
                    datenum(this.dates(end)) < datenum(str2double(endDate)))
                Logger.getInstance.log(LogType.FATAL,'Insufficient length of stored timeseries to extract per dates');
                TSObj = [];
            else

                idStart = find(this.dates == str2double(startDate), 1, 'first');
                idEnd = find(this.dates == str2double(endDate), 1, 'first');
                dates_n = this.dates(idStart:idEnd,:);
                values_n = this.values(idStart:idEnd,:);
                name_n = [this.name,'_', char(this.type),'_from:_', startDate, '_to:_', endDate];
                TSObj = TimeSeries(name_n,dates_n,values_n,this.type);
            end
        end

        % getLastValue
        function val = getLastValue(this)
            val = this.values;
            val = val(end);
        end

        % getLastDate
        function dat = getLastDate(this)
            dat = this.dates;
            dat = dat(end);
        end
    end

    %% Private methods
    methods (Access = private)
        % buildReturnTimeSeries
        function ReturnTSObj = buildReturnTimeseries(this,method,multFact)

            if (strcmp(method,'periodic') || strcmp(method,'continuous'))
                % dev: check supported price2ret methods
                prices = this.values;
                % dev: calculate returns
                values_r = price2ret(prices,'method',method) * sqrt(multFact);
                dates_r = this.dates(1:end-1);
                name_r = [this.name,'_', char(TimeseriesType.RETURN)];
                type_r = TimeseriesType.RETURN;
                % dev: create new TimeSeries instance and return
                ReturnTSObj = TimeSeries(name_r,dates_r,values_r,type_r);
            else
                Logger.getInstance.log(LogType.FATAL,['Return calculation method: ', method, ' not supported']);
                ReturnTSObj = [];
            end
        end
    end

    %% Getters and Setters
    methods

        %getName
        function name = getName(this)
            name = this.name;
        end

        %setName
        function setName(this,name)
            this.name = name;
        end

        %getDates
        function dates = getDates(this)
            dates = this.dates;
        end

        %setDates
        function setDates(this,dates)
            this.dates = dates;
        end
        %getValues
        function values = getValues(this)
            values = this.values;
        end

        %setValues
        function setValues(this,values)
            this.values = values;
        end

        %getType
        function type = getType(this)
            type = this.type;
        end

        %setType
        function setType(this,type)
            this.type = type;
        end
    end
end

