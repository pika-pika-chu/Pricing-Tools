classdef OptionPricer < handle

    %% Protected properties
    properties (Access = protected)
        underlying;
        iRCurve;
    end

    %% Constructor
    methods
        function this = OptionPricer(underlying)
            if (isa(underlying,'Coin'))
                this.underlying = underlying;
                Logger.getInstance.log(LogType.INFO,...
                    'Underlying instruments parsed to EuropeanOptionPricer');

            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Incorrect underlying instruments parsed to EuropeanOptionPricer');
            end

            settle = this.underlying.getPriceTS.getLastDate;
            currency = this.underlying.getBaseCurreny;

            this.iRCurve = sourceIRCurve(this,settle,currency);
        end

    end

    %% Public methods
    methods (Access = public)

        %getUnderlying
        function underlying = getUnderlying(this)
            underlying = this.underlying;
        end

        %setUnderlying
        function setUnderlying(this,underlying)
            this.underlying = underlying;
        end
    end

    %% Private methods
    methods (Access = private)

        % sourceIRCurve
        function irCurve = sourceIRCurve(this, settle,currency)
                irCurve = buildIRCurve(datestr(settle),currency);
        end

    end
end

