classdef IRCurveObj < handle
    % dev: IRCurve object with all curve functionalities

    %% Private properties
    properties (Access = private)
        iRCurve;
        settle;
        currency;
        ticker;
    end

    %% Constructor
    methods
        function this = IRCurveObj(ticker,settle,currency,dates,iRValues,rateType, varargin)

            this.settle = settle;
            this.ticker = ticker;

            if isa(currency,'Currency')
                this.currency = currency;
                Logger.getInstance.log(LogType.INFO,...
                    ['IR Curve currency set to ', char(currency)]);
            else
                Logger.getInstance.log(LogType.WARN,...
                    'Currency parsed to IR Curve currency is not of Currency ENUM type');
            end

            if (length(dates) == length(iRValues))

                [compounding,basis,interpMethod] = getVararginValues(this,varargin);

                if isa(rateType,'IRCurveType')

                    if (rateType == IRCurveType.ZERO)
                        this.iRCurve = IRDataCurve('zero',settle,dates,iRValues,...
                            'Compounding',compounding,'basis',basis,'InterpMethod',interpMethod);
                        Logger.getInstance.log(LogType.INFO,...
                            ['Zero rate curve built for currency ', char(currency)]);

                    elseif (rateType == IRCurveType.FORWARD)
                        this.iRCurve = IRDataCurve('forward',settle,dates,iRValues,...
                            'Compounding',compounding,'basis',basis,'InterpMethod',interpMethod);
                        Logger.getInstance.log(LogType.INFO,...
                            ['Forward rate curve built for currency ', char(currency)]);

                    elseif (rateType == IRCurveType.DISCOUNT)
                        this.iRCurve = IRDataCurve('discount',settle,dates,iRValues,...
                            'Compounding',compounding,'basis',basis,'InterpMethod',interpMethod);
                        Logger.getInstance.log(LogType.INFO,...
                            ['Discount rate curve built for currency ', char(currency)]);

                    end
                else
                    Logger.getInstance.log(LogType.FATAL,...
                        'Incorrect IRCurveType ENUM parsed, cannot build IRCurve!');
                end

            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Number of dates and ir values parsed do not match, cannot build IRCurve!');
            end
        end
    end


    %% Public methods
    methods (Access = public)
        function rate = getZeroRates(this,date)
            rate = this.iRCurve.getZeroRates(date);
        end

        %getSettle
        function settle = getSettle(this)
            settle = this.settle;
        end

        %setSettle
        function setSettle(this,settle)
            this.settle = settle;
        end

        %getCurrency
        function currency = getCurrency(this)
            currency = this.currency;
        end

        %setCurrency
        function setCurrency(this,currency)
            this.currency = currency;
        end

    end
    %% Private methods
    methods (Access = private)
        % getVararginValues
        function [compounding,basis,interpMethod] = getVararginValues(this,varargin)

            % set default values before checking
            compounding = -1; % dev: matlab ir curve obj defualt
            basis = 0; % dev: matlab ir curve obj defualt
            interpMethod = 'linear'; % dev: matlab ir curve obj defualt

            % check if name value arg is present in varargin
            compoundingIdx = find(strcmp(varargin,'compounding'),1);
            basisIdx = find(strcmp(varargin,'basis'),1);
            interpMethodIdx = find(strcmp(varargin,'interpMethod'),1);

            % set name value args if present in varargin
            if ~isempty(compoundingIdx)
                compounding = varargin{compoundingIdx + 1};
            end

            if ~isempty(basisIdx)
                basis = varargin{basisIdx + 1};
            end

            if ~isempty(interpMethodIdx)
                interpMethod = varargin{interpMethodIdx + 1};
            end

        end
    end
end

