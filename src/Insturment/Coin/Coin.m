classdef Coin < Instrument
    % @dev: abstract class to parent/ child instruments modelled

    %% Private Properties
    properties (Access = private)
        type; % dev: token or coin-> CoinType ENUM
        baseCurrency;% Currency ENUM
        priceTS; % Timeseries
        divYield;
    end

    %% Constructor
    methods
        function this = Coin(varargin)
            this = this@Instrument(varargin{:});
            Logger.getInstance.log(LogType.INFO,'Coin initalised');
        end

    end
    %% Getters and Setters
    methods (Access = public)
        %getType
        function type = getType(this)
            type = this.type;
        end

        %setType
        function setType(this,type)
            % dev: check if parsed type isa CoinType ENUM
            if isa(type,'CoinType')
                this.type = type;
                Logger.getInstance.log(LogType.INFO,...
                    'Coin type provided and set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Parse type not of CoinType ENUM');
            end
        end
        %getPriceTS
        function priceTS = getPriceTS(this)
            priceTS = this.priceTS;
        end

        %setPriceTS
        function setPriceTS(this,priceTS)
            % dev: check if parsed time-series isa TimeSeries type
            % check type = TimeseriesType.PRICE
            if (isa(priceTS,'TimeSeries') && ...
                    priceTS.getType == TimeseriesType.PRICE)
                this.priceTS = priceTS;
                Logger.getInstance.log(LogType.INFO,...
                    'Price timeseries provided and set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Parse timeseries not of Timeseries.PRICE type');
            end
        end

        %getBaseCurreny
        function baseCurrency = getBaseCurrency(this)
            baseCurrency = this.baseCurrency;
        end

        %setBaseCurreny
        function setBaseCurrency(this,baseCurrency)
            % dev: check if parsed type isa CoinType ENUM
            if isa(baseCurrency,'Currency')
                this.baseCurrency = baseCurrency;
                Logger.getInstance.log(LogType.INFO,...
                    'BaseCurrency type provided and set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Parse type not of Currency ENUM');
            end
        end

        %getDivYield
        function divYield = getDivYield(this)
            divYield = this.divYield;
        end

        %setDivYield
        function setDivYield(this,divYield)
            this.divYield = divYield;
        end

    end
end

