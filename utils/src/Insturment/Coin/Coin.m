classdef Coin < Instrument
    % @dev: abstract class to parent/ child instruments modelled

    %% Private Properties
    properties (Access = private)
        type; % dev: token or coin-> CoinType ENUM
        baseCurreny;% Currency ENUM
        priceTS; % Timeseries
    end

    %% Constructor
    methods
        function this = Coin(varargin)
            this = this@Instrument(varargin{:});
             Logger.getInstance.log(LogType.INFO,'Coin initalised');
        end

    end
    %% Getters and Setters
    methods
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
        function baseCurreny = getBaseCurreny(this)
            baseCurreny = this.baseCurreny;
        end

        %setBaseCurreny
        function setBaseCurreny(this,baseCurreny)
            % dev: check if parsed type isa CoinType ENUM
            if isa(baseCurreny,'Currency')
                this.baseCurreny = baseCurreny;
                Logger.getInstance.log(LogType.INFO,...
                    'BaseCurrency type provided and set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Parse type not of Currency ENUM');
            end
        end

    end
end

