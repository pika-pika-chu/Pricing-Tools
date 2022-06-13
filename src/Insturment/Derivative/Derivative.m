classdef Derivative < Instrument
    % @dev: abstract class to parent/ child instruments modelled

    %% Private Properties
    properties (Access = private)
        maturity;
        underlying;
        price;
        valuationDate;
        baseCurrency;% Currency ENUM
    end

    %% Constructor
    methods
        function this = Derivative(varargin)
            this = this@Instrument(varargin{:});
        end

    end
    %% Getters and Setters
    methods

        %getMaturity
        function maturity = getMaturity(this)
            maturity = this.maturity;
        end

        %setMaturity
        function setMaturity(this,maturity)
            if isempty(maturity)
                Logger.getInstance.log(LogType.FATAL,...
                    'Maturity not parsed to EuropeanOption, Terminating');
            else
                this.maturity = maturity;
            end
        end

        %getUnderlying
        function underlying = getUnderlying(this)
            underlying = this.underlying;
        end

        %setUnderlying
        function setUnderlying(this,underlying)
            this.underlying = underlying;
        end


        %getPrice
        function price = getPrice(this)
            price = this.price;
        end

        %setPrice
        function setPrice(this,price)
            this.price = price;
        end

        %getValuationDate
        function valuationDate = getValuationDate(this)
            valuationDate = this.valuationDate;
        end

        %setValuationDate
        function setValuationDate(this,valuationDate)
            this.valuationDate = valuationDate;
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

    end
end

