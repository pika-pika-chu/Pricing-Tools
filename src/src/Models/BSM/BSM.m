classdef BSM < handle
    % @dev: class for Black Scholes Metron

    %% Private Properties
    properties (Access = private)
        underlying;

    end

    %% Constructor
    methods
        function this = BSM(coin)
            this.underlying = coin;
        end


    end


    %% Public methods
    methods (Access = public)

        % getCallPrice
        function price = getCallPrice(this,strike,timeToMaturity,interestRate,volatility,varargin)
            
            if ~isempty(varargin)   
            [sInitial,div] = getVararginValues(this,varargin{:});
            
            else
            sInitial = this.underlying.getPriceTS.getLastValue;
            div = 0;
            
            end

            price = calculateCallPrice(this,strike,timeToMaturity,interestRate,volatility,sInitial,div);
        end

        % getPutPrice
        function price = getPutPrice(this,strike,timeToMaturity,interestRate,volatility,varargin)

            if ~isempty(varargin)               
            [sInitial,div] = getVararginValues(varargin{:});

            else
            sInitial = this.underlying.getPriceTS.getLastValue;
            div = 0;

            end

            price = calculatePutPrice(this,strike,timeToMaturity,interestRate,volatility,sInitial,div);
        end
        
        % getImpliedVol
        function impVol = getImpliedVol(this,price, strike,timeToMaturity,interestRate, type, varargin)

             if ~isempty(varargin)   
            [sInitial,div] = getVararginValues(this,varargin{:});
            
            else
            sInitial = this.underlying.getPriceTS.getLastValue;
            div = 0;
            
            end

            if strcmp(type,'call')
                impVol = calculateCallImpliedVol(this,price,sInitial,strike,timeToMaturity,interestRate,div);

            elseif strcmp(type,'put')
                impVol = calculatePutImpliedVol(this,price,sInitial,strike,timeToMaturity,interestRate,div);

            else
                Logger.getInstance.log(LogType.FATAL,...
                    ['Type unknown for id = ', num2str(option.getId),...
                    ' name = ',option.getName]);

            end
        end
    end

    %% Private methods
    methods (Access = private)

        % calculateCallPrice
        function callPrice = calculateCallPrice(this,strike,timeToMaturity,interestRate,volatility,sInitial,div)
            
            d1 = (log(sInitial/strike) + ((interestRate - div) + 0.5*volatility.^2).*timeToMaturity ) ./ ...
                (volatility * sqrt(timeToMaturity));
            d2 = d1 - volatility * sqrt(timeToMaturity);

            callPrice = sInitial * normcdf(d1) - strike *exp(-(interestRate - div) * (timeToMaturity))*normcdf(d2);
        end

        % calculatePutPrice
        function callPrice = calculatePutPrice(this,strike,timeToMaturity,interestRate,volatility,sInitial,div)

            d1 = (log(sInitial/strike) + ((interestRate - div) + 0.5*volatility.^2).*(timeToMaturity) ) ./ ...
                (sigma.*sqrt(timeToMaturity));
            d2 = d1 - volatility * sqrt(timeToMaturity);

            callPrice = normcdf(-d2) * strike * exp(-(interestRate - div)*(timeToMaturity))-normcdf(-d1) * sInitial;
        end

        %getVararginValues
        function [sInitial,div] = getVararginValues(this,varargin)
            % set empty values before checking
            sInitial = [];
            div = [];
            % check if name value arg is present in varargin
            sInitialIdx = find(strcmp(varargin,'sInitial'),1);
            divIdx = find(strcmp(varargin,'dividend'),1);
               % set name value args if present in varargin
               if ~isempty(sInitialIdx)
                    sInitial = varargin{sInitialIdx + 1};  
               end
            
               if ~isempty(divIdx)
                    div = varargin{divIdx + 1};         
               end

        end
        
        % calculateCallImpliedVol
        function impliedVol = calculateCallImpliedVol(this,price,sInitial,strike,timeToMaturity,interestRate,div)

            % herf: https://www.researchgate.net/publication/245065192_A_Simple_Formula_to_Compute_the_Implied_Standard_Deviation

            fun = @(x, price, sInitial, strike, timeToMaturity, interestRate, div) ...
                price - calculateCallPrice(this, strike, timeToMaturity, interestRate, x, sInitial,div);

            initialVol = sqrt(2 * pi / timeToMaturity) * (price / sInitial);
            impliedVol = max(fzero( @(x)fun(x,price, sInitial, strike, timeToMaturity, interestRate, div), initialVol),0);
        
        end

        % calculatePutImpliedVol
        function impliedVol = calculatePutImpliedVol(this,price,sInitial,strike,timeToMaturity,interestRate,div)

            % herf: https://www.researchgate.net/publication/245065192_A_Simple_Formula_to_Compute_the_Implied_Standard_Deviation

            fun = @(x, price, sInitial, strike, timeToMaturity, interestRate, div) ...
                price - calculatePutPrice(this, strike, timeToMaturity, interestRate, x, sInitial,div);

            initialVol = sqrt(2 * pi / timeToMaturity) * (price / sInitial);
            impliedVol = max(fzero( @(x)fun(x,price, sInitial, strike, timeToMaturity, interestRate, div), initialVol),0);
        
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
end

