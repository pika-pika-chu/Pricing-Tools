classdef BSM < handle
    % @dev: class for Black Scholes Metron

    %% Private Properties
    properties (Access = private)
        underlying;

    end

    %% Constructor
    methods
        function this = BSM()
        end


    end


    %% Public methods
    methods (Access = public)

        % getCallPrice
        function price = getCallPrice(this,underlying,strike,timeToMaturity,interestRate,volatility,sInitial,div)

            % check the underlying and call correct pricing function
            if isa(underlying,'Coin')
                price = calculateBSCallPrice(this,strike,timeToMaturity,interestRate,volatility,sInitial,div);

            elseif isa(underlying,'Future')
                price = calculateBlackCallPrice(this,strike,timeToMaturity,interestRate,volatility,sInitial);
            end

        end

        % getPutPrice
        function price = getPutPrice(this,underlying,strike,timeToMaturity,interestRate,volatility,sInitial,div)

            % check the underlying and call correct pricing function
            if isa(underlying,'Coin')
                price = calculateBSPutPrice(this,underlyingstrike,timeToMaturity,interestRate,volatility,sInitial,div);

            elseif isa(underlying,'Future')
                price = calculateBlackPutPrice(this,strike,timeToMaturity,interestRate,volatility,sInitial);
            end

        end

        % getImpliedVol
        function impVol = getImpliedVol(this,underlying,price, strike,timeToMaturity,interestRate, type, sInitial,div)

            if (ischar(type) && strcmp(type,'call')) || iscell(type) && (sum(strcmp(type, 'call')) == length(type))

                % check the underlying and call correct pricing function
                if isa(underlying,'Coin')
                    impVol = calculateBSCallImpliedVol(this,price,sInitial,strike,timeToMaturity,interestRate,div);

                elseif isa(underlying,'Future')
                    impVol = calculateBlackImpVol(this,price,sInitial,strike,timeToMaturity,interestRate,'call');

                end

            elseif (ischar(type) && strcmp(type,'put')) || (sum(strcmp(type, 'put')) == length(type))

                % check the underlying and call correct pricing function
                if isa(this.underlying,'Coin')
                    impVol = calculateBSPutImpliedVol(this,price,sInitial,strike,timeToMaturity,interestRate,div);

                elseif isa(this.underlying,'Future')
                    impVol = calculateBlackImpVol(this,price,sInitial,strike,timeToMaturity,interestRate,'put');

                end

            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Calls and Put option type parsed at the same time, only one can be handeled');

            end
        end
    end

    %% Private methods
    methods (Access = private)

        % calculateBSCallPrice
        function callPrice = calculateBSCallPrice(this,strike,timeToMaturity,interestRate,volatility,sInitial,div)

            d1 = (log(sInitial./strike) + ((interestRate - div) + 0.5*volatility.^2).*timeToMaturity ) ./ ...
                (volatility .* sqrt(timeToMaturity));
            d2 = d1 - volatility .* sqrt(timeToMaturity);

            callPrice = sInitial .* normcdf(d1) - strike .* exp(-(interestRate - div) .* (timeToMaturity)) .* normcdf(d2);
        end

        % calculateBSPutPrice
        function putPrice = calculateBSPutPrice(this,strike,timeToMaturity,interestRate,volatility,sInitial,div)

            d1 = (log(sInitial./strike) + ((interestRate - div) + 0.5*volatility.^2).*(timeToMaturity) ) ./ ...
                (sigma.*sqrt(timeToMaturity));
            d2 = d1 - volatility .* sqrt(timeToMaturity);

            putPrice = normcdf(-d2) .* strike .* exp(-(interestRate - div).*(timeToMaturity))-normcdf(-d1) .* sInitial;
        end

        % calculateBlackCallPrice
        function callPrice = calculateBlackCallPrice(this,strike,timeToMaturity,interestRate,volatility,sInitial)

            callPrice = blkprice(sInitial,strike,interestRate,timeToMaturity,volatility);
        end

        % calculateBlackPutPrice
        function putPrice = calculateBlackPutPrice(this,strike,timeToMaturity,interestRate,volatility,sInitial)

            [~,putPrice] = blkprice(sInitial,strike,interestRate,timeToMaturity,volatility);
        end

        % calculateBlackImpVol
        function impliedVol = calculateBlackImpVol(this,price,sInitial,strike,timeToMaturity,interestRate,type)

            if strcmp(type,'call')
                impliedVol = blkimpv(sInitial,strike,interestRate,timeToMaturity,price,'Class', 'call');

            elseif strcmp(type,'put')
                impliedVol = blkimpv(sInitial,strike,interestRate,timeToMaturity,price,'Class', 'put');

            end
        end


        % calculateBSCallImpliedVol
        function impliedVol = calculateBSCallImpliedVol(this,price,sInitial,strike,timeToMaturity,interestRate,div)
            impliedVol = nan(length(strike),1);
            % herf: https://www.researchgate.net/publication/245065192_A_Simple_Formula_to_Compute_the_Implied_Standard_Deviation

            fun = @(x, price, sInitial, strike, timeToMaturity, interestRate, div) ...
                price - calculateBSCallPrice(this, strike, timeToMaturity, interestRate, x, sInitial,div);

            for i = 1:length(strike)
                initialVol = sqrt(2 * pi / timeToMaturity(i)) * (price(i) / sInitial(i));
                impliedVol(i) = max(fzero( @(x) fun(x,price(i), sInitial(i),...
                    strike(i), timeToMaturity(i), interestRate(i), div(i)), initialVol),0);
            end
        end

        % calculateBSPutImpliedVol
        function impliedVol = calculateBSPutImpliedVol(this,price,sInitial,strike,timeToMaturity,interestRate,div)
            impliedVol = nan(length(strike),1);
            % herf: https://www.researchgate.net/publication/245065192_A_Simple_Formula_to_Compute_the_Implied_Standard_Deviation

            fun = @(x, price, sInitial, strike, timeToMaturity, interestRate, div) ...
                price - calculateBSPutPrice(this, strike, timeToMaturity, interestRate, x, sInitial,div);

            for i = 1:length(strike)
                initialVol = sqrt(2 * pi / timeToMaturity(i)) * (price(i) / sInitial(i));
                impliedVol(i) = max(fzero( @(x) fun(x,price(i), sInitial(i),...
                    strike(i), timeToMaturity(i), interestRate(i), div(i)), initialVol),0);
            end
        end
    end
end

