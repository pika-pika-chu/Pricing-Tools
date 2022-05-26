classdef Option < Instrument
    % @dev: abstract class to parent/ child instruments modelled

    %% Private Properties
    properties (Access = private)
        type;
        strike;
        maturity;
        exerciseStyle; % dev: set in child class constructor
        price; % dev: getters and setters in child class
        underlying;
    end

    %% Constructor
    methods
        function this = Option(varargin)
            this = this@Instrument(varargin{:});
             Logger.getInstance.log(LogType.INFO,'Option initalised');
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
            if isempty(type)
                Logger.getInstance.log(LogType.WARN,...
                    'Type not parsed to EuropeanOption, set to default = call');
                this.type = 'call';
            else
                this.type = type;
            end
        end

        %getStrike
        function strike = getStrike(this)
            strike = this.strike;
        end

        %setStrike
        function setStrike(this,strike)
            if isempty(strike)
                Logger.getInstance.log(LogType.FATAL,...
                    'Strike not parsed to EuropeanOption, Terminating...');
            else
                this.strike = strike;
            end
        end

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

        %getExerciseStyle
        function exerciseStyle = getExerciseStyle(this)
            exerciseStyle = this.exerciseStyle;
        end

        %setExerciseStyle
        function setExerciseStyle(this,exerciseStyle)
            this.exerciseStyle = exerciseStyle;
        end

                %getPrice
        function price = getPrice(this)
            price = this.price;
        end

        %setPrice
        function setPrice(this,price)
            this.price = price;
        end

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

