classdef Option < Instrument
    % @dev: abstract class to parent/ child instruments modelled

    %% Private Properties
    properties (Access = private)
        type;
        strike;
        maturity;
        exerciseStyle; % dev: set in child class constructor
        price; % dev: getters and setters in child class
    end

    %% Constructor
    methods
        function this = Option(varargin)
            this = this@Instrument(varargin);
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
            this.type = type;
        end

        %getStrike
        function strike = getStrike(this)
            strike = this.strike;
        end

        %setStrike
        function setStrike(this,strike)
            this.strike = strike;
        end

        %getMaturity
        function maturity = getMaturity(this)
            maturity = this.maturity;
        end

        %setMaturity
        function setMaturity(this,maturity)
            this.maturity = maturity;
        end

        %getExerciseStyle
        function exerciseStyle = getExerciseStyle(this)
            exerciseStyle = this.exerciseStyle;
        end

        %setExerciseStyle
        function setExerciseStyle(this,exerciseStyle)
            this.exerciseStyle = exerciseStyle;
        end
    end
end

