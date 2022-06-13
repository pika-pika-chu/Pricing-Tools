classdef Option < Derivative
    % @dev: abstract class to parent/ child instruments modelled

    %% Private Properties
    properties (Access = protected)
        type;
        strike;
        exerciseStyle; % dev: set in child class constructor

    end

    %% Constructor
    methods
        function this = Option(varargin)
            this = this@Derivative(varargin{:});
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

