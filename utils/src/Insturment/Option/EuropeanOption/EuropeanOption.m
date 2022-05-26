classdef EuropeanOption < Option
    % @dev: abstract class to parent/ child instruments modelled

    %% Private Properties
    properties (Access = private)
    end

    %% Constructor
    methods
        function this = EuropeanOption(varargin)
            this = this@Option(varargin{:});
            this.setExerciseStyle(ExerciseStyle.EUROPEAN);
            Logger.getInstance.log(LogType.INFO,'EuropeanOption initalised');
        end

    end
    %% Getters and Setters
    methods (Access = public)


    end
end

