classdef EuropeanOption < Option
    % @dev: abstract class to parent/ child instruments modelled

    %% Private Properties
    properties (Access = private)
        impliedVol;
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

        %getImpliedVol
        function impliedVol = getImpliedVol(this)
            impliedVol = this.impliedVol;
        end

        %setImpliedVol
        function setImpliedVol(this,impliedVol)
            this.impliedVol = impliedVol;
        end

    end
end

