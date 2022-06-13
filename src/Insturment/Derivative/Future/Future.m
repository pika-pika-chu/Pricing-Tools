classdef Future < Derivative
    % @dev: abstract class to parent/ child instruments modelled

    %% Private Properties
    properties (Access = protected)

    end

    %% Constructor
    methods
        function this = Future(varargin)
            this = this@Derivative(varargin{:});
             Logger.getInstance.log(LogType.INFO,'Future initalised');
        end

    end
    %% Getters and Setters
    methods
        
    end
end

