classdef Exotic < Option
    % @dev: abstract class to parent/ child instruments modelled

    %% Private Properties
    properties (Access = private)
     
    end

    %% Constructor
    methods
        function this = Exotic(varargin)
            this = this@Option(varargin{:});
        end

    end
    %% Getters and Setters
    methods
    end
end


