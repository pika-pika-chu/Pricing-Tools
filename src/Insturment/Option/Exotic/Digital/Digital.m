classdef Digital < Exotic
    % @dev: abstract class to parent/ child instruments modelled

    %% Private Properties
    properties (Access = private)

    end

    %% Constructor
    methods
        function this = Digital(varargin)
            this = this@Exotic(varargin{:});
        end

    end
    %% Getters and Setters
    methods (Access = public)

        %getPrice
        function price = getPrice(this)
            price = this.price;
        end

        %setPrice
        function setPrice(this,price)
            this.price = price;
        end
    end
end


