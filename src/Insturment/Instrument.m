classdef Instrument < handle 
    % @dev: abstract class to parent instruments modelled

    %% Private Properties
    properties (Access = private)
        name; % dev: char
        id; % dev: char
        ticker; % dev: char
    end

    %% Constructor
    methods
        function this = Instrument(varargin)
        end
    end


    %% Getters and Setters
    methods
        %getId
        function id = getId(this)
            id = this.id;
        end
        %setId
        function setId(this,id)
            this.id = id;
        end

        %getName
        function name = getName(this)
            name = this.name;
        end
        %setName
        function setName(this,name)
            this.name = name;
        end
        %getTicker
        function ticker = getTicker(this)
            ticker = this.ticker;
        end
        %setTicker
        function setTicker(this,ticker)
            this.ticker = ticker;
        end
    end
end
