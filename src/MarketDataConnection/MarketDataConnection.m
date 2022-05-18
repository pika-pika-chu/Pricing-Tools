classdef MarketDataConnection < handle
    % @dev: Purpose of this call is to isolate all connection to market
    % data (internal or external)
    % @dev: this is to ensure that code has uniform inputs
    
    %% Private Properties
    properties (Access = private)
        conn; %
    end
    
    %% COnstructor
    methods
        function this = MarketDataConnection()
        end
       
    end
end

