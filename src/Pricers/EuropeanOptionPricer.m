classdef EuropeanOptionPricer < handle
    % @dev: abstract class to implement pricing routine for
    % european calls/ puts

    %% Private Properties
    properties (Access = private)
        model;
        instrumentList;
    end

    %% Constructor
    methods
        function this = EuropeanOptionPricer(europeanOptions,model)
            if (isa(europeanOptions,'EuropeanOption'))
                this.instrumentList = europeanOptions;
                Logger.getInstance.log(LogType.INFO,...
                    'Instruments parsed to EuropeanOptionPricer');

            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Incorrect instruments parsed to EuropeanOptionPricer');
            end

            %set model
            setModel(this,model);
        end

    end

    %% Public methods
    methods (Access = public)

        % calculatePrice
        function calculatePrice(this)
            % dev: apply payoff and update price in the instrument list
        end
    end

    %% Private methods
    methods (Access = private)

        %ApplyCallPayoff
        function price = ApplyCallPayoff(this,europeanCall)
            price = 0;
        end

        %ApplyPutPayoff
        function price = ApplyPutPayoff(this,europeanPut)
           price = 0;
        end
    end

    %% Getters and Setters
    methods (Access = public)

        %getModel
        function model = getModel(this)
            model = this.model;
        end

        %setModel
        function setModel(this,model)

            if (model.isCalibrated == ModelCalibration.SUCCESS)
                this.model = model;
                Logger.getInstance.log(LogType.INFO,...
                    'Model set to EuropeanOptionPricer');
            else
                Logger.getInstance.log(LogType.WARN,...
                    'Potentially non-clibrated model set to EuropeanOptionPricer');
            end

        end
    end
end

