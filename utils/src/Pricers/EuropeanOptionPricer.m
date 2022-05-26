classdef EuropeanOptionPricer < handle
    % @dev: abstract class to implement pricing routine for
    % european calls/ puts

    %% Private Properties
    properties (Access = private)
        model;
        instrumentList;
        underlying;
    end

    %% Constructor
    methods
        function this = EuropeanOptionPricer(underlying,instrumentList,model)
            if (isa(underlying,'Coin'))
                this.underlying = underlying;
                Logger.getInstance.log(LogType.INFO,...
                    'Underlying instruments parsed to EuropeanOptionPricer');

            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Incorrect underlying instruments parsed to EuropeanOptionPricer');
            end

            for i = 1:length(instrumentList)
                instrument = instrumentList(i);
                if (isa(instrument,'EuropeanOption'))
                    Logger.getInstance.log(LogType.INFO,...
                        ['Instrument with Id. ',num2str(instrument.getId), ' is valid for EuropeanOptionPricer']);

                else
                    Logger.getInstance.log(LogType.FATAL,...
                        ['Incorrect instrument parsed for Id. ',num2str(instrument.getId), ' to EuropeanOptionPricer']);
                end
            end
            this.instrumentList = instrumentList;

            Logger.getInstance.log(LogType.INFO,...
                'Model parsed to EuropeanOptionPricer, trying to calibrate');

            %set model
            setModel(this,model);

        end

    end

    %% Public methods
    methods (Access = public)

        % calculatePrice
        function calculatePrice(this)
            % dev: apply payoff and update price in the instrument list
            for i = 1:length(this.instrumentList)
                instrument = this.instrumentList(i);

                % dev: get instrument and log
                Logger.getInstance.log(LogType.INFO,...
                    ['Pricing instruement with id = ', num2str(instrument.getId),...
                    ' name = ',instrument.getName,' and type = ',instrument.getType]);

                price = computePrice(this,instrument);
                this.instrumentList(i).setPrice(price);
            end
        end
    end

    %% Private methods
    methods (Access = private)

        %ApplyCallPayoff
        function price = computePrice(this,option)

            [sInitial,t0] = getInitialPriceAndT0(this);
            % get necessary properties
            strike = option.getStrike;

            % dev: get the number of periods that match maturity date
            maturityInPeriods = abs(days365(datenum(t0),datenum(option.getMaturity)));
            % dev: check if returns are simulated and stored to model
            if isempty(this.model.getSimulatedPrices)
                this.model.simulate;
            end
            % dev: get the simulated returns
            simulatedPrices = this.model.getSimulatedPrices;
            % dev: check if simulated returns are longer than or equal to
            % maturity in periods + 1
            % since, P(:,1) = sinitial
            if width(simulatedPrices) < maturityInPeriods + 1
                Logger.getInstance.log(LogType.INFO,...
                    ['Not enough prices simulated for instrument with id = ', num2str(option.getId),...
                    ' name = ',option.getName,' and type = ',option.getType]);
                Logger.getInstance.log(LogType.INFO,...
                    ['simulated return count = ', num2str(width(simulatedPrices)),...
                    ' required = ',num2str(maturityInPeriods),' for maturity = ',option.getMaturity]);

            else

                if strcmp(option.getType,'call')
                    price = applyCallPayoff(this,simulatedPrices,maturityInPeriods,strike);
                    Logger.getInstance.log(LogType.INFO,...
                        ['Price calculated and stored for instrument with id = ', num2str(option.getId),...
                        ' name = ',option.getName,' and type = ',option.getType]);

                elseif strcmp(option.getType,'put')
                    price = applyPutPayoff(this,simulatedPrices,maturityInPeriods,strike);
                    Logger.getInstance.log(LogType.INFO,...
                        ['Price calculated and stored for instrument with id = ', num2str(option.getId),...
                        ' name = ',option.getName,' and type = ',option.getType]);

                else
                    Logger.getInstance.log(LogType.FATAL,...
                        ['Type unknown for id = ', num2str(option.getId),...
                        ' name = ',option.getName]);
                end
            end
        end

        %applyPutPayoff
        function price = applyCallPayoff(this,simulatedPrices,maturityInPeriods,strike)
            terminalPrices = simulatedPrices(:,maturityInPeriods+1);
            % max(s-k,0)
            terminalPayoffs = max(terminalPrices - strike,0);
            % dev: TEMP DISCOUNTING
            price = mean(terminalPayoffs * exp(-0.03 * maturityInPeriods/365));
        end

        %applyPutPayoff
        function price = applyPutPayoff(this,simulatedPrices,maturityInPeriods,strike)

            terminalPrices = simulatedPrices(:,maturityInPeriods+1);
            % max(k-s,0)
            terminalPayoffs = max(strike -terminalPrices ,0);

            price = mean(terminalPayoffs);
        end

        %getInitialPriceAndT0
        function [sInitial,t0] = getInitialPriceAndT0(this)
            sInitial = this.underlying.getPriceTS.getLastValue;
            t0 = datestr(this.underlying.getPriceTS.getLastDate);
        end

        %setModel
        function setModel(this,model)
            if (strcmp(model,'SVCJ')) % dev: check model name

                if isempty(this.getModel)  % dev: if empty create and calibrate
                    mod = SVCJ(this.underlying);
                    mod.calibrate();
                    if (mod.isCalibrated == ModelCalibration.SUCCESS)       
                        Logger.getInstance.log(LogType.INFO,...
                            'Model is calibrated and set to EuropeanOptionPricer');
                    else
                        Logger.getInstance.log(LogType.WARN,...
                            'Potentially non-clibrated model set to EuropeanOptionPricer');
                    end
                    this.model = mod;
                elseif (isa(this.model,'SVCJ') && this.model.isCalibrated == ModelCalibration.NOT_CALIBRATED)
                    this.model.calibrate();
                    if (this.model.isCalibrated == ModelCalibration.SUCCESS)       
                        Logger.getInstance.log(LogType.INFO,...
                            'Model is calibrated and set to EuropeanOptionPricer');
                    else
                        Logger.getInstance.log(LogType.WARN,...
                            'Potentially non-clibrated model set to EuropeanOptionPricer');
                    end

                else
                    Logger.getInstance.log(LogType.INFO,...
                            ['SVCJ model already set to EuropeanOptionPricer with calibration status: ', char(this.model.isCalibrated)]);
                end

            else
                Logger.getInstance.log(LogType.FATAL,...
                    [model, ' not implemented']);

            end



        end
    end

    %% Getters and Setters
    methods (Access = public)

        %getModel
        function model = getModel(this)
            model = this.model;
        end

        %getInstrumentList
        function instrumentList = getInstrumentList(this)
            instrumentList = this.instrumentList;
        end

        %setInstrumentList
        function setInstrumentList(this,instrumentList)
            for i = 1:length(instrumentList)
                instrument = instrumentList(i);
                if (isa(instrument,'EuropeanOption'))
                    Logger.getInstance.log(LogType.INFO,...
                        ['Instrument with Id. ',num2str(instrument.getId), ' is valid for EuropeanOptionPricer']);

                else
                    Logger.getInstance.log(LogType.FATAL,...
                        ['Incorrect instrument parsed for Id. ',num2str(instrument.getId), ' to EuropeanOptionPricer']);
                end
            end
            this.instrumentList = instrumentList;
        end
    end
end

