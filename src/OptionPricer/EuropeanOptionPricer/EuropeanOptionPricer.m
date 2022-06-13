classdef EuropeanOptionPricer < OptionPricer
    % @dev: abstract class to implement pricing routine for
    % european calls/ puts

    %% Protected Properties
    properties (Access = protected)
        model;
        instrumentList;
    end

    %% Constructor
    methods
        function this = EuropeanOptionPricer(valuationDate,instrumentList,model)
            this = this@OptionPricer(valuationDate)

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
        function calculatePrice(this,varargin)
            % dev: apply payoff and update price in the instrument list
            for i = 1:length(this.instrumentList)
                instrument = this.instrumentList(i);

                % dev: get instrument and log
                Logger.getInstance.log(LogType.INFO,...
                    ['Pricing instruement with id = ', num2str(instrument.getId),...
                    ' name = ',instrument.getName,' and type = ',instrument.getType]);

                % check if IR curve is loaded
                checkAndLoadIRCurve(this,instrument)

                if (isa(this.model,'SVCJ'))
                    price = computeSVCJPrice(this,instrument);

                elseif (isa(this.model,'BSM'))
                    price = computeBSMPrice(this,instrument,varargin{:});

                elseif (isa(this.model,'Heston'))
                    price = computeHestonPrice(this,instrument,varargin{:});

                elseif (isa(this.model,'Bates'))
                    price = computeBatesPrice(this,instrument,varargin{:});

                else
                    Logger.getInstance.log(LogType.FATAL,...
                        'Non-implemented model is set to EuropeanOptionPricer, cannot calculate price');
                end
                this.instrumentList(i).setPrice(price)
            end
        end

        % calculateImpliedVol
        function calculateImpliedVol(this,varargin)
            for i = 1:length(this.instrumentList)
                instrument = this.instrumentList(i);
                price = instrument.getPrice;
                strike = instrument.getStrike;
                [sInitial,valDate,div] = getUnderyingParam(this,instrument);
                % dev: get the number of periods that match maturity date
                maturityInYears = abs(days365(datenum(valDate),datenum(instrument.getMaturity)))/365;

                % check if IR curve is loaded
                checkAndLoadIRCurve(this,instrument)
                interestRate = this.irCurve.getZeroRates(instrument.getMaturity);

                type = instrument.getType;
                underlying = instrument.getUnderlying;
                if (isa(this.model,'BSM'))
                    impVol = this.model.getImpliedVol(underlying,price, strike,maturityInYears,interestRate, type, sInitial,div);

                else
                    Logger.getInstance.log(LogType.FATAL,...
                        'Non-implemented model to calculate implied vol set to EuropeanOptionPricer, cannot imply vol');
                end
                this.instrumentList(i).setImpliedVol(impVol);
            end

        end
    end

    %% Private methods
    methods (Access = private)

        % computeHestonPrice
        function price = computeHestonPrice(this,option,varargin)

            % get necessary properties
            strike = option.getStrike;
            [sInitial,valDate,div] = getUnderyingParam(this,option);
            % dev: get the number of periods that match maturity date and
            % turn into years
            maturityInYears = abs(days365(datenum(valDate),datenum(option.getMaturity)))/365;
            interestRate = this.irCurve.getZeroRates(option.getMaturity);

            if strcmp(option.getType,'call')
                price = this.model.getCallPrice(strike,maturityInYears,interestRate,sInitial,div);

            elseif strcmp(option.getType,'put')
                price = this.model.getPutPrice(strike,maturityInYears,interestRate,sInitial,div);

            else
                Logger.getInstance.log(LogType.FATAL,...
                    ['Type unknown for id = ', num2str(option.getId),...
                    ' name = ',option.getName]);

            end

        end

        % computeBatesPrice
        function price = computeBatesPrice(this,option,varargin)

            % get necessary properties
            strike = option.getStrike;
            [sInitial,valDate,div] = getUnderyingParam(this,option);
            % dev: get the number of periods that match maturity date and
            % turn into years
            maturityInYears = abs(days365(datenum(valDate),datenum(option.getMaturity)))/365;
            interestRate = this.irCurve.getZeroRates(option.getMaturity);

            if strcmp(option.getType,'call')
                price = this.model.getCallPrice(strike,maturityInYears,interestRate,sInitial,div);

            elseif strcmp(option.getType,'put')
                %price = this.model.getPutPrice(strike,maturityInPeriods,interestRate,varargin{:});
                Logger.getInstance.log(LogType.FATAL,...
                    'Bates pricing function for Puts not implemented yet');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    ['Type unknown for id = ', num2str(option.getId),...
                    ' name = ',option.getName]);

            end

        end

        % computeBSMPrice
        function price = computeBSMPrice(this,option,varargin)

            % get necessary properties
            strike = option.getStrike;
            [sInitial,valDate,div] = getUnderyingParam(this,option);
            % dev: get the number of periods that match maturity date
            maturityInYears = abs(days365(datenum(valDate),datenum(option.getMaturity)))/365;
            volatility = option.getImpliedVol;
            interestRate = this.irCurve.getZeroRates(option.getMaturity);
            underlying = option.getUnderlying;

            if strcmp(option.getType,'call')
                price = this.model.getCallPrice(underlying,strike,maturityInYears,interestRate,volatility,sInitial,div);

            elseif strcmp(option.getType,'put')
                price = this.model.getPutPrice(underlying,strike,maturityInYears,interestRate,volatility,sInitial,div);

            else
                Logger.getInstance.log(LogType.FATAL,...
                    ['Type unknown for id = ', num2str(option.getId),...
                    ' name = ',option.getName]);

            end

        end

        % computeSVCJPrice
        function price = computeSVCJPrice(this,option)

            [sInitial,t0] = getUnderyingParam(this);
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

        % applyCallPayoff
        function price = applyCallPayoff(this,simulatedPrices,maturityInPeriods,strike)
            terminalPrices = simulatedPrices(:,maturityInPeriods+1);
            % max(s-k,0)
            terminalPayoffs = max(terminalPrices - strike,0);
            % dev: TEMP DISCOUNTING
            price = mean(terminalPayoffs * exp(-0.03 * maturityInPeriods/365));
        end

        % applyPutPayoff
        function price = applyPutPayoff(this,simulatedPrices,maturityInPeriods,strike)

            terminalPrices = simulatedPrices(:,maturityInPeriods+1);
            % max(k-s,0)
            terminalPayoffs = max(strike -terminalPrices ,0);

            price = mean(terminalPayoffs);
        end

        % getUnderyingParam
        function [sInitial,valDate,div] = getUnderyingParam(this,option)

            if isa(option.getUnderlying,'Coin')
                sInitial = option.getUnderlying.getPriceTS.getLastValue;
                valDate = datestr(option.getUnderlying.getPriceTS.getLastDate);
                div = option.getUnderlying.getDivYield;

            elseif isa(option.getUnderlying,'Future')
                sInitial = option.getUnderlying.getPrice;
                valDate = datestr(option.getUnderlying.getValuationDate);
                div = option.getUnderlying.getUnderlying.getDivYield;
            end
        end

        % setModel
        function setModel(this,model)
            if (strcmp(model,'SVCJ')) % dev: check model name
                setSVCJModel(this);

            elseif (strcmp(model,'BSM'))
                mod = BSM();
                this.model = mod;
                Logger.getInstance.log(LogType.INFO,...
                    'BSM model is set to EuropeanOptionPricer');

            elseif (strcmp(model,'Heston'))
                setHestonModel(this)

            elseif (strcmp(model,'Bates'))
                setBatesModel(this)

            else
                Logger.getInstance.log(LogType.FATAL,...
                    [model, ' not implemented']);

            end
        end

        % setSVCJModel
        function setSVCJModel(this)
            if isempty(this.getModel)  % dev: if empty create and calibrate
                mod = SVCJ(this.underlying);
                mod.calibrate();
                if (mod.isCalibrated == ModelCalibration.SUCCESS)
                    Logger.getInstance.log(LogType.INFO,...
                        'SCVJ model is calibrated and set to EuropeanOptionPricer');
                else
                    Logger.getInstance.log(LogType.WARN,...
                        'Potentially non-clibrated SCVJ model set to EuropeanOptionPricer');
                end
                this.model = mod;
            elseif (isa(this.model,'SVCJ') && this.model.isCalibrated == ModelCalibration.NOT_CALIBRATED)
                this.model.calibrate();
                if (this.model.isCalibrated == ModelCalibration.SUCCESS)
                    Logger.getInstance.log(LogType.INFO,...
                        'SCVJ model is calibrated and set to EuropeanOptionPricer');
                else
                    Logger.getInstance.log(LogType.WARN,...
                        'Potentially non-clibrated SCVJ model set to EuropeanOptionPricer');
                end

            else
                Logger.getInstance.log(LogType.INFO,...
                    ['SVCJ model already set to EuropeanOptionPricer with calibration status: ', char(this.model.isCalibrated)]);
            end
        end

        % setHestonModel
        function setHestonModel(this)

            if isempty(this.getModel)  % dev: if empty create and calibrate
                mod = Heston();
                this.model = mod;
                Logger.getInstance.log(LogType.INFO,...
                    'Heston model parsed to EuropeanOptionPricer but not calibrated');
                % dev: get market data from instrument list and calibrate
                % model
                [marketPrice, strike, maturityInYears,interestRate,type,sInitial] = buildOptionMarketDataForCalibration(this);
                this.model.calibrate(marketPrice, strike, maturityInYears, interestRate,type, sInitial,0);

                if (this.model.isCalibrated == ModelCalibration.SUCCESS)
                    Logger.getInstance.log(LogType.INFO,...
                        'Heston model is calibrated and set to EuropeanOptionPricer');
                else
                    Logger.getInstance.log(LogType.WARN,...
                        'Potentially non-clibrated Heston model set to EuropeanOptionPricer');
                end

            elseif this.model.isCalibrated == ModelCalibration.NOT_CALIBRATED

                [marketPrice, strike, maturityInYears,interestRate,type,sInitial] = buildOptionMarketDataForCalibration(this);
                this.model.calibrate(marketPrice, strike, maturityInYears, interestRate,type, sInitial,0);
                if (this.model.isCalibrated == ModelCalibration.SUCCESS)
                    Logger.getInstance.log(LogType.INFO,...
                        'Heston model is calibrated and set to EuropeanOptionPricer');
                else
                    Logger.getInstance.log(LogType.WARN,...
                        'Potentially non-clibrated Heston model set to EuropeanOptionPricer');
                end

            else
                Logger.getInstance.log(LogType.INFO,...
                    ['Heston model already set to EuropeanOptionPricer with calibration status: ', char(this.model.isCalibrated)]);
            end
        end

        % setBatesModel
        function setBatesModel(this)

            if isempty(this.getModel)  % dev: if empty create and calibrate
                mod = Bates();
                this.model = mod;
                Logger.getInstance.log(LogType.INFO,...
                    'Bates model parsed to EuropeanOptionPricer but not calibrated');
                % dev: get market data from instrument list and calibrate
                % model
                [marketPrice, strike, maturityInYears,interestRate, type, sInitial] = buildOptionMarketDataForCalibration(this);
                this.model.calibrate(marketPrice, strike, maturityInYears, interestRate,type, sInitial,0);

                if (this.model.isCalibrated == ModelCalibration.SUCCESS)
                    Logger.getInstance.log(LogType.INFO,...
                        'Bates model is calibrated and set to EuropeanOptionPricer');
                else
                    Logger.getInstance.log(LogType.WARN,...
                        'Potentially non-clibrated Bates model set to EuropeanOptionPricer');
                end

            elseif this.model.isCalibrated == ModelCalibration.NOT_CALIBRATED

                [marketPrice, strike, maturityInYears,interestRate,type, sInitial] = buildOptionMarketDataForCalibration(this);
                this.model.calibrate(marketPrice, strike, maturityInYears, interestRate,type, sInitial,0);
                if (this.model.isCalibrated == ModelCalibration.SUCCESS)
                    Logger.getInstance.log(LogType.INFO,...
                        'Bates model is calibrated and set to EuropeanOptionPricer');
                else
                    Logger.getInstance.log(LogType.WARN,...
                        'Potentially non-clibrated Bates model set to EuropeanOptionPricer');
                end

            else
                Logger.getInstance.log(LogType.INFO,...
                    ['Bates model already set to EuropeanOptionPricer with calibration status: ', char(this.model.isCalibrated)]);
            end
        end

        % buildOptionMarketDataForCalibration
        function [marketPrice, strike, maturityInYears,interestRate,type, sInitial] = buildOptionMarketDataForCalibration(this)

            nosOfOptions = length(this.instrumentList);
            strike = nan(nosOfOptions,1);
            marketPrice = nan(nosOfOptions,1);
            maturityInYears = nan(nosOfOptions,1);
            interestRate = nan(nosOfOptions,1);
            sInitial = nan(nosOfOptions,1);

            for i = 1:nosOfOptions
                instrument = this.instrumentList(i);
                [underlyingPrice,t0] = getUnderyingParam(this,instrument);

                strike(i) = instrument.getStrike;
                marketPrice(i) = instrument.getPrice;

                % check if IR curve is loaded
                checkAndLoadIRCurve(this,instrument)

                interestRate(i) = this.irCurve.getZeroRates(instrument.getMaturity);
                type{i,1} = instrument.getType;
                sInitial(i) = underlyingPrice;

                % dev: get the number of periods that match maturity date
                maturityInYears(i) = abs(days365(datenum(t0),datenum(instrument.getMaturity)))/365;

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

