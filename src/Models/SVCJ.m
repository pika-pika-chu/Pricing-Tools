classdef SVCJ < handle
    % @dev: class for SVCJ

    %% Private Properties
    properties (Access = private)
        priceTS;
        returnTS;
        modelConfig = struct(...
            'returnMethod',[],...
            'annulisationFactor',[],...
            'NumRuns',[],...
            'burnInPeriod',[]);
        calibratedParmeters = struct(...
            'mu',[],...
            'muY',[],...
            'sigmaY',[],...
            'alpha',[],...
            'beta',[],...
            'rho',[],...
            'sigmaV',[],...
            'rhoJ',[],... 
            'muV',[],...
            'lambda',[]);
        isCalibrated = ModelCalibration.NOT_CALIBRATED;
        simulatedReturnsTS; % dev: getter only
        simulatedVolatilityTS; % dev: getter only
    end

    %% Constructor
    methods
        function this = SVCJ(coin,modelConfig)

            % dev: store the price TS from coin
            this.priceTS = coin.getPriceTS;

            % dev: cehck and load modelConfig to properties
            checkAndLoadModelConfig(this,modelConfig);

            % dev: source the return timeseries and store
            sourceReturnTimeseries(this)
        end

    end

    %% Public methods
    methods (Access = public)
        %calibrate
        function calibrate(this)
            calibratedParams = calibrateModel(this);
            this.calibratedParmeters = calibratedParams;
        end

        % simulateReturns
        function simulateReturns(this)

        end
    end

    %% Private methods
    methods (Access = private)
        % calibrateModel
        function calibratedParams = calibrateModel(this)
             % get initial distribution values    
            distVal = fetchDistributionValues(this);
            
            % dev:
            % | Return process (Y) |
            % | Volatility process (V) |

            % dev: store muY
            % muY = mean of return process (Y)
            muY = distVal.a; muYSum = 0 ;muYSqSum = 0;

            % dev: store kappa
            kappa = 0; kappaSum = 0 ; kappaSqSum = 0;

            % dev: store alpha
            b = distVal.b; alpha = b(1); alphaSum = 0; alphaSqSum = 0;

            % dev: store beta
            beta = b(2); betaSum = 0 ;betaSqSum = 0;

            % dev: store sigmaV
            % sigmaV = vol of vol process (V)
            C = distVal.C; c = distVal.c;
            sigmaV = C/(c - 2); sigmaVSum = 0; sigmaVSqSum = 0;

            % dev: store rho
            % rho = corr between W1 and W2 process used for return process (Y) 
            % and vol process (V)
            rho = 0; rhoSum = 0; rhoSqSum = 0; 
            
            % dev: store muV_J
            % muV_J = exp param for volatility (process V) jump distribution
            D = distVal.D; d = distVal.d;
            muV_J = D/(d-2); muV_JSum = 0; muV_JSqSum = 0;

            % dev: store muY_J
            % muY_J = mean of jump size in return process (Y) 
            e = distVal.e;
            muY_J = e; muY_JSum = 0; muY_JSqSum = 0;

            % dev: store sigmaY_J
            % sigmaY_J = vol of jump in return (process Y)
            F = distVal.F; f = distVal.f;
            sigmaY_J = F/(f - 2); sigmaY_JSum = 0; sigmaY_JSqSum = 0;

            % dev: store rhoJ
            % rhoJ = corr between return jump process and vol jump process
            g = distVal.g;
            rhoJ = g; rhoJSum = 0; rhoJSqSum = 0;

            % dev: store lambda
            % lambda = jump intensity
            lambda = 0; lambdaSum = 0; lambdaSqSum = 0;

            % get returns as ret
            ret = this.returnTS.getValues;
            numOfRet = length(ret);
            
            % initial values for variance_t;
            V = 0.1*(ret - mean(ret)).^2 + 0.9*var(ret); 
            % estimate jumps
            J = abs(ret) - mean(ret) > 2 * std(ret);
            JSum = 0;

            % get jump size in volatility (V), zV
            zV = exprnd(muV_J,numOfRet,1); 
            zVSum = 0;
            
            % the jump size in return (Y), zY
            zY = mvnrnd((muY_J + zV * rhoJ), sigmaY_J); 
            zYSum = 0;
            
            % ===== CHECK WHERE THEY ARE COMING FROM 
            % 
            stdevrho = 0.01;
            dfrho = 6.5;
            stdevV = 0.9;
            dfV = 4.5;
            acceptsumV = zeros(size(V));
            acceptsumrho = 0;
            acceptsums2V = 0;
            Z = ones(T,1);

            numRuns = this.modelConfig.NumRuns;
            numBurnIn = this.modelConfig.burnInPeriod;




            if (true)
                this.isCalibrated = ModelCalibration.SUCCESS;
            else
                this.isCalibrated = ModelCalibration.FAILURE;
            end
        end


        % checkAndLoadModelConfig
        function checkAndLoadModelConfig(this,modelConfig)
            % returnMethod
            if isfield(modelConfig,'returnMethod')
                this.modelConfig.returnMethod = modelConfig.returnMethod;
                Logger.getInstance.log(LogType.INFO,...
                    'modelConfig provided, returnMethod set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Incomplete modelConfig provided, returnMethod missing');
            end

            % annulisationFactor
            if isfield(modelConfig,'annulisationFactor')
                this.modelConfig.annulisationFactor = modelConfig.annulisationFactor;
                Logger.getInstance.log(LogType.INFO,...
                    'modelConfig provided, annulisationFactor set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Incomplete modelConfig provided, annulisationFactor missing');
            end

            % NumRuns
            if isfield(modelConfig,'NumRuns')
                this.modelConfig.NumRuns = modelConfig.NumRuns;
                Logger.getInstance.log(LogType.INFO,...
                    'modelConfig provided, NumRuns set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Incomplete modelConfig provided, NumRuns missing');
            end

            % burnInPeriod
            if isfield(modelConfig,'burnInPeriod')
                this.modelConfig.burnInPeriod = modelConfig.burnInPeriod;
                Logger.getInstance.log(LogType.INFO,...
                    'modelConfig provided, burnInPeriod set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Incomplete modelConfig provided, burnInPeriod missing');
            end
        end

        % sourceReturnTimeseries
        function sourceReturnTimeseries(this)
            this.returnTS = this.priceTS.getReturnTimeseries(this.modelConfig.returnMethod,'annulisationFactor',...
                this.modelConfig.annulisationFactor);
        end

        % fetchDistributionValues
        function iniDistVal = fetchDistributionValues(this)

            % dev: prior for mu
            iniDistVal.a = 0; iniDistVal.A = 25;
            % dev: prior for alpha and beta
            b = zeros(2,1); iniDistVal.b = b; 
            B = 1*eye(length(b)); iniDistVal.B = B;
            % dev: prior for sigmaV
            iniDistVal.c = 2.5; iniDistVal.C = 0.1;
            % dev: prior for muV
            iniDistVal.d = 10; iniDistVal.D = 20;
            % dev: prior for muY
            iniDistVal.e = 0; iniDistVal.E = 100;
            % dev: prior for sigmaY
            iniDistVal.f = 10; iniDistVal.F = 40;
            % dev: prior for rhoY
            iniDistVal.g = 0; iniDistVal.G = 4;
            % dev: prior for lambda
            iniDistVal.k = 2; iniDistVal.K = 40;

        end
    end

    %% Getters and Setters
    methods (Access = public)

        %getPriceTS
        function priceTS = getPriceTS(this)
            priceTS = this.priceTS;
        end

        %setPriceTS
        function setPriceTS(this,priceTS)
            % dev: check if parsed time-series isa TimeSeries type
            % check type = TimeseriesType.PRICE
            if (isa(priceTS,'TimeSeries') && ...
                    priceTS.getType == TimeseriesType.PRICE)
                this.priceTS = priceTS;
                Logger.getInstance.log(LogType.INFO,...
                    'Price timeseries provided and set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Parse timeseries not of Timeseries.PRICE type');
            end

        end

        %getReturnTS
        function returnTS = getReturnTS(this)
            returnTS = this.returnTS;
        end

        %setReturnTS
        function setReturnTS(this,returnTS)
             % dev: check if parsed time-series isa TimeSeries type
             % check type = TimeseriesType.RETURN
            if (isa(returnTS,'TimeSeries') && ...
                    returnTS.getType == TimeseriesType.RETURN)
                this.returnTS = returnTS;
                Logger.getInstance.log(LogType.INFO,...
                    'Price timeseries provided and set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Parse timeseries not of Timeseries.PRICE type');
            end
        end

        %getModelConfig
        function modelConfig = getModelConfig(this)
            modelConfig = this.modelConfig;
        end

        %setModelConfig
        function setModelConfig(this,modelConfig)
            this.modelConfig = modelConfig;
        end

        %getSimulatedReturnsTS
        function simulatedReturnsTS = getSimulatedReturnsTS(this)
            simulatedReturnsTS = this.simulatedReturnsTS;
        end

        %getSimulatedVolatilityTS
        function simulatedVolatilityTS = getSimulatedVolatilityTS(this)
            simulatedVolatilityTS = this.simulatedVolatilityTS;
        end
    end
end