classdef Bates < handle
    % @dev: class for Bates 1996 implementation

    %% Private Properties
    properties (Access = private)
        underlying;
        calibrated = ModelCalibration.NOT_CALIBRATED;

        initialParamConfig = struct(...
            'v0',[],...
            'vT',[],...
            'rho',[],...
            'kappa',[],...
            'sigma',[],...
            'lambda',[],...
            'muJ',[],...
            'sigmaJ',[]);

        calibratedParam = struct(...
            'v0',[],...
            'vT',[],...
            'rho',[],...
            'kappa',[],...
            'sigma',[],...
            'lambda',[],...
            'muJ',[],...
            'sigmaJ',[]);
    end

    %% Constructor
    methods
        function this = Bates(coin)
            this.underlying = coin;

            % dev: check and load inital param for calibration
            checkAndLoadModelConfig(this);
        end
    end

    %% Public methods
    methods (Access = public)
        % getCallPrice
        function price = getCallPrice(this,strike,timeToMaturity,interestRate,varargin)

            if (this.isCalibrated == ModelCalibration.NOT_CALIBRATED ||...
                    this.isCalibrated == ModelCalibration.FAILURE)
                Logger.getInstance.log(LogType.FATAL,...
                    'Model not calibrated, calibrate before calling pricing methods!');

            end

            if ~isempty(varargin)
                [sInitial,div] = getVararginValues(this,varargin{:});

            else
                sInitial = this.underlying.getPriceTS.getLastValue;
                div = 0;

            end

            price = calculateCallPrice(this,strike,timeToMaturity,interestRate,sInitial,div,...
                this.calibratedParam.v0,...
                this.calibratedParam.vT,...
                this.calibratedParam.rho,...
                this.calibratedParam.kappa,...
                this.calibratedParam.sigma,...
                this.calibratedParam.lambda,...
                this.calibratedParam.muJ,...
                this.calibratedParam.sigmaJ);
        end

        % calibrate
        function calibrate(this,marketPrice, strike, timeToMaturity, interestRate,type, varargin)

            if (this.isCalibrated == ModelCalibration.SUCCESS)
                Logger.getInstance.log(LogType.INFO,...
                    'Model already calibrated!');

            else
                if ~isempty(varargin)
                    [sInitial,div] = getVararginValues(this,varargin{:});

                else
                    sInitial = this.underlying.getPriceTS.getLastValue;
                    div = 0;

                end
                calibrateModel(this,marketPrice,strike,timeToMaturity,interestRate,type,sInitial,div);
                this.calibrated = ModelCalibration.SUCCESS;
            end
        end

    end


    %% Public methods
    methods (Access = public)

        %getUnderlying
        function underlying = getUnderlying(this)
            underlying = this.underlying;
        end

        %setUnderlying
        function setUnderlying(this,underlying)
            this.underlying = underlying;
        end

        %isCalibrated
        function calibrationStatus = isCalibrated(this)
            calibrationStatus = this.calibrated;
        end
    end

    %% Private methods
    methods (Access = private)

        % calculateCallPrice
        function prices = calculateCallPrice(this,strike,timeToMaturity,interestRate,sInitial,div,...
                v0,vT,rho,kappa,sigma,lambda,muJ,sigmaJ)

            % check total num of inputs
            nosOfOptions = length(strike);

            % create nan vec to store results
            prices = nan(nosOfOptions,1);

            for i = 1:nosOfOptions

                vP1 = 0.5 + 1/pi * quadl(@this.P1, 0, 200, [], [], sInitial, strike(i),...
                    timeToMaturity(i), interestRate(i), div,...
                    v0 , vT, rho, kappa, sigma, lambda, muJ, sigmaJ);

                vP2 = 0.5 + 1/pi * quadl(@this.P2, 0, 200, [], [], sInitial, strike(i),...
                    timeToMaturity(i), interestRate(i), div,...
                    v0, vT, rho, kappa, sigma, lambda, muJ, sigmaJ);

                prices(i) = exp(-div * timeToMaturity(i)) * sInitial *...
                    vP1 - exp(-interestRate(i) * timeToMaturity(i)) * strike(i) * vP2;
            end

        end

        % P1
        function p = P1(this, om, sInitial, strike, tau, interestRate, div, v0, vT, rho, kappa, sigma, lambda, muJ, sigmaJ)
            i = 1i;
            p = real(exp(-i * log(strike) * om).*...
                this.cfBates(om-i, sInitial, tau, interestRate, div, v0, vT, rho,...
                kappa, sigma, lambda, muJ, sigmaJ)./(i * om * sInitial * exp((interestRate-div) * tau)));
        end

        % P2
        function p = P2(this, om, sInitial, strike, tau, interestRate, div, v0, vT, rho, kappa, sigma, lambda, muJ, sigmaJ)
            i = 1i;
            p = real(exp(-i * log(strike) * om).*...
                this.cfBates(om, sInitial, tau, interestRate, div, v0, vT, rho,...
                kappa, sigma, lambda, muJ, sigmaJ)./(i * om));
        end

        % cfBates
        function cf = cfBates(this, om, sInitial, tau, interestRate, div, v0, vT, rho, kappa, sigma, lambda, muJ, sigmaJ)

            d = sqrt((rho * sigma * 1i * om - kappa).^2 + sigma^2 * (1i * om + om.^2));

            g2 = (kappa - rho * sigma * 1i * om - d)./(kappa - rho * sigma * 1i * om + d);

            cf1 = 1i * om.*(log(sInitial) + (interestRate - div)* tau);

            cf2 = vT * kappa/(sigma^2) * ((kappa - rho * sigma * 1i * om - d) * ...
                tau - 2 * log((1 - g2.*exp(-d * tau))./(1 - g2)));

            cf3 = v0/sigma^2 * (kappa - rho * sigma * 1i * om - d).*(1 - exp(-d * tau))./(1 - g2.*exp(-d * tau));

            %jump
            cf4 = -lambda * muJ * 1i * tau * om + lambda * tau * ((1 + muJ).^(1i * om).*exp(sigmaJ * (1i * om/2).*(1i * om-1)) - 1);

            cf = exp(cf1 + cf2 + cf3 + cf4);
        end
    
        % calibrateModel
        function calibrateModel(this,marketPrice,strike,timeToMaturity,interestRate,type,sInitial,div)

            if (sum(strcmp(type, type{1})) ~= length(type))

                Logger.getInstance.log(LogType.FATAL,...
                    'Use identical type of options to calibrate, either all calls or puts');
            end

            iniParams = [this.initialParamConfig.v0,...
                this.initialParamConfig.vT,...
                this.initialParamConfig.rho,...
                this.initialParamConfig.kappa,...
                this.initialParamConfig.sigma,...
                this.initialParamConfig.lambda,...
                this.initialParamConfig.muJ,...
                this.initialParamConfig.sigmaJ];
          
            options = optimoptions('lsqnonlin','Algorithm','levenberg-marquardt', 'Display', 'iter');

            if (sum(strcmp(type, 'call')) == length(type))
                paramsOut = lsqnonlin(@(x) this.calculateCallPrice(strike,timeToMaturity,interestRate,sInitial,div,x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8)) - marketPrice,...
                    iniParams,...
                    [eps eps -1+eps eps eps eps eps eps ], ... % dev:LB calibrated params
                    [Inf Inf 1-eps Inf  Inf Inf Inf Inf ], ...  % dev:UB for calibrated params
                    options);
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Puts pricing routine not implemented for Bates, cannot calibrate model to put prices');

            end
            this.calibratedParam.v0 = paramsOut(1);
            this.calibratedParam.vT = paramsOut(2);
            this.calibratedParam.rho = paramsOut(3);
            this.calibratedParam.kappa = paramsOut(4);
            this.calibratedParam.sigma = paramsOut(5);
            this.calibratedParam.lambda = paramsOut(6);
            this.calibratedParam.muJ = paramsOut(7);
            this.calibratedParam.sigmaJ = paramsOut(8);
        end
    end
    %% Private methods
    methods (Access = private)

        % checkAndLoadModelConfig
        function checkAndLoadModelConfig(this)
            modelConfig = getConfig('Bates','calibration');

            % v0
            if isfield(modelConfig,'v0')
                this.initialParamConfig.v0 = modelConfig.v0;
                Logger.getInstance.log(LogType.INFO,...
                    'v0 provided, initial vol set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'v0 not provided, initial vol missing');
            end

            % vT
            if isfield(modelConfig,'vT')
                this.initialParamConfig.vT = modelConfig.vT;
                Logger.getInstance.log(LogType.INFO,...
                    'vT provided, vT set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'vT not provided, vT missing');
            end

            % rho
            if isfield(modelConfig,'rho')
                this.initialParamConfig.rho = modelConfig.rho;
                Logger.getInstance.log(LogType.INFO,...
                    'rho provided, rho set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'rho not provided, rho missing');
            end

            % kappa
            if isfield(modelConfig,'kappa')
                this.initialParamConfig.kappa = modelConfig.kappa;
                Logger.getInstance.log(LogType.INFO,...
                    'kappa provided, kappa set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'kappa not provided, kappa missing');
            end

            % sigma
            if isfield(modelConfig,'sigma')
                this.initialParamConfig.sigma = modelConfig.sigma;
                Logger.getInstance.log(LogType.INFO,...
                    'sigma provided, vol of vol set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'sigma not provided, vol of vol missing');
            end

            % lambda
            if isfield(modelConfig,'lambda')
                this.initialParamConfig.lambda = modelConfig.lambda;
                Logger.getInstance.log(LogType.INFO,...
                    'lambda provided, lambda set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'lambda not provided, lambda missing');
            end

            % muJ
            if isfield(modelConfig,'muJ')
                this.initialParamConfig.muJ = modelConfig.muJ;
                Logger.getInstance.log(LogType.INFO,...
                    'muJ provided, muJ set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'muJ not provided, muJ missing');
            end

            % sigmaJ
            if isfield(modelConfig,'sigmaJ')
                this.initialParamConfig.sigmaJ = modelConfig.lambda;
                Logger.getInstance.log(LogType.INFO,...
                    'sigmaJ provided, sigmaJ set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'sigmaJ not provided, sigmaJ missing');
            end
        end

        % getVararginValues
        function [sInitial,div] = getVararginValues(this,varargin)
            % set empty values before checking
            sInitial = [];
            div = [];
            % check if name value arg is present in varargin
            sInitialIdx = find(strcmp(varargin,'sInitial'),1);
            divIdx = find(strcmp(varargin,'dividend'),1);
            % set name value args if present in varargin
            if ~isempty(sInitialIdx)
                sInitial = varargin{sInitialIdx + 1};
            end

            if ~isempty(divIdx)
                div = varargin{divIdx + 1};
            end

        end
    end
end

