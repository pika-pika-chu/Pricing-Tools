classdef Heston < handle
    % @dev: class for Heston 1993 implementation

    %% Private Properties
    properties (Access = private)

        initialParamConfig = struct(...
            'v0',[],...
            'theta',[],...
            'rho',[],...
            'kappa',[],...
            'sigma',[],...
            'alpha0',[]);
        calibrated = ModelCalibration.NOT_CALIBRATED;

        calibratedParam = struct(...
            'v0',[],...
            'theta',[],...
            'rho',[],...
            'kappa',[],...
            'sigma',[],...
            'alpha',[]);

    end

    %% Constructor
    methods
        function this = Heston()
            % dev: check and load inital param for calibration
            checkAndLoadModelConfig(this);
        end
    end


    %% Public methods
    methods (Access = public)

        % getCallPrice
        function [prices,alphas] = getCallPrice(this,strike,timeToMaturity,interestRate,sInitial,div)

            if (this.isCalibrated == ModelCalibration.NOT_CALIBRATED ||...
                    this.isCalibrated == ModelCalibration.FAILURE)
                Logger.getInstance.log(LogType.FATAL,...
                    'Model not calibrated, calibrate before calling pricing methods!');

            end

            [prices,alphas] = calculateCallPrice(this,strike,timeToMaturity,interestRate,sInitial,div,...
                this.calibratedParam.v0,...
                this.calibratedParam.theta,...
                this.calibratedParam.rho,...
                this.calibratedParam.kappa,...
                this.calibratedParam.sigma,...
                this.initialParamConfig.alpha0);
        end

        % getPutPrice
        function [prices,alphas] = getPutPrice(this,strike,timeToMaturity,interestRate,sInitial,div)

            if (isCalibrated == ModelCalibration.NOT_CALIBRATED ||...
                    isCalibrated == ModelCalibration.FAILURE)
                Logger.getInstance.log(LogType.FATAL,...
                    'Model not calibrated, calibrate before calling pricing methods!');

            end


            [prices,alphas] = calculatePutPrice(this,strike,timeToMaturity,interestRate,sInitial,div,...
                this.calibratedParam.v0,...
                this.calibratedParam.theta,...
                this.calibratedParam.rho,...
                this.calibratedParam.kappa,...
                this.calibratedParam.sigma,...
                this.initialParamConfig.alpha0);
        end

        % calibrate
        function calibrate(this,marketPrice, strike, timeToMaturity, interestRate,type, sInitial,div)

            if (this.isCalibrated == ModelCalibration.SUCCESS)
                Logger.getInstance.log(LogType.INFO,...
                    'Model already calibrated!');

            else

                calibrateModel(this,marketPrice,strike,timeToMaturity,interestRate,type,sInitial,div);
                this.calibrated = ModelCalibration.SUCCESS;
            end
        end
    end

    %% Private methods
    methods (Access = private)

        % calibrateModel
        function calibrateModel(this,marketPrice,strike,timeToMaturity,interestRate,type,sInitial,div)

            if (sum(strcmp(type, type{1})) ~= length(type))

                Logger.getInstance.log(LogType.FATAL,...
                    'Use identical type of options to calibrate, either all calls or puts');
            end

            iniParams = [this.initialParamConfig.v0,...
                this.initialParamConfig.theta,...
                this.initialParamConfig.rho,...
                this.initialParamConfig.kappa,...
                this.initialParamConfig.sigma];
            alpha0 = this.initialParamConfig.alpha0;
            options = optimoptions('lsqnonlin','Display', 'iter'); %'Algorithm','levenberg-marquardt'

            if (sum(strcmp(type, 'call')) == length(type))
                    
                 paramsOut = fmincon(@(x) sum((marketPrice - ...
                    this.calculateCallPrice(strike,timeToMaturity,interestRate,sInitial,div,x(1),x(2),x(3),x(4),x(5),alpha0)).^2),...
                    iniParams,[],[],[],[],...
                    [0 0 -0.9 0 0], ... % dev:LB calibrated params
                    [1 5 0.9 20 5]);
%                 paramsOut = lsqnonlin(@(x) (marketPrice - ...
%                     this.calculateCallPrice(strike,timeToMaturity,interestRate,sInitial,div,x(1),x(2),x(3),x(4),x(5),alpha0)),...
%                     iniParams,...
%                     [0 0 -0.9 0 0], ... % dev:LB calibrated params
%                     [1 3 0.9 100 2], ...  % dev:UB for calibrated params
%                     options);
            else
                paramsOut = fmincon(@(x) sum((marketPrice - ...
                    this.calculatePutPrice(strike,timeToMaturity,interestRate,sInitial,div,x(1),x(2),x(3),x(4),x(5),alpha0)).^2),...
                    iniParams,[],[],[],[],...
                    [0 0 -0.9 0 0], ... % dev:LB calibrated params
                    [1 3 0.9 100 2]);
                
                
%                 lsqnonlin(@(x) sum(abs(marketPrice - ...
%                     this.calculatePutPrice(strike,timeToMaturity,interestRate,sInitial,div,x(1),x(2),x(3),x(4),x(5),alpha0))./marketPrice),...
%                     iniParams,...
%                     [0 0 -0.9 0 0], ... % dev:LB calibrated params
%                     [1 1 0.9 100 1], ...  % dev:UB for calibrated params
%                     options);

            end
            this.calibratedParam.v0 = paramsOut(1);
            this.calibratedParam.theta = paramsOut(2);
            this.calibratedParam.rho = paramsOut(3);
            this.calibratedParam.kappa = paramsOut(4);
            this.calibratedParam.sigma = paramsOut(5);
        end

        % calculateCallPrice
        function [prices,alphas] = calculateCallPrice(this,strike,timeToMaturity,interestRate,sInitial,div,...
                v0,theta,rho,kappa,sigma,alpha)

            % check total num of inputs
            nosOfOptions = length(strike);

            % ctrate nan vec to store results
            prices = nan(nosOfOptions,1);
            alphas = nan(nosOfOptions,1);

            % calculate forward S
            mu = interestRate - div;
            fwds = sInitial.*exp(mu.*timeToMaturity);

            for i = 1:nosOfOptions
                try
                    % using fzero here instead of fminsearch
                    alphaPrime = fzero( @(a) this.psi(a,strike(i), fwds(i), kappa, theta, rho, sigma, timeToMaturity(i), v0), alpha);
                    alphas(i) = alphaPrime;
                catch
                    alphas(i) = alpha;

                end
                prices(i) =  this.Ralpha(fwds(i), strike(i), alphas(i))+1/pi* ...
                    integral(@(x) this.phi(x, strike(i), alphas(i), fwds(i), kappa, theta, rho, sigma, timeToMaturity(i), v0) , 0, Inf);
            end

        end

        % calculatePutPrice
        function [prices,alphas] = calculatePutPrice(this,strike,timeToMaturity,interestRate,sInitial,div,...
                v0,theta,rho,kappa,sigma,alpha)

            % check total num of inputs
            nosOfOptions = length(strike);

            % ctrate nan vec to store results
            prices = nan(nosOfOptions,1);
            alphas = nan(nosOfOptions,1);

            % calculate forward S
            mu = interestRate - div;
            fwds = sInitial.*exp(mu.*timeToMaturity);

            for i = 1:nosOfOptions
                try
                    % using fzero here instead of fminsearch
                    alphaPrime = fzero( @(a) this.psi(a,strike(i), fwds(i), kappa, theta, rho, sigma, timeToMaturity(i), v0), alpha);
                    alphas(i) = alphaPrime;
                catch
                    alphas(i) = alpha;

                end
                prices(i) =  this.Ralpha(fwds(i), strike(i), alphas(i))+1/pi* ...
                    integral(@(x) this.phi(x, strike(i), alphas(i), fwds(i), kappa, theta, rho, sigma, timeToMaturity(i), v0) , 0, Inf);

                prices(i) = prices(i) + strike(i) * exp(-interestRate(i) * ...
                    timeToMaturity(i)) - sInitial(i) * exp(-div(i)*timeToMaturity(i));
            end

        end
    end


    %% Private methods
    methods (Access = private)

        % psi
        function p = psi(this,alpha, K, F, kappa, theta, rho, sigma, tau, v0)
            k = log(K);
            p = -alpha*k+0.5*log(this.phi(-(alpha+1)*1i, K, alpha, F, kappa, theta, rho, sigma, tau, v0)^2);
        end

        % Ralpha
        function r = Ralpha(this,F, K, alpha)
            r = F*(alpha<=0)-K*(alpha<=-1)-0.5*(F*(alpha==0)-K*(alpha==-1));
        end

        % phi
        function y = phi(this,v, K, alpha, F, kappa, theta, rho, sigma, tau, v0)
            k = log(K);
            y = real(exp(-1i*(v-1i*alpha)*k).*( this.cf(v-1i*(alpha+1), F, kappa, theta, rho, sigma, tau, v0)./(-(v-1i*(alpha+1)).*(v-1i*alpha))));
        end

        % cf
        function c = cf(this,u, F, kappa, theta, rho, sigma, tau, v0)
            f = log(F);
            c = exp(1i*u*f+ this.A(u, kappa, theta, rho, sigma, tau) + this.Bv(u, rho, sigma, kappa, tau)*v0);
        end

        % Bv
        function b = Bv(this,u, rho, sigma, kappa, tau)
            b = ((this.beta(u,rho,sigma,kappa) - ...
                this.D(u, rho, sigma, kappa)).*(1-exp(-this.D(u, rho, sigma, kappa)*tau)))./...
                (sigma.^2 * ( 1 - this.G(u, rho, sigma, kappa).* exp(- this.D(u, rho, sigma, kappa)*tau)));
        end

        % A
        function a = A(this,u, kappa, theta, rho, sigma, tau)
            a = (kappa * theta * ((this.beta(u,rho,sigma,kappa) - this.D(u, rho, sigma, kappa)) *...
                tau-2*log(this.phi2(u, rho, sigma, kappa, tau))))/sigma.^2;
        end

        % phi2
        function p = phi2(this,u, rho, sigma, kappa, tau)
            p = (this.G(u, rho, sigma, kappa).*exp(-this.D(u, rho, sigma, kappa)*tau)-1)./...
                (this.G(u, rho, sigma, kappa)-1);
        end

        % G
        function g = G(this,u, rho, sigma, kappa)
            g = (this.beta(u,rho,sigma,kappa) - this.D(u, rho, sigma, kappa))./...
                (this.beta(u,rho,sigma,kappa) + this.D(u, rho, sigma, kappa));
        end

        % D
        function d = D(this,u, rho, sigma, kappa)
            d = sqrt(this.beta(u,rho,sigma,kappa).^2-4 * this.alphahat(u) * this.gamma(sigma));
        end

        % alphahat
        function a = alphahat(this,u)
            a = -0.5*u.*(1i+u);
        end

        % beta
        function b = beta(this,u,rho,sigma,kappa)
            b = kappa-rho*sigma*u*1i;
        end

        % gamma
        function y = gamma(this,sigma)
            y = 0.5*sigma.^2;
        end

        % checkAndLoadModelConfig
        function checkAndLoadModelConfig(this)
            modelConfig = getConfig('Heston','calibration');

            % v0
            if isfield(modelConfig,'v0')
                this.initialParamConfig.v0 = modelConfig.v0;
                Logger.getInstance.log(LogType.INFO,...
                    'v0 provided, initial vol set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'v0 not provided, initial vol missing');
            end

            % theta
            if isfield(modelConfig,'theta')
                this.initialParamConfig.theta = modelConfig.theta;
                Logger.getInstance.log(LogType.INFO,...
                    'theta provided, theta set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'theta not provided, theta missing');
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

            % alpha0
            if isfield(modelConfig,'alpha0')
                this.initialParamConfig.alpha0 = modelConfig.alpha0;
                Logger.getInstance.log(LogType.INFO,...
                    'alpha0 provided, alpha0 set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'alpha0 not provided, alpha0 missing');
            end
        end

        %getVararginValues
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

    %% Public methods
    methods (Access = public)

        %isCalibrated
        function calibrationStatus = isCalibrated(this)
            calibrationStatus = this.calibrated;
        end

        %getCalibratedParam
        function calibratedParam = getCalibratedParam(this)
            calibratedParam = this.calibratedParam;
        end
    end
end

