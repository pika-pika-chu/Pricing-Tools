classdef SVCJ < handle
    % @dev: class for SVCJ

    %% Private Properties
    properties (Access = private)
        priceTS;
        returnTS;
        calibrationConfig = struct(...
            'returnMethod',[],...
            'annulisationFactor',[],...
            'NumRuns',[],...
            'burnInPeriod',[]);
        calibratedParameters = struct(...
            'muY',[],... // mean of return process
            'muY_J',[],... // mean on return jump process
            'sigmaY_J',[],... // vol of return jump process
            'alpha',[],...
            'beta',[],...
            'kappa',[],...
            'rho',[],...
            'sigmaV',[],...// vol of vol
            'rhoJ',[],... // corr between jump process
            'muV_J',[],... // exp param for variance jump
            'lambda',[]);
        calibrated = ModelCalibration.NOT_CALIBRATED;
        simulatedReturns; % dev: getter only
        simulatedVolatility; % dev: getter only
        simulatedJumpsInReturns; % dev: getter only
        simulatedJumpsInVol; % dev: getter only
        simulatedPrices; % dev: getter only

        traceParamSet= struct(...% dev: to store estimated parameter evolution
            'muY',[],... // mean of return process
            'muY_J',[],... // mean on return jump process
            'sigmaY_J',[],... // vol of return jump process
            'alpha',[],...
            'beta',[],...
            'kappa',[],...
            'rho',[],...
            'sigmaV',[],...// vol of vol
            'rhoJ',[],... // corr between jump process
            'muV_J',[],... // exp param for variance jump
            'lambda',[],...
            'iterCount',[]);

        errMeasureParam = struct(...
            'Jsum',[],... // jump accpeted values summed
            'zYsum',[],... // return process accepted values summed
            'zVsum',[],... // vol process accepted values summed
            'varianceSum',[],...
            'iterCount',[]); % sample size
    end

    %% Constructor
    methods
        function this = SVCJ(coin)

            % dev: store the price TS from coin
            this.priceTS = coin.getPriceTS;

            % dev: check and load calibration config to properties
            checkAndLoadModelConfig(this);

            % dev: source the return timeseries and store
            sourceReturnTimeseries(this)
        end

    end

    %% Public methods
    methods (Access = public)
        %calibrate
        function calibrate(this)
            calibrateModel(this);

        end

        % simulateReturns
        function simulate(this)
            simulatePrices(this);
        end
    end

    %% Private methods
    methods (Access = private)
        % calibrateModel
        function calibrateModel(this)
            Logger.getInstance.log(LogType.INFO,...
                'SVCJ Model calibration started, this may take time!!!');

            numRuns = this.calibrationConfig.NumRuns;
            numBurnIn = this.calibrationConfig.burnInPeriod;

            % get initial distribution values
            distVal = fetchDistributionValues(this);
            tmp=[];
            % dev:
            % | Return process (Y) |
            % | Volatility process (variance) |

            % dev: store muY
            % muY = mean of return process (Y)
            muY = distVal.a; muYSum = 0 ;muYSqSum = 0;
            muYEvolve = nan(numRuns-numBurnIn,1);
            muYMeanEvolve = nan(numRuns-numBurnIn,1);

            % dev: store kappa
            kappa = 0; kappaSum = 0 ; kappaSqSum = 0;
            kappaEvolve = nan(numRuns-numBurnIn,1);
            kappaMeanEvolve = nan(numRuns-numBurnIn,1);

            % dev: store alpha
            b = distVal.b; alpha = b(1); alphaSum = 0; alphaSqSum = 0;
            alphaEvolve = nan(numRuns-numBurnIn,1);
            alphaMeanEvolve = nan(numRuns-numBurnIn,1);

            % dev: store beta
            beta = b(2); betaSum = 0 ;betaSqSum = 0;
            betaEvolve = nan(numRuns-numBurnIn,1);
            betaMeanEvolve = nan(numRuns-numBurnIn,1);

            % dev: store sigmaV
            % sigmaV = vol of vol process (variance)
            C = distVal.C; c = distVal.c;
            sigmaV = C/(c - 2); sigmaVSum = 0; sigmaVSqSum = 0;
            sigmaVEvolve = nan(numRuns-numBurnIn,1);
            sigmaVMeanEvolve = nan(numRuns-numBurnIn,1);

            % dev: store rho
            % rho = corr between W1 and W2 process used for return process (Y)
            % and vol process (variance)
            rho = 0; rhoSum = 0; rhoSqSum = 0;
            rhoEvolve = nan(numRuns-numBurnIn,1);
            rhoMeanEvolve = nan(numRuns-numBurnIn,1);

            % dev: store muV_J
            % muV_J = exp param for volatility (process variance) jump distribution
            D = distVal.D; d = distVal.d;
            muV_J = D/(d-2); muV_JSum = 0; muV_JSqSum = 0;
            muV_JEvolve = nan(numRuns-numBurnIn,1);
            muV_JMeanEvolve = nan(numRuns-numBurnIn,1);

            % dev: store muY_J
            % muY_J = mean of jump size in return process (Y)
            e = distVal.e;
            muY_J = e; muY_JSum = 0; muY_JSqSum = 0;
            muY_JEvolve = nan(numRuns-numBurnIn,1);
            muY_JMeanEvolve = nan(numRuns-numBurnIn,1);


            % dev: store sigmaY_J
            % sigmaY_J = vol of jump in return (process Y)
            F = distVal.F; f = distVal.f;
            sigmaY_J = F/(f - 2); sigmaY_JSum = 0; sigmaY_JSqSum = 0;
            sigmaY_JEvolve = nan(numRuns-numBurnIn,1);
            sigmaY_JMeanEvolve = nan(numRuns-numBurnIn,1);

            % dev: store rhoJ
            % rhoJ = corr between return jump process and vol jump processXVsum
            g = distVal.g;
            rhoJ = g; rhoJSum = 0; rhoJSqSum = 0;
            rhoJEvolve = nan(numRuns-numBurnIn,1);
            rhoJMeanEvolve = nan(numRuns-numBurnIn,1);

            % dev: store lambda
            % lambda = jump intensity
            lambda = 0; lambdaSum = 0; lambdaSqSum = 0;
            lambdaEvolve = nan(numRuns-numBurnIn,1);
            lambdaMeanEvolve = nan(numRuns-numBurnIn,1);

            % get returns as ret
            ret = this.returnTS.getValues;
            numOfRet = length(ret);

            % initial values for variance_t;
            variance = 0.1*(ret - mean(ret)).^2 + 0.9*var(ret);
            varianceSum = 0; varianceSqSum = 0;

            % estimate jumps
            J = abs(ret) - mean(ret) > 2 * std(ret);
            JSum = 0;

            % get jump size in volatility (variance), zV
            zV = exprnd(muV_J,numOfRet,1);
            zVSum = 0;

            % the jump size in return (Y), zY
            zY = mvnrnd((muY_J + zV * rhoJ), sigmaY_J);
            zYSum = 0;

            % ===== CHECK ref article WHERE THEY ARE COMING FROM
            % moments of prioire non standard distribution
            stdevrho = 0.2666; %0.01;
            dfrho = 8; %6/5;
            stdevV = 0.9;
            dfV = 4.5;
            acceptsumV = zeros(size(variance));
            acceptsumrho = 0;
            acceptsums2V = 0;
            Z = ones(length(ret),1);
            test = zeros(length(ret),10);%%matrix for params

            Logger.getInstance.log(LogType.INFO,...
                'Calibration routine is estimating parameters...');
            for i = 1:numRuns
                Rho = 1 / (1 - rho^2);
                variance0 = variance(1);

                % Draw muY(i+1)
                Q = ( ret - zY.*J - rho/sigmaV^0.5.*...
                    ( variance - [variance0;variance(1:end-1)]*(1+beta)-alpha - J.*zV ) )./ ([variance0;variance(1:end-1)] ).^0.5;
                %W = [1./([variance(1);variance(1:end-1)]).^0.5, [0;ret(1:end-1)]./([variance(1);variance(1:end-1)]).^0.5];
                W = 1./([variance0;variance(1:end-1)]).^0.5;

                As = inv( inv(distVal.A) + 1 / (1-rho^2) * W'*W );
                as = As*( inv(distVal.A) * distVal.a + 1 / ( 1 - rho^2 ) * W'*Q );

                muY = mvnrnd(as',As)';

                if i > numBurnIn
                    muYSum = muYSum + muY;
                    muYEvolve(i-numBurnIn) = muYSum/(i-numBurnIn); % for trace
                    %muYMeanEvolve(i-numBurnIn) = sum(muYEvolve(1:i-numBurnIn))/(i-numBurnIn); % for trace
                    muYSqSum = muYSqSum + muY.^2;
                end

                % (alpha, beta)
                eY = ret - Z*muY - zY.*J; %expected return
                eV = variance - [variance0;variance(1:end-1)] - zV.*J; %expected variance
                Q = (eV - rho*sqrt(sigmaV)*eY)./[variance0;variance(1:end-1)].^0.5;
                W = [1./[variance0;variance(1:end-1)].^0.5 [variance0;variance(1:end-1)].^0.5];
                Bs = inv(inv(distVal.B) + (Rho/sigmaV)*W'*W);
                bs = Bs*(inv(distVal.B)*distVal.b + (Rho/sigmaV)*(W'*Q));
                temp = mvnrnd(bs,Bs);
                alpha = temp(1);
                beta = temp(2);
                kappa = - alpha/beta;
                if i > numBurnIn
                    alphaSum = alphaSum + alpha;alphaSqSum = alphaSqSum + alpha^2;
                    alphaEvolve(i-numBurnIn) = alphaSum/(i-numBurnIn); % for trace
                    %alphaMeanEvolve(i-numBurnIn) = sum(alphaEvolve(1:i-numBurnIn))/(i-numBurnIn); % for trace

                    betaSum = betaSum + beta; betaSqSum = betaSqSum + beta^2;
                    betaEvolve(i-numBurnIn) = betaSum/(i-numBurnIn); % for trace
                    %betaMeanEvolve(i-numBurnIn) = sum(betaEvolve(1:i-numBurnIn))/(i-numBurnIn); % for trace

                    kappaSum = kappaSum + kappa; kappaSqSum = kappaSqSum + kappa^2;
                    kappaEvolve(i-numBurnIn) = kappaSum/(i-numBurnIn); % for trace
                    %kappaMeanEvolve(i-numBurnIn) = sum(kappaEvolve(1:i-numBurnIn))/(i-numBurnIn); % for trace
                end

                % sigmaV
                cs = c + length(ret);
                Cs = C + sum( ((variance - [variance0;variance(1:end-1)] - alpha - beta*[variance0;variance(1:end-1)] - zV.*J).^2)./[variance0;variance(1:end-1)] );
                sigmaVTemp = iwishrnd(Cs,cs);

                q = exp(-0.5*sum((variance - [variance0;variance(1:end-1)]*(1+beta) - alpha - J.*zV).^2./(sigmaVTemp*[variance0;variance(1:end-1)])-...
                    (variance - [variance0;variance(1:end-1)]*(1+beta) - alpha - J.*zV).^2./(sigmaV*[variance0;variance(1:end-1)])));

                p = exp(-0.5*sum((variance - [variance0;variance(1:end-1)]*(1+beta) - alpha - J.*zV - rho*sigmaVTemp^0.5*(ret-Z*muY-J.*zY)).^2./...
                    ((1-rho^2)*sigmaVTemp*[variance0;variance(1:end-1)])-...
                    (variance - [variance0;variance(1:end-1)]*(1+beta) - alpha - J.*zV -rho*sigmaV^0.5*(ret-Z*muY-J.*zY)).^2./...
                    ((1-rho^2)*sigmaV*[variance0;variance(1:end-1)])));

                x = min(p/q,1);
                u = rand(1);
                if x > u
                    sigmaV = sigmaVTemp;
                    if i > numBurnIn; acceptsums2V = acceptsums2V + 1;end
                end

                if i > numBurnIn
                    sigmaVSum = sigmaVSum + sigmaV;
                    sigmaVEvolve(i-numBurnIn) = sigmaVSum/(i-numBurnIn); % for trace
                    %sigmaVMeanEvolve(i-numBurnIn) = sum(sigmaVEvolve(1:i-numBurnIn))/(i-numBurnIn); % for trace
                    sigmaVSqSum = sigmaVSqSum + sigmaV^2;
                end

                % rho
                %Draw a candidate for rho(i+1)
                rhoprop = rho + stdevrho * trnd(dfrho); %draw rhoc from a t distribution with 8 df and std  of 0.2666
                if abs(rhoprop) < 1
                    p = (sqrt( 1 - rho^2 )/ sqrt( 1 - rhoprop^2 )).^length(ret) * exp( sum( - 1 / ( 2 * ( 1 - rhoprop^2 ) ) *...
                        ( ret - Z*muY - J.*zY - rhoprop / sigmaV^0.5 * ( variance - alpha - [variance0;variance(1:end-1)] * ( 1 + beta ) - J.*zV ) ).^2./...
                        [variance0;variance(1:end-1)] + 1 / ( 2 * ( 1 - rho^2 ) ) *...
                        ( ret - Z*muY - J.*zY - rho / sigmaV^0.5 * ( variance(1:end) - alpha - [variance0;variance(1:end-1)] * ( 1 + beta ) - J.*zV ) ).^2./...
                        [variance0;variance(1:end-1)] ) );
                    u = rand(1);
                    x = min(p,1);
                    if x > u
                        rho = rhoprop;
                        if i > numBurnIn; acceptsumrho = acceptsumrho +1;end
                    end
                end
                if i > numBurnIn
                    rhoSum = rhoSum + rho;
                    rhoEvolve(i-numBurnIn) = rhoSum/(i-numBurnIn); % for trace
                    %rhoMeanEvolve(i-numBurnIn) = sum(rhoEvolve(1:i-numBurnIn))/(i-numBurnIn); % for trace
                    rhoSqSum = rhoSqSum + rho^2;
                end

                % muV_J
                ds = d + 2*length(ret);
                Ds = D + 2*sum(zV);
                muV_J = iwishrnd(Ds,ds);
                if i > numBurnIn
                    muV_JSum = muV_JSum + muV_J;
                    muV_JEvolve(i-numBurnIn) = muV_JSum/(i-numBurnIn); % for trace
                    %muV_JMeanEvolve(i-numBurnIn) = sum(muV_JEvolve(1:i-numBurnIn))/(i-numBurnIn); % for trace
                    muV_JSqSum = muV_JSqSum + muV_J^2;
                end

                % muY_J
                Es = 1/(length(ret)/sigmaY_J + 1/distVal.E);
                es = Es * (sum( (zY - zV*rhoJ)/sigmaY_J ) + distVal.e/distVal.E);
                muY_J = normrnd(es,Es^0.5);
                if i > numBurnIn
                    muY_JSum = muY_JSum + muY_J;
                    muY_JEvolve(i-numBurnIn) = muY_JSum/(i-numBurnIn); % for trace
                    %muY_JMeanEvolve(i-numBurnIn) = sum(muY_JEvolve(1:i-numBurnIn))/(i-numBurnIn); % for trace
                    muY_JSqSum = muY_JSqSum + muY_J^2;
                end

                % s2YSqSum
                fs = f + length(ret);
                Fs = F + sum((zY - muY_J - rhoJ*zV).^2);
                sigmaY_J = iwishrnd(Fs,fs);
                if i > numBurnIn
                    sigmaY_JSum = sigmaY_JSum + sigmaY_J;
                    sigmaY_JEvolve(i-numBurnIn) = sigmaY_JSum/(i-numBurnIn); % for trace
                    %sigmaY_JMeanEvolve(i-numBurnIn) = sum(sigmaY_JEvolve(1:i-numBurnIn))/(i-numBurnIn); % for trace
                    sigmaY_JSqSum = sigmaY_JSqSum + sigmaY_J^2;
                end

                % rhoJ
                Gs = inv(sum(zV.^2)/sigmaY_J + 1/distVal.G);
                gs = Gs * (sum((zY - muY_J).*zV)/sigmaY_J + distVal.g/distVal.G);
                rhoJ = normrnd(gs,Gs^0.5);
                if i > numBurnIn
                    rhoJSum = rhoJSum + rhoJ; rhoJSqSum = rhoJSqSum + rhoJ^2;
                    rhoJEvolve(i-numBurnIn) = rhoJSum/(i-numBurnIn); % for trace
                    %rhoJMeanEvolve(i-numBurnIn) = sum(rhoJEvolve(1:i-numBurnIn))/(i-numBurnIn); % for trace
                end

                % lambda
                ks = distVal.k + sum(J);
                Ks = distVal.K + length(ret) - sum(J);
                lambda = betarnd(ks,Ks);
                if i > numBurnIn
                    lambdaSum = lambdaSum + lambda;
                    lambdaEvolve(i-numBurnIn) = lambdaSum/(i-numBurnIn); % for trace
                    %lambdaMeanEvolve(i-numBurnIn) = sum(lambdaEvolve(1:i-numBurnIn))/(i-numBurnIn); % for trace % for trace
                    lambdaSqSum = lambdaSqSum + lambda^2;
                end

                % J
                eY1 = ret - Z*muY - zY;
                eY2 = ret - Z*muY;
                eV1 = variance - [variance0;variance(1:end-1)] - alpha - beta*[variance0;variance(1:end-1)] - zV;
                eV2 = variance - [variance0;variance(1:end-1)] - alpha - beta*[variance0;variance(1:end-1)];
                p1 = lambda*exp( -0.5 * ( ((eY1 - (rho/sqrt(sigmaV))*eV1).^2)./((1-rho^2)*[variance0;variance(1:end-1)]) + (eV1.^2)./(sigmaV*[variance0;variance(1:end-1)]) ) );
                p2 = (1 - lambda) * exp( -0.5 * ( ((eY2 - (rho/sqrt(sigmaV))*eV2).^2)./((1-rho^2)*[variance0;variance(1:end-1)]) + (eV2.^2)./(sigmaV*[variance0;variance(1:end-1)]) ) );
                p = p1./(p1 + p2);
                tmp=[tmp;[p1 p2 p]];

                u = rand(length(ret),1);
                J = double(u < p);
                if i > numBurnIn; JSum = JSum + J; end

                Jindex = find(J == 1);

                % zV
                zV(logical(~J)) = exprnd(muV_J,length(ret)-sum(J),1);
                if ~isempty(Jindex)
                    if Jindex(1) == 1
                        t = 1;
                        eV = variance(1) - variance0 - alpha - beta*variance(1);
                        eY = ret(1) - Z(1,:)*muY - zY(1);
                        H = inv( 1 /((1 - rho^2)*sigmaV*variance0) + rhoJ^2/sigmaY_J );
                        h = H * ((eV-rho*sqrt(sigmaV)*eY)/((1 - rho^2)*sigmaV*variance0) + rhoJ*(zY(1) - muY_J)/sigmaY_J - 1/muV_J);
                        if h+5*sqrt(H) > 0; zV(1) = normt_rnd(h,H,0,h+5*sqrt(H)); else zV(1) = 0; end
                        if zV(1) == Inf || isnan(zV(1)) ; zV(1) = 0; end
                    else
                        t = Jindex(1);
                        eV = variance(t) - variance(t-1) - alpha - beta*variance(t-1);
                        eY = ret(t) - Z(t,:)*muY - zY(t);
                        H = inv( 1 /((1 - rho^2)*sigmaV*variance(t-1)) + rhoJ^2/sigmaY_J );
                        h = H * ((eV-rho*sqrt(sigmaV)*eY)/((1 - rho^2)*sigmaV*variance(t-1)) + rhoJ*(zY(t) - muY_J)/sigmaY_J - 1/muV_J);
                        if h+5*sqrt(H) > 0; zV(t) = normt_rnd(h,H,0,h+5*sqrt(H)); else zV(t) = 0; end
                        if zV(t) == Inf || isnan(zV(t)); zV(t) = 0; end
                    end
                    if length(Jindex) > 1
                        for t = Jindex(2:end)'
                            eV = variance(t) - variance(t-1) - alpha - beta*variance(t-1);
                            eY = ret(t) - Z(t,:)*muY - zY(t);
                            H = inv( 1 /((1 - rho^2)*sigmaV*variance(t-1)) + rhoJ^2/sigmaY_J );
                            h = H * ((eV-rho*sqrt(sigmaV)*eY)/((1 - rho^2)*sigmaV*variance(t-1)) + rhoJ*(zY(t) - muY_J)/sigmaY_J - 1/muV_J);
                            if h+5*sqrt(H) > 0; zV(t) = normt_rnd(h,H,0,h+5*sqrt(H)); else zV(t) = 0; end
                            if zV(t) == Inf || isnan(zV(t)); zV(t) = 0; end
                        end
                    end
                end
                if i > numBurnIn; zVSum = zVSum + zV; end

                % zY
                zY(logical(~J)) = normrnd(muY_J + rhoJ*zV(logical(~J)),sigmaY_J^0.5)';
                if ~isempty(Jindex)
                    if Jindex(1) == 1
                        t = 1;
                        eV = variance(1) - variance0 - alpha - beta*variance(1) - zV(1);
                        eY = ret(1) - Z(1,:)*muY;
                        L = inv(1/((1 - rho^2)*variance0) + 1/sigmaY_J);
                        l = L * ( (eY - (rho/sqrt(sigmaV))*eV)/((1 - rho^2)*variance0) + (muY_J + rhoJ*zV(1))/sigmaY_J );
                        zY(1) = normrnd(l,sqrt(L));
                    else
                        t = Jindex(1);
                        eV = variance(t) - variance(t-1) - alpha - beta*variance(t-1) - zV(t);
                        eY = ret(t) - Z(t,:)*muY;
                        L = inv(1/((1 - rho^2)*variance(t-1)) + 1/sigmaY_J);
                        l = L * ( (eY - (rho/sqrt(sigmaV))*eV)/((1 - rho^2)*variance(t-1)) + (muY_J + rhoJ*zV(t))/sigmaY_J );
                        zY(t) = normrnd(l,sqrt(L));
                    end
                    if length(Jindex) > 1
                        for t = Jindex(2:end)'
                            eV = variance(t) - variance(t-1) - alpha - beta*variance(t-1) - zV(t);
                            eY = ret(t) - Z(t,:)*muY;
                            L = inv(1/((1 - rho^2)*variance(t-1)) + 1/sigmaY_J);
                            l = L * ( (eY - (rho/sqrt(sigmaV))*eV)/((1 - rho^2)*variance(t-1)) + (muY_J + rhoJ*zV(t))/sigmaY_J );
                            zY(t) = normrnd(l,sqrt(L));
                        end
                    end
                end
                if i > numBurnIn; zYSum = zYSum + zY;end

                % Draw variance
                epsilon = trnd(dfV,length(ret),1);
                [mv, v] = tstat(dfV);
                epsilon = (stdevV/sqrt(v)) * epsilon;
                if i == floor(numBurnIn / 2)
                    Vindex1 = find(varianceSqSum > prctile(varianceSqSum, 92.5));
                    Vindex2 = find(varianceSqSum > prctile(varianceSqSum, 75) & ...
                        varianceSqSum < prctile(varianceSqSum, 92.5));
                    Vindex3 = find(varianceSqSum < prctile(varianceSqSum, 25) & ...
                        varianceSqSum > prctile(varianceSqSum, 2.5));
                    Vindex4 = find(varianceSqSum < prctile(varianceSqSum, 2.5));
                end
                if i > floor(numBurnIn / 2) - 1
                    epsilon(Vindex1) = 1.35 * epsilon(Vindex1);
                    epsilon(Vindex2) = 1.25 * epsilon(Vindex2);
                    epsilon(Vindex3) = 0.75 * epsilon(Vindex3);
                    epsilon(Vindex4) = 0.65 * epsilon(Vindex4);
                end
                j = 1;
                Vprop = variance + epsilon;
                p1 = max(0,exp( -0.5 * ( ( ret(j+1) - Z(j+1,1)*muY(1)  - J(j+1)*zY(j+1) - rho / sigmaV^0.5 *(variance(j+1) - Vprop(j) - alpha - Vprop(j) * beta - J(j+1)*zV(j+1) ) )^2/( (1 - rho^2) * Vprop(j) ) +...
                    ( ret(j) - Z(j,1)*muY(1)  - J(j)*zY(j) - rho / sigmaV^0.5 *(Vprop(j) - variance(1) - alpha - variance0 * beta - J(j)*zV(j)))^2/( (1 - rho^2) * variance0 ) +...
                    ( variance(j+1) - Vprop(j) - alpha - Vprop(j) * beta - J(j+1)*zV(j+1) )^2/( sigmaV * Vprop(j) ) +...
                    ( Vprop(j) - variance0 - alpha - variance0 * beta - J(j)*zV(j))^2/( sigmaV * variance0 ) ) ) / Vprop(j));
                p2 = max(0,exp( -0.5 * ( ( ret(j+1) - Z(j+1,1)*muY(1)  - J(j+1)*zY(j+1) - rho / sigmaV^0.5 *(variance(j+1) - variance(j) - alpha - variance(j) * beta - J(j+1)*zV(j+1) ) )^2/( (1 - rho^2) * variance(j) ) +...
                    ( ret(j) - Z(j,1)*muY(1)  - J(j)*zY(j) -rho / sigmaV^0.5 *(variance(j) - variance0 - alpha - variance0 * beta - J(j)*zV(j)) )^2/( (1 - rho^2) * variance0 ) +...
                    ( variance(j+1) - variance(j) - alpha - variance(j) * beta - J(j+1)*zV(j+1))^2/( sigmaV * variance(j) ) +...
                    ( variance(j) - variance0 - alpha - variance0 * beta - J(j)*zV(j))^2/( sigmaV * variance0 ) ) ) / variance(j));
                if p2 ~= 0; acceptV = min(p1/p2, 1); elseif p1 > 0; acceptV = 1; else; acceptV = 0; end
                u = rand(length(ret),1);
                if u(j) < acceptV
                    variance(j) = Vprop(j);
                    if i > numBurnIn; acceptsumV(j) = acceptsumV(j) + 1; end
                end

                for j = 2:length(ret)-1
                    p1 = max(0,exp( -0.5 * ( ( ret(j+1) - Z(j+1,1)*muY(1)  - J(j+1)*zY(j+1) - rho / sigmaV^0.5 *(variance(j+1) - Vprop(j) - alpha - Vprop(j) * beta - J(j+1)*zV(j+1)) )^2/( (1 - rho^2) * Vprop(j) ) +...
                        ( ret(j) - Z(j,1)*muY(1)  - J(j)*zY(j) - rho / sigmaV^0.5 *(Vprop(j) - variance(j-1) - alpha - variance(j-1) * beta - J(j)*zV(j) ) )^2/( (1 - rho^2) * variance(j-1) ) +...
                        ( variance(j+1) - Vprop(j) - alpha - Vprop(j) * beta - J(j+1)*zV(j+1))^2/( sigmaV * Vprop(j) ) +...
                        ( Vprop(j) - variance(j-1) - alpha - variance(j-1) * beta - J(j)*zV(j))^2/( sigmaV * variance(j-1) ) ) ) / Vprop(j));
                    p2 = max(0,exp( -0.5 * ( ( ret(j+1) - Z(j+1,1)*muY(1)  - J(j+1)*zY(j+1) - rho / sigmaV^0.5 *(variance(j+1) - variance(j) - alpha - variance(j) * beta - J(j+1)*zV(j+1)) )^2/( (1 - rho^2) * variance(j) ) +...
                        ( ret(j) - Z(j,1)*muY(1) - J(j)*zY(j) - rho / sigmaV^0.5 *(variance(j) - variance(j-1) - alpha - J(j)*zV(j) - variance(j-1) * beta ) )^2/( (1 - rho^2) * variance(j-1) ) +...
                        ( variance(j+1) - variance(j) - alpha - variance(j) * beta - J(j+1)*zV(j+1))^2/( sigmaV * variance(j) ) +...
                        ( variance(j) - variance(j-1) - alpha - variance(j-1) * beta - J(j)*zV(j))^2/( sigmaV * variance(j-1) ) ) ) / variance(j));
                    if p2 ~= 0; acceptV = min(p1/p2, 1); elseif p1 > 0; acceptV = 1; else; acceptV = 0; end

                    if u(j) < acceptV
                        variance(j) = Vprop(j);
                        if i > numBurnIn; acceptsumV(j) = acceptsumV(j) + 1; end
                    end
                end
                j = length(ret);
                p1 = max(0,exp( -0.5 * ( ( ret(j) - Z(j,1)*muY(1)  - J(j)*zY(j) - rho / sigmaV^0.5 *(Vprop(j) - variance(j-1) - alpha - variance(j-1) * beta - J(j)*zV(j)) )^2/( (1 - rho^2) * variance(j-1) ) +...
                    ( Vprop(j) - variance(j-1) - alpha - variance(j-1) * beta - J(j)*zV(j))^2/( sigmaV * variance(j-1) ) ) ) / Vprop(j)^0.5);
                p2 = max(0,exp( -0.5 * ( ( ret(j) - Z(j,1)*muY(1)  - J(j)*zY(j) - rho / sigmaV^0.5 *(variance(j) - variance(j-1) - alpha - variance(j-1) * beta - J(j)*zV(j)) )^2/( (1 - rho^2) * variance(j-1) ) +...
                    ( variance(j) - variance(j-1) - alpha - variance(j-1) * beta  - J(j)*zV(j))^2/( sigmaV * variance(j-1) ) ) ) / variance(j)^0.5);
                if p2 ~= 0; acceptV = min(p1/p2, 1); elseif p1 > 0; acceptV = 1; else; acceptV = 0; end

                if u(j) < acceptV
                    variance(j) = Vprop(j);
                    if i > numBurnIn; acceptsumV(j) = acceptsumV(j) + 1; end
                end

                if i > numBurnIn; varianceSum = varianceSum + variance;end
                if i > floor(numBurnIn / 2) - 100 || i < floor(numBurnIn / 2); varianceSqSum = varianceSqSum + variance; end
                test(i,:) = [muY muY_J sigmaY_J lambda alpha beta  rho sigmaV rhoJ muV_J];

            end
            Logger.getInstance.log(LogType.INFO,...
                'Model calibration completed!');
            sampleSize = (numRuns-numBurnIn);
            movingWindow = 5;
            % mean of return process
            this.calibratedParameters.muY = muYSum/sampleSize;
            this.traceParamSet.muY.evolve = muYEvolve;
            this.traceParamSet.muY.meanEvolve = movmean(muYEvolve,movingWindow);

            % mean on return jump process
            this.calibratedParameters.muY_J = muY_JSum/sampleSize;
            this.traceParamSet.muY_J.evolve = muY_JEvolve;
            this.traceParamSet.muY_J.meanEvolve = movmean(muY_JEvolve,movingWindow);

            % vol of return jump process
            this.calibratedParameters.sigmaY_J = sigmaY_JSum/sampleSize;
            this.traceParamSet.sigmaY_J.evolve = sigmaY_JEvolve;
            this.traceParamSet.sigmaY_J.meanEvolve = movmean(sigmaY_JEvolve,movingWindow);

            % alpha
            this.calibratedParameters.alpha = alphaSum/sampleSize;
            this.traceParamSet.alpha.evolve = alphaEvolve;
            this.traceParamSet.alpha.meanEvolve = movmean(alphaEvolve,movingWindow);

            % beta
            this.calibratedParameters.beta = betaSum/sampleSize;
            this.traceParamSet.beta.evolve = betaEvolve;
            this.traceParamSet.beta.meanEvolve = movmean(betaEvolve,movingWindow);

            % rho
            this.calibratedParameters.rho = rhoSum/sampleSize;
            this.traceParamSet.rho.evolve = rhoEvolve;
            this.traceParamSet.rho.meanEvolve = movmean(rhoEvolve,movingWindow);

            % sigmaV vol of vol
            this.calibratedParameters.sigmaV = sigmaVSum/sampleSize;
            this.traceParamSet.sigmaV.evolve = sigmaVEvolve;
            this.traceParamSet.sigmaV.meanEvolve = movmean(sigmaVEvolve,movingWindow);

            % kappa
            this.calibratedParameters.kappa = kappaSum/sampleSize;
            this.traceParamSet.kappa.evolve = kappaEvolve;
            this.traceParamSet.kappa.meanEvolve = movmean(kappaEvolve,movingWindow);

            % rhoJ
            this.calibratedParameters.rhoJ = rhoJSum/sampleSize;
            this.traceParamSet.rhoJ.evolve = rhoJEvolve;
            this.traceParamSet.rhoJ.meanEvolve = movmean(rhoJEvolve,movingWindow);

            % muV_J
            this.calibratedParameters.muV_J = muV_JSum/sampleSize;
            this.traceParamSet.muV_J.evolve = muV_JEvolve;
            this.traceParamSet.muV_J.meanEvolve = movmean(muV_JEvolve,movingWindow);

            % lambda
            this.calibratedParameters.lambda = lambdaSum/sampleSize;
            this.traceParamSet.lambda.evolve = lambdaEvolve;
            this.traceParamSet.lambda.meanEvolve = movmean(lambdaEvolve,movingWindow);

            % store itercount to trace param for plotting
            this.traceParamSet.iterCount = sampleSize;

            % store err measurement params
            this.errMeasureParam.Jsum = JSum;
            this.errMeasureParam.zYsum = zYSum;
            this.errMeasureParam.zVsum = zVSum;
            this.errMeasureParam.varianceSum = varianceSum;
            this.errMeasureParam.iterCount = sampleSize;
            %             disp(['muY: mean of returns (Y)'])
            %             disp([muYSum/(numRuns-numBurnIn) (muYSqSum/(numRuns-numBurnIn)-(muYSum/(numRuns-numBurnIn)).^2).^0.5])
            %             disp(['muY_J: mean of jump size in return process (Y)'])
            %             disp([muY_JSum/(numRuns-numBurnIn) (muY_JSqSum/(numRuns-numBurnIn)-(muY_JSum/(numRuns-numBurnIn))^2)^0.5])
            %             disp(['sigmaY_J: Vol of jump in return process (Y)'])
            %             disp([sigmaY_JSum/(numRuns-numBurnIn) (sigmaY_JSqSum/(numRuns-numBurnIn)-(sigmaY_JSum/(numRuns-numBurnIn))^2)^0.5])
            %             disp(['lambda'])
            %             disp([lambdaSum/(numRuns-numBurnIn) (lambdaSqSum/(numRuns-numBurnIn)-(lambdaSum/(numRuns-numBurnIn))^2)^0.5])
            %             disp(['alpha'])
            %             disp([alphaSum/(numRuns-numBurnIn) (alphaSqSum/(numRuns-numBurnIn)-(alphaSum/(numRuns-numBurnIn))^2)^0.5])
            %             disp(['beta'])
            %             disp([betaSum/(numRuns-numBurnIn) (betaSqSum/(numRuns-numBurnIn)-(betaSum/(numRuns-numBurnIn))^2)^0.5])
            %             disp(['kappa'])
            %             disp([kappaSum/(numRuns-numBurnIn) (kappaSqSum/(numRuns-numBurnIn)-(kappaSum/(numRuns-numBurnIn))^2)^0.5])
            %             disp(['rho: corr between W1 and W2'])
            %             disp([rhoSum/(numRuns-numBurnIn) (rhoSqSum/(numRuns-numBurnIn)-(rhoSum/(numRuns-numBurnIn))^2)^0.5])
            %             disp(['sigmaV: vol of vol process (V)'])
            %             disp([sigmaVSum/(numRuns-numBurnIn) (sigmaVSqSum/(numRuns-numBurnIn)-(sigmaVSum/(numRuns-numBurnIn))^2)^0.5])
            %             disp(['rhoJ: corr between ZY and ZV, the 2 jump processes'])
            %             disp([rhoJSum/(numRuns-numBurnIn) (rhoJSqSum/(numRuns-numBurnIn)-(rhoJSum/(numRuns-numBurnIn))^2)^0.5])
            %             disp(['muV_J: Exp param for Jump in Vol process'])
            %             disp([muV_JSum/(numRuns-numBurnIn) (muV_JSqSum/(numRuns-numBurnIn)-(muV_JSum/(numRuns-numBurnIn))^2)^0.5])


            %             if (true)
            this.calibrated = ModelCalibration.SUCCESS;
            %             else
            %                 this.isCalibrated = ModelCalibration.FAILURE;
            %             end
        end

        %simulate
        function simulatePrices(this)
            % mean of return process
            muY = this.calibratedParameters.muY;
            % mean on return jump process
            muY_J = this.calibratedParameters.muY_J;
            % vol of return jump process
            sigmaY_J = this.calibratedParameters.sigmaY_J;
            % alpha
            alpha = this.calibratedParameters.alpha;
            % beta
            beta = this.calibratedParameters.beta;
            % rho
            rho = this.calibratedParameters.rho;
            % sigmaV vol of vol
            sigmaV = this.calibratedParameters.sigmaV;
            % kappa
            kappa = this.calibratedParameters.kappa;
            % rhoJ
            rhoJ = this.calibratedParameters.rhoJ;
            % muV_J
            muV_J = this.calibratedParameters.muV_J;
            % lambda
            lambda = this.calibratedParameters.lambda;
            % theta
            theta = alpha/kappa;

            % get simulated confing and load
            simulationConfig = getConfig('SVCJ','simulation');
            numOfIter = simulationConfig.numOfIter;
            numOfPeriod = simulationConfig.numOfPeriod;

            % dev: create empty vectors for process variables
            Y = nan(numOfIter,numOfPeriod); % log return
            V = nan(numOfIter,numOfPeriod); % vol process
            P = nan(numOfIter,numOfPeriod); % simulated Prices
            % set V0 for all parths to mean vol muV_J
            V(:,1) = muV_J;
            % set initalPrice
            iniPrice = this.priceTS.getLastValue;
            P(:,1) = iniPrice;
            Jy = nan(numOfIter,numOfPeriod); % jumps in return process
            Jv = nan(numOfIter,numOfPeriod); % jumps in vol process
            % initial price


            for i = 2:numOfPeriod
                z = mvnrnd([0;0],[1,rho;rho,1]); % multivariate for W1 W2
                z1 = z(:,1); % W1
                z2 = z(:,2); % W2
                rr = rand(1,numOfIter);
                j = (rr<lambda)';  % bernoulli for jump
                XV = exprnd(1/muV_J,[numOfIter,1]); % exp jump distribution for vol process
                X = normrnd(muY_J + rhoJ * XV, sigmaY_J,[numOfIter,1]); % normrandom jump for return process

                V(:,i) = kappa * theta + (1 - kappa) * V(:,i- 1) + sigmaV*sqrt(V(:,i-1))*z2 + XV.*j; % vol process
                %             zz = alpha + beta * V(:,i- 1) + sigmaV*sqrt(V(:,i-1))*z2 + XV.*j;
                Y(:,i) = muY + sqrt(V(:,i-1))*z1 + X.*j;  % return process
                Jv(:,i) = XV.*j;  % Jumps in volatilty (0 in case of no jump)
                Jy(:,i) = X.*j;  % Jumps in log return (0 in case of no jump)
                P(:,i) = (Y(:,i)/sqrt(this.calibrationConfig.annulisationFactor)+1).*P(:,i-1);

            end
            this.simulatedReturns = Y;
            this.simulatedVolatility = V;
            this.simulatedJumpsInReturns = Jy;
            this.simulatedJumpsInVol = Jv;
            this.simulatedPrices = P;

        end

        % checkAndLoadModelConfig
        function checkAndLoadModelConfig(this)
            modelConfig = getConfig('SVCJ','calibration');
            % returnMethod
            if isfield(modelConfig,'returnMethod')
                this.calibrationConfig.returnMethod = modelConfig.returnMethod;
                Logger.getInstance.log(LogType.INFO,...
                    'modelConfig provided, returnMethod set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Incomplete modelConfig provided, returnMethod missing');
            end

            % annulisationFactor
            if isfield(modelConfig,'annulisationFactor')
                this.calibrationConfig.annulisationFactor = modelConfig.annulisationFactor;
                Logger.getInstance.log(LogType.INFO,...
                    'modelConfig provided, annulisationFactor set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Incomplete modelConfig provided, annulisationFactor missing');
            end

            % NumRuns
            if isfield(modelConfig,'NumRuns')
                this.calibrationConfig.NumRuns = modelConfig.NumRuns;
                Logger.getInstance.log(LogType.INFO,...
                    'modelConfig provided, NumRuns set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Incomplete modelConfig provided, NumRuns missing');
            end

            % burnInPeriod
            if isfield(modelConfig,'burnInPeriod')
                this.calibrationConfig.burnInPeriod = modelConfig.burnInPeriod;
                Logger.getInstance.log(LogType.INFO,...
                    'modelConfig provided, burnInPeriod set');
            else
                Logger.getInstance.log(LogType.FATAL,...
                    'Incomplete modelConfig provided, burnInPeriod missing');
            end
        end

        % sourceReturnTimeseries
        function sourceReturnTimeseries(this)
            this.returnTS = this.priceTS.getReturnTimeseries(this.calibrationConfig.returnMethod,'annulisationFactor',...
                this.calibrationConfig.annulisationFactor);
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
        function calibrationConfig = getCalibrationConfig(this)
            calibrationConfig = this.calibrationConfig;
        end

        %setModelConfig
        function setCalibrationConfig(this,calibrationConfig)
            this.calibrationConfig = calibrationConfig;
        end

        %getSimulatedReturns
        function simulatedReturns = getSimulatedReturns(this)
            simulatedReturns = this.simulatedReturns;
        end

        %getSimulatedVolatility
        function simulatedVolatility = getSimulatedVolatility(this)
            simulatedVolatility = this.simulatedVolatility;
        end

        %getSimulatedJumpsInReturns
        function simulatedJumpsInReturns = getSimulatedJumpsInReturns(this)
            simulatedJumpsInReturns = this.simulatedJumpsInReturns;
        end

        %getSimulatedJumpsInVol
        function simulatedJumpsInVol = getSimulatedJumpsInVol(this)
            simulatedJumpsInVol = this.simulatedJumpsInVol;
        end

        %getSimulatedPrices
        function simulatedPrices = getSimulatedPrices(this)
            simulatedPrices = this.simulatedPrices;
        end

        %isCalibrated
        function calibrationStatus = isCalibrated(this)
            calibrationStatus = this.calibrated;
        end

        %getCalibratedParameters
        function calibratedParmeters = getCalibratedParameters(this)
            calibratedParmeters = this.calibratedParameters;
        end

        %getTraceParamSet
        function traceParamSet = getTraceParamSet(this)
            traceParamSet = this.traceParamSet;
        end

        %getErrMeasureparam
        function errMeasureParam = getErrMeasureParam(this)
            errMeasureParam = this.errMeasureParam;
        end

    end
end