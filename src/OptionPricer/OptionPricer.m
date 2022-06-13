classdef OptionPricer < handle

    %% Protected properties
    properties (Access = protected)
        valuationDate;
        irCurve;
    end

    %% Constructor
    methods
        function this = OptionPricer(valuationDate)
            this.valuationDate = valuationDate;
        end

    end

    %% Public methods
    methods (Access = public)

        %getValuationDate
        function valuationDate = getValuationDate(this)
            valuationDate = this.valuationDate;
        end

        %setValuationDate
        function setValuationDate(this,valuationDate)
            this.valuationDate = valuationDate;
        end

        %getIrCurve
        function irCurve = getIrCurve(this)
            irCurve = this.irCurve;
        end

        %setIrCurve
        function setIrCurve(this,irCurve)
            this.irCurve = irCurve;
        end
    end

    %% Private methods
    methods (Access = protected)

        % checkAndLoadIRCurve(this,instrument)
        function checkAndLoadIRCurve(this,instrument)
            if isempty(this.irCurve)
                currency = instrument.getBaseCurrency;
                this.irCurve = sourceIRCurve(this, this.valuationDate,currency);
                Logger.getInstance.log(LogType.INFO,...
                    ['IR Curve loaded for currency', char(currency), ' and date: ',char(this.valuationDate)]);

            else
                if strcmp(this.irCurve.getSettle,this.valuationDate) && ...
                        strcmp(char(this.irCurve.getCurrency),char(instrument.getBaseCurrency))

                    Logger.getInstance.log(LogType.INFO,...
                        ['IR Curve already loaded for currency', char(instrument.getBaseCurrency), ' and date: ',char(this.valuationDate)]);

                else
                    currency = instrument.getBaseCurrency;
                    this.irCurve = sourceIRCurve(this, this.valuationDate,currency);
                    Logger.getInstance.log(LogType.INFO,...
                        ['IR Curve loaded for currency', char(currency), ' and date: ',char(this.valuationDate)]);

                end

            end
        end

        % sourceIRCurve
        function irCurve = sourceIRCurve(this, settle,currency)
            irCurve = buildIRCurve(datestr(settle),currency);
        end

    end
end

