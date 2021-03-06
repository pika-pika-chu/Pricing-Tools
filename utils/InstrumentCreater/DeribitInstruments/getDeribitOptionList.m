function list = getDeribitOptionList(folderLocation,valuationDate,bidAsk)


futurePriceLoc = [folderLocation,'\FuturesPrices.txt'];
futsPriceTable = readtable(futurePriceLoc,"FileType","text","Delimiter",",");


files = dir(folderLocation);
id = 1;

for i = 1: length(files)
    file = files(i);

    if ~(strcmp(file.name,'.') || strcmp(file.name,'..') || strcmp(file.name,'FuturesPrices.txt'))
        % read table
        tt = readtable([file.folder,'\',file.name]);

        tt= removevars(tt,{'Volume'});

        % remove Nans and remove -
        tt = tt(~any(ismissing(tt),2),:);
        tt = tt(~any(ismissing(tt,'-'),2),:);

        % replace - in instrument column to space
        tt.Instrument = replace(tt.Instrument,"-"," ");

        % convert Instrument to mat
        tt.Instrument = char(tt.Instrument);

        for k = 1:height(tt)
            % find inx of spaces
            idx = find(isspace(tt.Instrument(k,:)));

            % split instrument after each space into new columns
            Ticker{k,1} = tt.Instrument(k,1:idx(1)-1);
            Maturity{k,1} = tt.Instrument(k,idx(1)+1:idx(2)-1);
            Strike{k,1} = tt.Instrument(k,idx(2)+1:idx(3)-1);
            Type{k,1} = tt.Instrument(k,idx(3)+1:end);
        end

        tt.Ticker = Ticker;
        tt.Maturity = Maturity;
        tt.Strike = Strike;
        tt.Type = Type;

        clear Ticker
        clear Maturity
        clear Strike
        clear Type

        % remove puts
        tt = tt(~any(ismissing(tt.Type,'P'),2),:);
        tt = tt(~any(ismissing(tt.Type,' P'),2),:);
        tt = tt(~any(ismissing(tt.Type,'P '),2),:);

        % check if cell str and convert to double
        if iscellstr(tt.IV_Bid_)
            tt.IV_Bid_ = str2double(tt.IV_Bid_);
        end

        % check if cell str and convert to double
        if iscellstr(tt.IV_Ask_)
            tt.IV_Ask_ = str2double(tt.IV_Ask_);
        end

        % check if cell str and convert to double
        if iscellstr(tt.Bid)
            tt.Bid = str2double(tt.Bid);
        end

        % check if cell str and convert to double
        if iscellstr(tt.Ask)
            tt.Ask = str2double(tt.Ask);
        end

        % convert exp to datestr
        tt.Maturity = datestr(tt.Maturity);

        tt = sortrows(tt,"Strike");


        for j = 1:height(tt)

            % get Future Price
            [~,index] = ismember(tt.Maturity(j,:),futsPriceTable.date);
            if ~isempty(index)
                sInitial = futsPriceTable.sInitial(index);

            else
                Logger.getInstance.log(LogType.FATAL,...
                    ['Future price not found for date: ',tt.Maturity(j,:)]);
            end
            strike = str2double(tt.Strike(j,:));
            if strike >=  sInitial
                euOption = EuropeanOption();
                euOption.setId(id);
                euOption.setTicker(tt.Instrument(j,:));
                euOption.setType('call');
                euOption.setName(tt.Instrument(j,:));
                euOption.setMaturity((tt.Maturity(j,:)));
                euOption.setStrike(str2double(tt.Strike(j,:)));
                euOption.setBaseCurrency(Currency.USD);

                % calculte IV_MID
                askVol = tt.IV_Ask_(j)/100;
                bidVol = tt.IV_Bid_(j)/100;
                

                if strcmp(bidAsk,'ask')
                    euOption.setImpliedVol(askVol);
                    euOption.setPrice(askVol * sInitial);

                elseif strcmp(bidAsk,'bid')
                    euOption.setImpliedVol(bidVol);
                    euOption.setPrice(bidVol * sInitial);

                elseif strcmp(bidAsk,'mid')
                    midVol = (askVol + bidVol)/2;
                    euOption.setImpliedVol(midVol);
                    euOption.setPrice(midVol * sInitial);
                end

                % build Coin
                coinName = tt.Ticker{j};
                coin = Coin();
                coin.setBaseCurrency(Currency.USD);
                coin.setId(char(id));
                coin.setName(coinName);
                coin.setTicker(coinName);
                coin.setType(CoinType.COIN);
                coin.setDivYield(0);

                % build Future
                futsMaturity = tt.Maturity(j,:);
                futsPrice = sInitial;
                future = Future();
                future.setUnderlying(coin);
                future.setMaturity(futsMaturity);
                future.setPrice(futsPrice);
                future.setValuationDate(valuationDate)
                future.setBaseCurrency(coin.getBaseCurrency);
                euOption.setUnderlying(future);


                list(id) = euOption;
                id = id + 1;
            end
        end
    end
end

