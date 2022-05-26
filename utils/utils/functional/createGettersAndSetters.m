function createGettersAndSetters(paramName)
% @dev: function prints getters and setters for class propeties

% @dev: check if input is char
if ~ischar(paramName)
    paramName = char(paramName);
end
    
% @dev: capatilise 1st char of input
ParamName = [upper(paramName(1)),paramName(2:end)];

%% Print getters
disp(['%get',ParamName])
disp (['function ', paramName ,' = get',ParamName, '(this)']);
disp ([paramName , ' = this.' ,paramName,';']);
disp (['end']);


%% Print setters
disp('   ')
disp(['%set',ParamName])
disp (['function set',ParamName, '(this,',paramName,')']);
disp (['this.' paramName ,' = ' paramName,';']);
disp (['end']);
end

