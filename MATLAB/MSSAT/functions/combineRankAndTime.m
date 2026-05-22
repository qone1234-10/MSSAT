function T = combineRankAndTime(jrank, time)

    fields = fieldnames(jrank); 
    numFields = numel(fields);   

    numRanks = numel(jrank.(fields{1}));  

    dataMatrix = zeros(numRanks + 1, numFields);  

    for i = 1:numFields
        fieldName = fields{i};
        
        dataMatrix(1:numRanks, i) = jrank.(fieldName);  
        
        dataMatrix(numRanks + 1, i) = time.(fieldName);
    end

    rowNames = [arrayfun(@(x) sprintf('# rank %d', x), 1:numRanks, 'UniformOutput', false), {'time'}];
    T = array2table(dataMatrix, 'RowNames', rowNames, 'VariableNames', fields);

end
