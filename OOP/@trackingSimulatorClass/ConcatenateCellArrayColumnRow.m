function finalConcat = ConcatenateCellArrayColumnRow(~,trackLengthsMatrix)
    % DESCRIPTION: 
    % 
    [n, m] = size(trackLengthsMatrix); 
    % Initialize the result cell array
    result = cell(1, m);
    % Concatenate each column
    for j = 1:m
        % Initialize an empty array for the concatenation of this column
        columnConcat = [];
    
        % Concatenate arrays in this column
        for i = 1:n
            columnConcat = [columnConcat, trackLengthsMatrix{i, j}];
        end
    
        % Store the concatenated array in the result cell array
        result{j} = columnConcat;
    end

    % 'result' now contains the concatenated arrays for each column

    % Initialize an empty array for the final concatenation
    finalConcat = [];
    
    % Concatenate each array in the 'result' array
    for j = 1:length(result)
        finalConcat = [finalConcat, result{j}];
    end
end 

