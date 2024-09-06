function adjustedArray = adjustArraySum(~,initialArray, targetSum)
    % adjustArraySum adjusts the sum of an nxm array to match a specified value,
    % ensuring no element falls below 0.

    % Make a copy of the initial array
    adjustedArray = initialArray;

    % Calculate the current sum of the array
    currentSum = sum(adjustedArray, 'all');

    % Loop until the sum of the array elements matches the target sum
    while currentSum ~= targetSum
        % Randomly select an index to modify
        [numRows, numCols] = size(adjustedArray);
        randRow = randi(numRows);
        randCol = randi(numCols);

        % Determine if we should add or subtract
        if (currentSum > targetSum && adjustedArray(randRow, randCol) > 0) || currentSum < targetSum
            if rand() < 0.5 && currentSum > targetSum
                % Subtract 1, but only if the element is greater than 0
                adjustedArray(randRow, randCol) = adjustedArray(randRow, randCol) - 1;
            else
                % Add 1
                adjustedArray(randRow, randCol) = adjustedArray(randRow, randCol) + 1;
            end
        end

        % Update the current sum
        currentSum = sum(adjustedArray, 'all');
    end
end
