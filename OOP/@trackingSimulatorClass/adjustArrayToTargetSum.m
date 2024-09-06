function adjustedArray = adjustArrayToTargetSum(~, inputArray, targetSum)
    % Get the size of the input array
    [n, m] = size(inputArray);

    % Create a copy of the input array to adjust
    adjustedArray = inputArray;

    % Calculate the current sum of the array elements
    currentSum = sum(adjustedArray, 'all');

    % Adjust the elements to meet the target sum
    while currentSum ~= targetSum
        % Generate random linear indices
        linearIndices = randperm(n * m);

        for idx = linearIndices
            % Convert linear index to row and column indices
            [row, col] = ind2sub([n, m], idx);

            if currentSum < targetSum
                % Increment the element and update the sum
                adjustedArray(row, col) = adjustedArray(row, col) + 1;
                currentSum = currentSum + 1;
            elseif currentSum > targetSum
                % Decrement the element and update the sum
                if adjustedArray(row, col) > 1 % Ensure the element doesn't go below 1
                    adjustedArray(row, col) = adjustedArray(row, col) - 1;
                    currentSum = currentSum - 1;
                end
            end

            % Check if the target sum is reached
            if currentSum == targetSum
                break;
            end
        end
    end
end
