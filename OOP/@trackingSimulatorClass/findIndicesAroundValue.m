function [higherIndex, lowerIndex] = findIndicesAroundValue(~, orderedArray, value)
    % Initialize indices to NaN to indicate 'not found'
    higherIndex = NaN;
    lowerIndex = NaN;

    % Check if the value is higher than the maximum of the array
    if value >= max(orderedArray)
        lowerIndex = 1; % Index of the highest element
        return;
    end

    % Check if the value is lower than the minimum of the array
    if value <= min(orderedArray)
        higherIndex = numel(orderedArray); % Index of the lowest element
        return;
    end

    % Find the index of the smallest value that is greater than the specified value
    higherIndex = find(orderedArray > value, 1, 'last');

    % Find the index of the largest value that is less than the specified value
    lowerIndex = find(orderedArray < value, 1);
end