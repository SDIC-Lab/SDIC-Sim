function adjustedArray = adjustArrayLength(~, originalArray, targetLength)
    % Calculate the current length of the array
    currentLength = length(originalArray);
    
    % Calculate how many zeros need to be added
    zerosToAdd = targetLength - currentLength;
    
    % Check if zeros need to be added
    if zerosToAdd > 0
        % Create an array of zeros
        zerosArray = zeros(1, zerosToAdd);
        
        % Append the zeros array to the original array
        adjustedArray = [originalArray, zerosArray];
    else
        % If the array is already the target length or longer, do not modify it
        adjustedArray = originalArray;
    end
end
