function modifiedArray = modifyRepeatingElements(~,array, decimalDigit)
    % Initialize the modified array with the original array
    modifiedArray = array;
    
    % Find unique elements and their counts
    [uniqueElements, ~, indices] = unique(array);
    counts = accumarray(indices, 1);
    
    % Identify the repeating elements
    repeatingElements = uniqueElements(counts > 1);
    
    % Scale factor based on the decimal digit
    scaleFactor = 10^decimalDigit;
    
    % Modify each repeating element
    for i = 1:length(repeatingElements)
        % Find indices of the current repeating element
        idx = find(array == repeatingElements(i));
        
        % Modify elements by adding n times at the specified decimal place
        % and place them back in their original positions
        for j = 1:length(idx)
            modifiedArray(idx(j)) = array(idx(j)) + j/scaleFactor;
        end
    end
end