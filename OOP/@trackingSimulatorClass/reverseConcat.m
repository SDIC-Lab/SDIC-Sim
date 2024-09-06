function trackLengthsMatrix = reverseConcat(obj, lengths, n, m)
    % Initialize the cell array with the original dimensions
    trackLengthsMatrix = cell(n, m);

    % Index to keep track of position in 'finalConcat'
    currentIndex = 1;
    finalConcat = obj.energyEmissionMap{end}.energy;

    % Iterate over each cell position
    for j = 1:m
        for i = 1:n
            % Length of the current array to extract
            currentLength = lengths(i, j);

            % Extract the array from 'finalConcat'
            trackLengthsMatrix{i, j} = finalConcat(currentIndex:(currentIndex+currentLength-1));

            % Update the index for the next iteration
            currentIndex = currentIndex + currentLength;
        end
    end
end
