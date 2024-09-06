function idx = findCellIndexWithChar(~,cellArray, charToFind)
    % Compare each element in the cell array with the character/string to find
    matches = strcmp(cellArray, charToFind);
    
    % Find the index of the first match
    idx = find(matches, 1);
    
end