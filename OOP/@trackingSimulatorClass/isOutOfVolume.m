function filteredStruct = isOutOfVolume(~, inputStruct, origin, vertex)
    % inputStruct: Array of structs with fields 'x', 'y', and 'z'.
    % origin: The origin of the rectangle [x, y, z].
    % vertex: The 3D dimensions of the rectangle [vx, vy, vz] 

   

    % Loop through each element in the struct array
    for i = 1:length(inputStruct.x)
        % Access the fields directly
        x = inputStruct.x(i);
        y = inputStruct.y(i);
        z = inputStruct.z(i);

        % Check if the coordinate is within the defined rectangle
        if x >= origin(1) && x <= origin(1) + vertex(1) && ...
           y >= origin(2) && y <= origin(2) + vertex(2) && ...
           z >= origin(3) && z <= origin(3) + vertex(3)
           
            % Add the coordinate to the filtered array
            filteredStruct.x(i) = inputStruct.x(i);
            filteredStruct.y(i) = inputStruct.y(i);
            filteredStruct.z(i) = inputStruct.z(i);
            filteredStruct.energy(i) = inputStruct.energy(i);  
        else
            filteredStruct.x(i) = 0;
            filteredStruct.y(i) = 0;
            filteredStruct.z(i) = 0;
            filteredStruct.energy(i) = 0;  
        end 
    end
end