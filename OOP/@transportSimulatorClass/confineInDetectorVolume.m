function electronPosition = confineInDetectorVolume(electronPosition, obj)
    % Extract initial positions
    initial_x = electronPosition.x0; 
    initial_y = electronPosition.y0; 
    initial_z = electronPosition.z0; 

    % Extract final positions
    final_x = electronPosition.x;  
    final_y = electronPosition.y; 
    final_z = electronPosition.z; 
    
    % Use arrayfun to process each coordinate pair
    [out_x, out_y, out_z] = arrayfun(@(ix, iy, iz, fx, fy, fz) ...
        findFurthestPointInsideVolume(obj, ix, iy, iz, fx, fy, fz), ...
        initial_x, initial_y, initial_z, final_x, final_y, final_z);
    
    % Update electronPosition with new coordinates
    electronPosition.x = out_x;
    electronPosition.y = out_y;
    electronPosition.z = out_z;
end