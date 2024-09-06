function newPosition = projectDetectorVolume(obj, initialPosition, planeNormal, planePoint, lifeTime, driftTime)
    
    initialPosition_x = initialPosition.x;
    initialPosition_y = initialPosition.y;
    initialPosition_z = initialPosition.z;

    %Define plane normal and plane point

    
    % Use arrayfun to process each coordinate pair
    [proj_x, proj_y, proj_z] = arrayfun(@(lifeTime, driftTime, x, y, z) ...
        projectPointToPlane(obj, planeNormal, planePoint, lifeTime, driftTime, x, y, z), ...
        lifeTime, driftTime, initialPosition_x, initialPosition_y, initialPosition_z);
    
    % Update electronPosition with new coordinates
    newPosition.x = proj_x;
    newPosition.y = proj_y;
    newPosition.z = proj_z;

    newPosition.x0 = initialPosition.x0;
    newPosition.y0 = initialPosition.y0;
    newPosition.z0 = initialPosition.z0;
end


