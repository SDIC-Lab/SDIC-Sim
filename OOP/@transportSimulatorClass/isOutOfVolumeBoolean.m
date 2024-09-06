function outofVolume = isOutOfVolumeBoolean(~, x, y, z, origin, vertex)


    % Check if the coordinate is within the defined 3D rectangle
    if x >= origin(1) && x <= origin(1) + vertex(1) && ...
       y >= origin(2) && y <= origin(2) + vertex(2) && ...
       z >= origin(3) && z <= origin(3) + vertex(3)
       
       outofVolume = false; 
    else
       outofVolume = true; 
    end 
end