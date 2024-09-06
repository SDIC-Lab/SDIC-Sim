function [furthestPoint_x, furthestPoint_y, furthestPoint_z] = findFurthestPointInsideVolume(obj, initial_x, initial_y, initial_z, final_x, final_y, final_z)
    
   origin = obj.detectorVolume{1,1}.origin; 
   v1 = [obj.detectorVolume{1,1}.vertexes(1) 0 0];
   v2 = [0 obj.detectorVolume{1,1}.vertexes(2) 0];
   v3 = [0 0 obj.detectorVolume{1,1}.vertexes(3)];

   initial = [initial_x initial_y initial_z];
   final = [final_x final_y final_z];

    % Bounds in each direction
    xMin = origin(1);
    xMax = origin(1) + v1(1);
    yMin = origin(2);
    yMax = origin(2) + v2(2);
    zMin = origin(3);
    zMax = origin(3) + v3(3);
    
    % Direction of the line
    dir = final - initial;
    
    % Parametric form of the line: initial + t*dir, solve for t where the line intersects each plane
    tVals = [];
    if dir(1) ~= 0
        tXMin = (xMin - initial(1)) / dir(1);
        tXMax = (xMax - initial(1)) / dir(1);
        tVals = [tVals, tXMin, tXMax];
    end
    if dir(2) ~= 0
        tYMin = (yMin - initial(2)) / dir(2);
        tYMax = (yMax - initial(2)) / dir(2);
        tVals = [tVals, tYMin, tYMax];
    end
    if dir(3) ~= 0
        tZMin = (zMin - initial(3)) / dir(3);
        tZMax = (zMax - initial(3)) / dir(3);
        tVals = [tVals, tZMin, tZMax];
    end
    
    % Filter t values to find those that are within the segment [0, 1]
    tVals = tVals(tVals >= 0 & tVals <= 1);
    
    % If no t values are within [0, 1], the final point is the furthest within the volume
    if isempty(tVals)
        furthestPoint = final;
    else
        % Find the furthest valid t value and calculate the corresponding point
        tMax = max(tVals);
        furthestPoint = initial + tMax * dir;
    end
    furthestPoint_x = furthestPoint(1);
    furthestPoint_y = furthestPoint(2);
    furthestPoint_z = furthestPoint(3);

end