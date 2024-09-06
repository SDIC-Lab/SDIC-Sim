function [projX, projY, projZ] = projectPointToPlane(obj,planeNormal, planePoint, lifeTime, driftTime, x, y, z)
    % x, y, z: The coordinates of the point to project
    % planeNormal: The normal vector of the plane, e.g., [0, 0, 1] for the xy-plane
    % planePoint: A point through which the plane passes, e.g., [0, 0, z0] for the xy-plane at z = z0

    % Create a vector for the input point
    point = [x y z];

    origin = obj.detectorVolume{1,1}.origin; 
    vertex = obj.detectorVolume{1,1}.vertexes;
    
    if lifeTime > driftTime ||  isOutOfVolumeBoolean(obj, x, y, z, origin, vertex)
        % Vector from planePoint to point
        pointVector = point - planePoint;
        
        % Distance from point to plane along the plane normal
        distance = dot(pointVector, planeNormal) / norm(planeNormal);
        
        % Projection of the point onto the plane
        projection = point - distance * planeNormal;
        
        % Extract projected coordinates
        projX = projection(1);
        projY = projection(2);
        projZ = projection(3);
   
    else

        projX = x;
        projY = y;
        projZ = z;

    end 

end
