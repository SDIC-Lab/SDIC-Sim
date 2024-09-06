function furthestPoint = findFurthestPointInsideCube(obj, initial, final)

   origin = obj.detectorVolume{1,1}.origin; 
    v1 = [obj.detectorVolume{1,1}.vertexes(1) 0 0];
    v2 = [0 obj.detectorVolume{1,1}.vertexes(2) 0];
    v3 = [0 0 obj.detectorVolume{1,1}.vertexes(3)];
    
   initial = [initial.x initial.y initial.x];
   final = [final.x final.y final.x];


    % Direction vector of the trajectory
    direction = final - initial;
    
    % Parametric line equation components
    p = initial;
    d = direction;
    
    % Initialize the furthest point as the initial point if it lies within the parallelepiped
    if isPointInsideParallelepiped(initial, origin, v1, v2, v3)
        furthestPoint = initial;
    else
        furthestPoint = []; % No part of the trajectory is inside the parallelepiped
        return;
    end
    
    % Check for intersections with each face of the parallelepiped and update the furthestPoint accordingly
    faces = defineParallelepipedFaces(origin, v1, v2, v3);
    maxDistance = 0;
    for i = 1:length(faces)
        [intersect, point] = checkIntersection(p, d, faces{i});
        if intersect
            distance = norm(point - initial);
            if distance > maxDistance && isPointInsideParallelepiped(point, origin, v1, v2, v3)
                maxDistance = distance;
                furthestPoint = point;
            end
        end
    end
end

function inside = isPointInsideParallelepiped(point, origin, v1, v2, v3)
    % This function checks if a given point is inside the parallelepiped
    % Placeholder for actual implementation
    inside = true; % Simplification, implement proper logic
end

function faces = defineParallelepipedFaces(origin, v1, v2, v3)
    % This function defines the faces of the parallelepiped for intersection checks
    % Placeholder for actual implementation
    faces = {}; % Simplification, implement proper logic to define faces
end

function [intersect, point] = checkIntersection(p, d, face)
    % This function checks if the line intersects with a given face and returns the intersection point
    % Placeholder for actual implementation
    intersect = false; % Simplification, implement proper logic
    point = [0, 0, 0]; % Placeholder
end
