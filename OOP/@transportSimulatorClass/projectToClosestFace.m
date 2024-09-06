function projPoint = projectToClosestFace(obj,position)

    point = [position.x position.y position.x];
    % Calculate the vertices of the parallelepiped
    origin = obj.detectorVolume{1,1}.origin; 
    v1 = [obj.detectorVolume{1,1}.vertexes(1) 0 0];
    v2 = [0 obj.detectorVolume{1,1}.vertexes(2) 0];
    v3 = [0 0 obj.detectorVolume{1,1}.vertexes(3)];

    A = origin;
    B = origin + v1;
    C = origin + v2;
    D = origin + v3;
    E = B + v2;
    F = B + v3;
    G = C + v3;
    H = E + v3;

    % Define the faces using the vertices
    % Each face is defined by a point and a normal vector
    faces = {
        {A, cross(v1, v2)}, {A, cross(v2, v3)}, {A, cross(v3, v1)}, ... % Faces adjacent to the origin
        {D, cross(v2, v1)}, {B, cross(v3, v2)}, {C, cross(v1, v3)}  % Opposite faces
    };

    % Initialize variables to track the closest projection
    minDist = inf;
    projPoint = [0, 0, 0];

    % Iterate through each face to find the closest projection
    for i = 1:length(faces)
        face = faces{i};
        facePoint = face{1};
        faceNormal = face{2};

        % Project the point onto the face
        proj = projectPointOnPlane(point, facePoint, faceNormal);

        % Calculate distance from the original point to the projection
        dist = norm(proj - point);

        % Update if this is the closest projection so far
        if dist < minDist
            minDist = dist;
            projPoint = proj;
        end
    end
end

function projPoint = projectPointOnPlane(point, planePoint, planeNormal)
    % Normalize the plane normal vector
    n = planeNormal / norm(planeNormal);
    % Calculate the projection of the point onto the plane
    projPoint = point + dot(planePoint - point, n) * n;
end
