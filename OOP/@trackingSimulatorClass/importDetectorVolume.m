function obj = importDetectorVolume(obj, GeometryClassObject)
    %DESCRIPTION: 
    %Import detector volume geometry from a geometryClass object 
    
    % Look for keyword detector volume of type "merge", look in the
    % mergelist and import the type "cube" geometry object 
    [~, size1] = size(GeometryClassObject.selections); 
    
    % import the detector volume 
    for sz = 1 : size1
        if strcmp(GeometryClassObject.selections(sz).name, 'detector_volume')
            for ssz = 1 : length(GeometryClassObject.selections(sz).mergelist)
                for sssz = 1 : size1
                    if strcmp(GeometryClassObject.volumes(sssz).name, GeometryClassObject.selections(sz).mergelist{ssz})
                        obj.detectorVolume{1} = GeometryClassObject.volumes(sssz); 
                    end 
                end 
            end 
        end 
    end 

    
end

