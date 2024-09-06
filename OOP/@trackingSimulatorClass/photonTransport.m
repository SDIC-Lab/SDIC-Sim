function obj = photonTransport(obj)
    %DESCRIPTION 
    %

    % Evaluate absorption interpolating XCOM X-ray cross-section data to
    % generate an array of sampled track length, in an array with the same
    % dimension of the EEM
    [trackLengths, obj] = sampleTrackLength(obj); 

    % Evolve photon position based on sampled track lengths taking into
    % account angles
    
    last = length(obj.energyDepositionMap);

    phi = obj.energyEmissionMap{end}.angles(1,:);
    theta = obj.energyEmissionMap{end}.angles(2,:);

    % Update the next index of the energy deposition map 
    obj.energyDepositionMap{last+1}.x = obj.energyEmissionMap{end}.x + trackLengths.*sin(phi).*cos(theta);
    obj.energyDepositionMap{last+1}.y = obj.energyEmissionMap{end}.y + trackLengths.*sin(phi).*sin(theta);
    obj.energyDepositionMap{last+1}.z = obj.energyEmissionMap{end}.z + trackLengths.*cos(theta);
    obj.energyDepositionMap{last+1}.energy = obj.energyEmissionMap{end}.energy; 
    
    % Discard out of volume events 
    obj.energyDepositionMap{end} = isOutOfVolume(obj, obj.energyDepositionMap{end}, obj.detectorVolume{1,1}.origin, obj.detectorVolume{1,1}.vertexes); 


end

