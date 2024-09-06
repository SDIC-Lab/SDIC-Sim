function obj = initializeEEM(obj,nEvents)
    % DESCRIPTION: 
    % Initialize Energy Emission Map (EEM), xyz position and energy of each event according to initial
    % distribution definition 
    
    obj.nEvents = nEvents; 

    if strcmp(obj.initialDistribution.type, "uniform rectangular")
            
        uniformSample = rand(3,nEvents); 
        obj.energyEmissionMap{1}.x = zeros(1,nEvents);
        obj.energyEmissionMap{1}.y = zeros(1,nEvents);
        obj.energyEmissionMap{1}.z = zeros(1,nEvents);
        
        obj.energyEmissionMap{1}.x = obj.initialDistribution.shape.origin(1) + uniformSample(1,:).*obj.initialDistribution.shape.vertexes(1);
        obj.energyEmissionMap{1}.y = obj.initialDistribution.shape.origin(2) + uniformSample(2,:).*obj.initialDistribution.shape.vertexes(2);
        obj.energyEmissionMap{1}.z = obj.initialDistribution.shape.origin(3) + uniformSample(3,:).*obj.initialDistribution.shape.vertexes(3);
        
        clear uniformSample; 
        
        % Populate the energy according to initial spectrum data 
        totalCounts = sum(obj.initialDistribution.spectrum.counts, 'all');
        normalizedCounts = obj.initialDistribution.spectrum.counts/totalCounts; 
        counts = round(normalizedCounts.*nEvents); 
       if sum(counts, 'all') ~= nEvents 
           diff = nEvents - sum(counts, 'all');
           rand_idx = randi(length(counts)); 
           counts(rand_idx) = abs(counts(rand_idx) + diff); 
       end 
       
       % 
       % Initialize an empty energy vector
       energyVector = [];
        
       % Create the vector by repeating elements based on their abundances
       for i = 1:length(obj.initialDistribution.spectrum.energy)
            energyVector = [energyVector, repmat(obj.initialDistribution.spectrum.energy(i), 1, counts(i))];
       end

       obj.energyEmissionMap{1}.energy = energyVector; 
       
       
       % Inizitialize the angles for each photon for initial distribution
       % definition 
       obj.energyEmissionMap{1}.angles(1,:) = repmat(obj.initialDistribution.angles(1), 1, obj.nEvents);
       obj.energyEmissionMap{1}.angles(2,:) = repmat(obj.initialDistribution.angles(2), 1 , obj.nEvents);

       obj.energyDepositionMap{1}.x = zeros(1,obj.nEvents); 
       obj.energyDepositionMap{1}.y = zeros(1,obj.nEvents); 
       obj.energyDepositionMap{1}.z = zeros(1,obj.nEvents); 
       obj.energyDepositionMap{1}.energy = zeros(1,obj.nEvents); 
    end 



    
end 

