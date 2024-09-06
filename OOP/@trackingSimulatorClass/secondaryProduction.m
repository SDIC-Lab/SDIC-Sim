function obj = secondaryProduction(obj)
%DESCRIPTION:
%

% assign to each energy in each element an absorption shell 
absorbedInShell = calculateAbsorption(obj); 

% based on the absorption shell evaluate the fluorescence yield and
% determine the fluorescence photon energy from IXAS data

fluorescenceEnergy = calculateFluorescenceEnergy(obj, absorbedInShell); 

% Subtract emitted energy 
obj.energyDepositionMap{end}.energy = obj.energyDepositionMap{end}.energy - fluorescenceEnergy; 

last = length(obj.energyEmissionMap);

obj.energyEmissionMap{last+1}.energy = fluorescenceEnergy;
obj.energyEmissionMap{last+1}.x = obj.energyDepositionMap{end}.x;
obj.energyEmissionMap{last+1}.y = obj.energyDepositionMap{end}.y;
obj.energyEmissionMap{last+1}.z = obj.energyDepositionMap{end}.z;



obj.energyEmissionMap{last+1}.angles(1,:) = randi([0, 360], size(obj.energyEmissionMap{last+1}.x));
obj.energyEmissionMap{last+1}.angles(2,:) = randi([0, 360], size(obj.energyEmissionMap{last+1}.x));



end

