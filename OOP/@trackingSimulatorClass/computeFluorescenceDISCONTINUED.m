function obj = computeFluorescence(obj, perElementPerEnergyAbsorbed)
%DESCRIPTION: 
%
%

% assign to each energy in each element an absorption shell 
absorbedInShell = calculateAbsorption(obj, perElementPerEnergyAbsorbed); 

% based on the absorption shell evaluate the fluorescence yield and
% determine the fluorescence photon energy from IXAS data

fluorescenceEnergy = calculateFluorescenceEnergy(obj, absorbedInShell, perElementPerEnergyAbsorbed); 

% Concatenate array and evaluate energy deposition
fluorescenceEnergy = ConcatenateCellArrayColumnRow(obj,fluorescenceEnergy); 


fluorescenceEnergy = adjustArrayLength(obj, fluorescenceEnergy, length(obj.energyDepositionMap{end}.energy));

if (obj.energyDepositionMap{end}.energy - fluorescenceEnergy) > 0
obj.energyDepositionMap{end}.energy = obj.energyDepositionMap{end}.energy - fluorescenceEnergy;
else 
    obj.energyDepositionMap{end}.energy = obj.energyDepositionMap{end}.energy;
end 

%Set to 0 out of volume events


last = length(obj.energyEmissionMap);

obj.energyEmissionMap{last+1}.energy = fluorescenceEnergy;
obj.energyEmissionMap{last+1}.x = obj.energyDepositionMap{end}.x;
obj.energyEmissionMap{last+1}.y = obj.energyDepositionMap{end}.y;
obj.energyEmissionMap{last+1}.z = obj.energyDepositionMap{end}.z;

obj.energyEmissionMap{last+1} = isOutOfVolume(obj, obj.energyEmissionMap{last+1}, obj.detectorVolume{1,1}.origin, obj.detectorVolume{1,1}.vertexes); 
obj.energyDepositionMap{end} = isOutOfVolume(obj, obj.energyDepositionMap{end}, obj.detectorVolume{1,1}.origin, obj.detectorVolume{1,1}.vertexes); 

obj.energyEmissionMap{last+1}.angles(1,:) = randi([0, 360], size(obj.energyEmissionMap{last+1}.x));
obj.energyEmissionMap{last+1}.angles(2,:) = randi([0, 360], size(obj.energyEmissionMap{last+1}.x));



end 