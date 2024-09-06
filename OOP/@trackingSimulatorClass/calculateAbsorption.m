

function absorbedInShell = calculateAbsorption(obj)
   
   nonzeroidx = find(obj.energyDepositionMap{1,end}.energy ~= 0); 
   absorbedEnergy = obj.energyDepositionMap{1,end}.energy(nonzeroidx); 
   absorbingElement = obj.absorbingElement{1,end}(nonzeroidx); 
   absorbingElement = absorbingElement(find(absorbingElement ~= 0)); 
   absorbingElementIndex = unique(absorbingElement);
   absorbedInShell = zeros(size(obj.energyDepositionMap{1,end}.energy)); 
   temp3absorbedInShell = zeros(size(nonzeroidx)); 
   if ~isempty(absorbingElement)
       for i = 1 : length(absorbingElementIndex)
           %Extract from IXAS data K edges and Jump Factor for absorbing material 
           idx = findCellIndexWithChar(obj,obj.materialProperties.IXASdata{1,i}(3:end,2), ' Energy (eV)');
           Kedges = str2double(obj.materialProperties.IXASdata{1,i}(3:idx,2));
           Jump = str2double(obj.materialProperties.IXASdata{1,i}(3:idx,5));
           matidx = find(absorbingElement == absorbingElementIndex(i)); 
           energyDepositionInMaterial = absorbedEnergy(matidx); 
           uniqueEnergies = unique(energyDepositionInMaterial);
    
           temp2AbsorbedinShell = zeros(size(energyDepositionInMaterial));
           % Find a sublist of available K-edges based on energy
           for j = 1 : length(uniqueEnergies)
              enidx = find(energyDepositionInMaterial == uniqueEnergies(j));  
              sameEnergyDepositedinMaterial = energyDepositionInMaterial(enidx); 
              kidx = find(Kedges < uniqueEnergies(j));
              % Perform a random sample to assign a sheel to each event 
              noc = length(kidx);     % numer of outcomes 
              % Construct the probability space based on Jump Factors
              availableShellsJump = Jump(kidx); 
              probabilities = zeros(size(availableShellsJump));
              probabilities(1,1) = (availableShellsJump(1,1)-1)/availableShellsJump(1,1);
              for k = 2 : length(probabilities)
                  probabilities(k,1) = prod(1./availableShellsJump(1:k-1,1)).*((availableShellsJump(k,1)-1)/availableShellsJump(k,1)); 
              end 
              probabilities(end,1) = 1-sum(probabilities);
              num_samples = length(sameEnergyDepositedinMaterial);
              restrictedShells = (randsample(noc, num_samples, true, probabilities)).';
              temp1AbsorbedinShell = (restrictedShells-1) + kidx(1,1); 
              temp2AbsorbedinShell(enidx) = temp1AbsorbedinShell;
           end 
           temp3absorbedInShell(matidx) = temp2AbsorbedinShell; 
       end 
   end 
   absorbedInShell(nonzeroidx) = temp3absorbedInShell;
end 

  
 