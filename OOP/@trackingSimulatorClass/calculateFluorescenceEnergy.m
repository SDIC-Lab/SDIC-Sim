function fluorescenceEnergies = calculateFluorescenceEnergy(obj, absorbedInShell)
    %DESCRIPTION: 
 
    nonzeroidx = find(obj.energyDepositionMap{1,end}.energy ~= 0); 
    absorbingElement = obj.absorbingElement{1,end}(nonzeroidx); 
    absorbingElement = absorbingElement(find(absorbingElement ~= 0)); 
    energies = obj.energyDepositionMap{1,end}.energy(nonzeroidx); 
    absorbingElementIndex = unique(absorbingElement);
    fluorescenceEnergies = zeros(size(absorbedInShell)); 
    temp3fluorescenceEnergies = zeros(size(nonzeroidx));
    temp1absorbedInShell = absorbedInShell(nonzeroidx);
    for i = 1 : length(absorbingElementIndex)
        matidx = find(absorbingElement == absorbingElementIndex(i));
        energiesperMaterial = energies(matidx); 
        absorbedInShellperElement = temp1absorbedInShell(matidx);
        uniqueShells = unique(absorbedInShellperElement); 
        temp2fluorescenceEnergies = zeros(size(absorbedInShellperElement)); 
        for j = 1 : length(uniqueShells)
            shellString = "("+obj.materialProperties.IXASdata{1,i}(2+uniqueShells(j),1);
            shellidx = find(absorbedInShellperElement == uniqueShells(j));
            energiesperShell = energiesperMaterial(shellidx); 
            startidx = find(contains(obj.materialProperties.IXASdata{1,i}(:,1),'Line'));
            strindx = find(contains(obj.materialProperties.IXASdata{1,i}(startidx+1:end,1),shellString)); 
            strindx = strindx + startidx;
            intensity = str2double(obj.materialProperties.IXASdata{1,i}(strindx,3));
            uniqueEnergiesperShell = unique(energiesperShell); 
            fluorescenceYield = str2double(obj.materialProperties.IXASdata{1,i}(2+uniqueShells(j),4)); 
            fluorescenceEnergyValues = str2double(obj.materialProperties.IXASdata{1,i}(strindx,2)); 
            temp0fluorescenceEnergies = zeros(size(shellidx)); 
            emitting = 2.*ones(size(shellidx));
            % Sample the indexes emitting a fluorescence photon 
            if fluorescenceYield > 0
                noc = 2;
                probabilities = [fluorescenceYield 1-fluorescenceYield];
                num_samples = length(shellidx);
                emitting = (randsample(noc, num_samples, true, probabilities)).';
            end 
            temp1fluorescenceEnergies = zeros(size(shellidx));
            % Sample the fluorescence energy 
            if ~isempty(intensity)
                noc = length(intensity);
                probabilities = intensity;
                num_samples = length(shellidx);
                fluorescenceEnergiesIndex = (randsample(noc, num_samples, true, probabilities)).'; 
                % Set to zero photons not emitting 
                for kk = 1 : length(temp1fluorescenceEnergies)
                    if emitting(1,kk) == 1
                        flag = 1;
                    elseif emitting(1,kk) == 2
                        flag = 0; 
                    end 
                    temp1fluorescenceEnergies(1,kk) = flag.*fluorescenceEnergyValues(round(fluorescenceEnergiesIndex(1,kk)));
                end
            end          
        temp2fluorescenceEnergies(shellidx) = temp1fluorescenceEnergies;    
        end 
        temp3fluorescenceEnergies(matidx) = temp2fluorescenceEnergies; 
    end 
    fluorescenceEnergies(nonzeroidx) = temp3fluorescenceEnergies; 
end
    