function obj = recordEnergySpectra(obj,transportSimulatorClass,energyThreshold)
    
    
   for j = 1 : length(transportSimulatorClass.CCEMap(:,1))
       for i = 1 : length(transportSimulatorClass.CCEMap(1,:))/5
           partialEnergy(i) = transportSimulatorClass.CCEMap(j,4+(i-1)*5)*transportSimulatorClass.CCEMap(j,5+(i-1)*5);
       end 
        recordedEnergy = sum(partialEnergy,'all'); 
        %Add Fano noise and electronic noise 
        eps = transportSimulatorClass.materialProperties.generationEnergy;
        F = transportSimulatorClass.materialProperties.Fano; 

        if recordedEnergy <0 
            recordedEnergy = 0; 
        end 

        sigmaElecronic = eps*obj.noiseParameters.ENC; 
        sigmaFano = sqrt(F*eps*recordedEnergy); 
        
        recordedEnergy = normrnd(0,sigmaFano,1) + normrnd(0,sigmaElecronic,1) + recordedEnergy; 
        
        nBin = length(obj.energySpectrum.energy) ;
        eMax = obj.energySpectrum.energy(end); 

        indx = round((recordedEnergy*nBin)/eMax); 
        
        if indx > round((energyThreshold*nBin)/eMax) && indx <= nBin   
            obj.energySpectrum.counts(indx) = obj.energySpectrum.counts(indx)  + 1; 
        end 
    end 
    
    

end