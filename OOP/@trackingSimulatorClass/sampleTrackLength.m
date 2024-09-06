function [trackLengths, obj] = sampleTrackLength(obj)
    %DESCRIPTION: 
    
    % log-log Interpolate XCOM data mass attenuation coefficient for each
    % spectrum energy and each element of the detector material 

    
    for i = 1 : length(obj.materialProperties.elementList)
        
        for j = 1 : length(obj.materialProperties.elementList)
            xdata = obj.materialProperties.XCOMdata(j).energy;
            xdata = modifyRepeatingElements(obj, xdata, 5); 
            ydata = obj.materialProperties.XCOMdata(j).photoelectric;
            queryEnergy = obj.energyEmissionMap{end}.energy;
            for k = 1: length(queryEnergy)
                if queryEnergy(k)> obj.materialProperties.productionCut 
                    interpValuePar(j,k) = exp(interp1(log(xdata),log(ydata),log(queryEnergy(k).*(10^-6)),"linear"));
                else 
                    interpValuePar(j,k) = 0; 
                end 
            end 
        end 
         %weighted sum of all elements cross-sections
         weigth = transpose(str2double(obj.materialProperties.abundanceList)./sum(str2double(obj.materialProperties.abundanceList), 'all')); 
         interpValue(i,:) = sum(times(interpValuePar, weigth),1); 
    end 
    
    %depending of relavtive abundace of each element and the photon energy
    %of the spectrum, sample the track length from a poisson distribution
    %with the interpolated mass attenuation coefficient

    % (str2double(obj.materialProperties.abundanceList)); 

    %For each photon energy upload the per element probability based o mass
    %attenuation coefficient data at its energy 
    
    perElementProbability = zeros(length((str2double(obj.materialProperties.abundanceList))), length(queryEnergy));
    
    fractionByWeigth = (str2double(obj.materialProperties.abundanceList)); 

   
   absorbingElement = zeros(1,length(queryEnergy)); 

    for i = 1 : length(queryEnergy)
        if queryEnergy(i) ~= 0
            for j = 1 : length((str2double(obj.materialProperties.abundanceList)))
                perElementProbability(j,i) = fractionByWeigth(1,j).*((interpValuePar(j,i))./interpValue(1,i)); 
            end 
            % number of outcomes 
            noc = length(fractionByWeigth); 
            probabilities = perElementProbability(:,i);
            num_samples = 1;
            if ~isempty(probabilities) && min(probabilities) > 0
                absorbingElement(1,i) = randsample(noc, num_samples, true, probabilities); 
            end 
        end 
    end 

   massAttenuationCoeff = zeros(1,length(absorbingElement)); 
   
   for i = 1 : length(absorbingElement)
       if queryEnergy(i) ~= 0 && absorbingElement(i) ~= 0
          massAttenuationCoeff(i) = interpValuePar(round(absorbingElement(1,i)),i); 
       end 
   end 
 
  
   trackLengths = computeTracks(obj,massAttenuationCoeff); 
   obj.absorbingElement{1,end+1} = absorbingElement; 
    

end
