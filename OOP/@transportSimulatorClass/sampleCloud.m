function [TrajCloud] = sampleCloud(obj, trackingSimulatorClass,nPseudoCarrier, driftTime, diffusionTerm, repulsionTerm, position, carrierType, nonzeroidx)
    % This function calculates the trajectory of electron clouds based on
    % diffusion and repulsion terms.
    
    if strcmp(carrierType, 'electrons')
        D = obj.materialProperties.mu_e*(obj.k*obj.materialProperties.T/obj.q_e);    
    elseif strcmp(carrierType, 'holes')  
        D = obj.materialProperties.mu_h*(obj.k*obj.materialProperties.T/obj.q_e);    
    end 

    % Computed the numberof generated electorn-hole pair based on the
    % energy of each deposited event

    N = []; 
    if obj.trackingSimulationFlag == 1
        for i = 2 : length(trackingSimulatorClass.energyDepositionMap)
            N = horzcat(N, trackingSimulatorClass.energyDepositionMap{1,i}.energy(nonzeroidx{i-1})./obj.materialProperties.ehE); 
        end 
    else 
        N = obj.userEnergyDepositionMap.E./obj.materialProperties.ehE; 
    end 

    N = repmat(N, nPseudoCarrier,1); 

    initialCloudTerm = diffusionTerm;
    
    photonEnergy = [];
    if obj.trackingSimulationFlag == 1 
        for i = 2 : length(trackingSimulatorClass.energyDepositionMap)
            photonEnergy = horzcat(photonEnergy, trackingSimulatorClass.energyDepositionMap{1,i}.energy(nonzeroidx{i-1})); 
        end 
    else
        photonEnergy = obj.userEnergyDepositionMap.E; 
    end 

    
    photonEnergy = repmat(photonEnergy, nPseudoCarrier,1); 
    
    A = obj.materialProperties.initialCloudRadius(1);
    B = obj.materialProperties.initialCloudRadius(2);
    C = obj.materialProperties.initialCloudRadius(3);

    sigmaInitialCloud = (A.*(photonEnergy/1000).*(1-B./(1+C.*(photonEnergy/1000))))*1e-6; 
    
    rInitialCloud = abs(initialCloudTerm.*sigmaInitialCloud);   

    eps = obj.materialProperties.epsr*obj.eps0; % Absolute permittivity of the detector material, expressed in F/m
    sigmaDiffusion = sqrt(2*D.*driftTime); % Analytical expression for the standard deviation of the diffusion-only solution
    if strcmp(carrierType, 'electrons')
        R_0 = ((3*obj.materialProperties.mu_e*obj.q_e.*N.*driftTime)./(4*pi*eps)).^(1/3); 
    elseif strcmp(carrierType, 'holes') 
        R_0 = ((3*obj.materialProperties.mu_h*obj.q_e.*N.*driftTime)./(4*pi*eps)).^(1/3); 
    end 

    rDiffElectron = abs(diffusionTerm.*sigmaDiffusion);
    rRepElectron = repulsionTerm.*R_0;

    rTrajCloud = rDiffElectron + rRepElectron + rInitialCloud;

    theta = rand(size(rTrajCloud))*360;
    phi = rand(size(rTrajCloud))*360;

    TrajCloud.x = rTrajCloud.*sin(theta).*cos(phi) + position.x;
    TrajCloud.y = rTrajCloud.*sin(theta).*sin(phi) + position.y;
    TrajCloud.z = rTrajCloud.*cos(theta) + position.z;

end
