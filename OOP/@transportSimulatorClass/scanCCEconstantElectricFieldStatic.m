function obj = scanCCEconstantElectricFieldStatic(obj, trackingSimulatorClass, FEMComsolSolverClass, electricField, nPseudoCarrier, userEnergyDepositionMap)
    
    % Import detector volume 

    % number of events to be transported  
    nEvents = length(userEnergyDepositionMap.x); 

    obj.userEnergyDepositionMap = userEnergyDepositionMap; 
    obj.trackingSimulationFlag = 0; 

    % Sample lifetimes for each carrier, based to mean carrier lifetime

    electronLifeTime = exprnd(obj.materialProperties.tau_e,nPseudoCarrier,nEvents);

    holeLifeTime = exprnd(obj.materialProperties.tau_h,nPseudoCarrier,nEvents);


    % Load energy deposition map from a trackingSimulatorClass object 

    initialPosition.x = [];
    initialPosition.y = [];
    initialPosition.z = []; 

    for i = 2 : length(userEnergyDepositionMap.x)
        initialPosition.x = userEnergyDepositionMap.x;
        initialPosition.y = userEnergyDepositionMap.y;
        initialPosition.z = userEnergyDepositionMap.z;
    end

    initialPosition.x = repmat(initialPosition.x, nPseudoCarrier,1); 
    initialPosition.y = repmat(initialPosition.y, nPseudoCarrier,1); 
    initialPosition.z = repmat(initialPosition.z, nPseudoCarrier,1); 

    electronPosition.x0 = initialPosition.x; 
    electronPosition.y0 = initialPosition.y; 
    electronPosition.z0 = initialPosition.z;

    holePosition.x0 = initialPosition.x; 
    holePosition.y0 = initialPosition.y; 
    holePosition.z0 = initialPosition.z;


    electronCloudDisplacement = obj.materialProperties.mu_e*nonzeros(electricField).*electronLifeTime; 
    holeCloudDisplacement = -obj.materialProperties.mu_h*nonzeros(electricField).*holeLifeTime; 

    displacementFlag = double(electricField ~= 0);

    electronPosition.x = electronPosition.x0 + electronCloudDisplacement.*displacementFlag(1); 
    electronPosition.y = electronPosition.y0 + electronCloudDisplacement.*displacementFlag(2); 
    electronPosition.z = electronPosition.z0 + electronCloudDisplacement.*displacementFlag(3); 

    electronPhantomTrajectory.x = electronPosition.x;
    electronPhantomTrajectory.y = electronPosition.y;
    electronPhantomTrajectory.z = electronPosition.z;
    

    holePosition.x = holePosition.x0 + holeCloudDisplacement.*displacementFlag(1); 
    holePosition.y = holePosition.y0 + holeCloudDisplacement.*displacementFlag(2); 
    holePosition.z = holePosition.z0 + holeCloudDisplacement.*displacementFlag(3);

    holePhantomTrajectory.x = holePosition.x;
    holePhantomTrajectory.y = holePosition.y;
    holePhantomTrajectory.z = holePosition.z;
    
    % Confine final transport position in detector volume  

    newPosition = confineInDetectorVolume(electronPosition, obj);

    electronPosition.x = newPosition.x;
    electronPosition.y = newPosition.y;
    electronPosition.z = newPosition.z;

    newPosition = confineInDetectorVolume(holePosition, obj);

    holePosition.x = newPosition.x;
    holePosition.y = newPosition.y;
    holePosition.z = newPosition.z;

    % Evaluate drift time for all cloud baricenters
    electronDriftTime = (displacementFlag(1).*abs(electronPosition.x - electronPosition.x0) + displacementFlag(2) ... 
        .*abs(electronPosition.y - electronPosition.y0) + displacementFlag(3).*abs(electronPosition.z - electronPosition.z0)) ... 
        ./(nonzeros(electricField)*obj.materialProperties.mu_e); 
    
     holeDirftTime = (displacementFlag(1).*abs(holePosition.x - holePosition.x0) + displacementFlag(2) ... 
        .*abs(holePosition.y - holePosition.y0) + displacementFlag(3).*abs(holePosition.z - holePosition.z0)) ... 
        ./(nonzeros(electricField)*obj.materialProperties.mu_h); 

        
    diffusionTerm = randn(nPseudoCarrier,nEvents);     % Gaussian distribution 
    repulsionTerm = rand(nPseudoCarrier,nEvents);      % Uniform distribution 
    
    electronCloud = sampleCloud(obj, trackingSimulatorClass,nPseudoCarrier, electronDriftTime, diffusionTerm, repulsionTerm, electronPosition, 'electrons'); 
    
    electronCloud.x0 = initialPosition.x;
    electronCloud.y0 = initialPosition.y;
    electronCloud.z0 = initialPosition.z;

    if nonzeros(electricField) > 0
        planeNormal = double(electricField ~= 0); 
        planePoint = obj.detectorVolume{1,1}.vertexes(find(electricField ~= 0));
    elseif  nonzeros(electricField) < 0
        planeNormal = double(electricField ~= 0); 
        planePoint = obj.detectorVolume{1,1}.origin(find(electricField ~= 0));
    end 
   
    electronCloud = projectDetectorVolume(obj, electronCloud, planeNormal, planePoint, electronLifeTime, electronDriftTime);

    electronCloud = confineInDetectorVolume(electronCloud, obj);

    holeCloud = sampleCloud(obj, trackingSimulatorClass,nPseudoCarrier, holeDirftTime, diffusionTerm, repulsionTerm, holePosition, 'holes'); 
    
    holeCloud.x0 = initialPosition.x;
    holeCloud.y0 = initialPosition.y;
    holeCloud.z0 = initialPosition.z;

    if nonzeros(electricField) > 0
        planeNormal = double(electricField ~= 0); 
        planePoint = obj.detectorVolume{1,1}.origin(find(electricField ~= 0));
    elseif  nonzeros(electricField) < 0
        planeNormal = double(electricField ~= 0); 
        planePoint = obj.detectorVolume{1,1}.vertexes(find(electricField ~= 0));
    end 
    

    holeCloud = projectDetectorVolume(obj, holeCloud, planeNormal, planePoint, holeLifeTime, holeDirftTime);
    
    holeCloud = confineInDetectorVolume(holeCloud, obj);
   

    coords(1,:) = electronCloud.x(:)';
    coords(2,:) = electronCloud.y(:)';
    coords(3,:) = electronCloud.z(:)';
    
    weigthingPotential = mphinterp(FEMComsolSolverClass.model,{'V'}, 'coord', coords);    
    weigthingPotentialElectron = reshape(weigthingPotential, size(electronCloud.x)); 
    
    coords(1,:) = holeCloud.x(:)';
    coords(2,:) = holeCloud.y(:)';
    coords(3,:) = holeCloud.z(:)';
    
    weigthingPotential = mphinterp(FEMComsolSolverClass.model,{'V'}, 'coord', coords);    
    weigthingPotentialHole = reshape(weigthingPotential, size(holeCloud.x)); 


    coords(1,:) = initialPosition.x(:)';
    coords(2,:) = initialPosition.y(:)';
    coords(3,:) = initialPosition.z(:)';
    
    weigthingPotential = mphinterp(FEMComsolSolverClass.model,{'V'}, 'coord', coords);    
    weigthingPotentialInitialPosition = reshape(weigthingPotential, size(initialPosition.x)); 

    % Evaluate elementary induced charge contribution
    
    electronCharge = (1./nPseudoCarrier).*(weigthingPotentialElectron-weigthingPotentialInitialPosition); 
    holeCharge = (1./nPseudoCarrier).*(-weigthingPotentialHole+weigthingPotentialInitialPosition); 
    % Generate CCE map 
    electronholecharge = electronCharge + holeCharge; 
    obj.CCE = sum(electronholecharge,1);  
            
    
   
                
                        
                        
 end