 function obj = constantElectricFieldStatic(obj, trackingSimulatorClass, FEMComsolSolverClass, nPseudoCarrier, electricField)
    
    obj.trackingSimulationFlag = 1; 
    % number of events to be transported  
    nEvents = 0; 
    for i = 2 : length(trackingSimulatorClass.energyDepositionMap)
        nonzeroidx{i-1} = find(trackingSimulatorClass.energyDepositionMap{1,i}.energy > 0);
        nEvents = nEvents + length(trackingSimulatorClass.energyDepositionMap{1,i}.energy(nonzeroidx{i-1}));
    end 

    % Sample lifetimes for each carrier, based to mean carrier lifetime

    electronLifeTime = exprnd(obj.materialProperties.tau_e,nPseudoCarrier,nEvents);

    holeLifeTime = exprnd(obj.materialProperties.tau_h,nPseudoCarrier,nEvents);


    % Load energy deposition map from a trackingSimulatorClass object 

    initialPosition.x = [];
    initialPosition.y = [];
    initialPosition.z = []; 

    for i = 2 : length(trackingSimulatorClass.energyDepositionMap)
        initialPosition.x = horzcat(initialPosition.x,trackingSimulatorClass.energyDepositionMap{1,i}.x(nonzeroidx{i-1}));
        initialPosition.y = horzcat(initialPosition.y,trackingSimulatorClass.energyDepositionMap{1,i}.y(nonzeroidx{i-1}));
        initialPosition.z = horzcat(initialPosition.z,trackingSimulatorClass.energyDepositionMap{1,i}.z(nonzeroidx{i-1}));
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
    
    electronCloud = sampleCloud(obj, trackingSimulatorClass,nPseudoCarrier, electronDriftTime, diffusionTerm, repulsionTerm, electronPosition, 'electrons',nonzeroidx); 
    
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

    holeCloud = sampleCloud(obj, trackingSimulatorClass,nPseudoCarrier, holeDirftTime, diffusionTerm, repulsionTerm, holePosition, 'holes',nonzeroidx); 
    
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
    coords = real(coords);
    weigthingPotential = mphinterp(FEMComsolSolverClass.model,{'V'}, 'coord', coords);    
    weigthingPotentialElectron = reshape(weigthingPotential, size(electronCloud.x)); 
    
    coords(1,:) = holeCloud.x(:)';
    coords(2,:) = holeCloud.y(:)';
    coords(3,:) = holeCloud.z(:)';
    coords = real(coords);
    weigthingPotential = mphinterp(FEMComsolSolverClass.model,{'V'}, 'coord', coords);    
    weigthingPotentialHole = reshape(weigthingPotential, size(holeCloud.x)); 


    coords(1,:) = initialPosition.x(:)';
    coords(2,:) = initialPosition.y(:)';
    coords(3,:) = initialPosition.z(:)';
    coords = real(coords);
    weigthingPotential = mphinterp(FEMComsolSolverClass.model,{'V'}, 'coord', coords);    
    weigthingPotentialInitialPosition = reshape(weigthingPotential, size(initialPosition.x)); 

    % Evaluate elementary induced charge contribution
    
    electronCharge = (1./nPseudoCarrier).*(weigthingPotentialElectron-weigthingPotentialInitialPosition); 
    holeCharge = (1./nPseudoCarrier).*(-weigthingPotentialHole+weigthingPotentialInitialPosition); 
    % Generate CCE map 
    electronholecharge = electronCharge + holeCharge; 
    obj.CCE = sum(electronholecharge,1);  
            
    
    clear obj.CCEmap
    
    obj.CCEMap = zeros(length(trackingSimulatorClass.energyDepositionMap{1,2}.x),5*(length(trackingSimulatorClass.energyDepositionMap)-1)); 
    cceidx = 1; 

    for i = 2 : length(trackingSimulatorClass.energyDepositionMap)
        obj.CCEMap(nonzeroidx{i-1},5*(i-2)+1) = trackingSimulatorClass.energyDepositionMap{1,i}.x(nonzeroidx{i-1});
        obj.CCEMap(nonzeroidx{i-1},5*(i-2)+2) = trackingSimulatorClass.energyDepositionMap{1,i}.y(nonzeroidx{i-1});
        obj.CCEMap(nonzeroidx{i-1},5*(i-2)+3) = trackingSimulatorClass.energyDepositionMap{1,i}.z(nonzeroidx{i-1});
        obj.CCEMap(nonzeroidx{i-1},5*(i-2)+4) = trackingSimulatorClass.energyDepositionMap{1,i}.energy(nonzeroidx{i-1});
        if i > 2
            cceidx = length(nonzeroidx{i-2})+cceidx; 
        end 
        obj.CCEMap(nonzeroidx{i-1},5*(i-2)+5) = obj.CCE(1,cceidx:((cceidx-1)+length(nonzeroidx{i-1}))); 
    end 
                
                        
                        
 end