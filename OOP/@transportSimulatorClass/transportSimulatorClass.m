classdef transportSimulatorClass
    %TRANSPORTSIMULATOR 
    %   Simulator of charge transport in the detector volume 
    %   Input parameters are:  
    %   
    
    properties
        materialProperties 
        simulationParameters
        detectorVolume
        depositedEnergy
        nPseudoCarrier
        CCE
        CCEMap
        transients
        SingleEventTransportReport
        MultipleEventTransportReport

        electronCloudTrajectory 
        holeCloudTrajectory
        
        electronCloudVelocity
        holeCloudVelocity

        userEnergyDepositionMap
        trackingSimulationFlag

        % Physical constants 
        q_e = 1.6e-19;              % Electron elemental charge, expressed in C  
        eps0 = 8.85e-12;            % Vacuum permittivity, expressed in F/m 
        k = 1.38e-23;               % Boltzmann constant, expressed n m^2*kg*s^-2*K^-1 
    end
    
    % Constructor method 
    methods
        function obj = transportSimulatorClass(ehE, mu_e, tau_e, mu_h, tau_h, T,trackingSimulatorClass, epsr, initialCloudRadius, Fano, generationEnergy) 
            % Build simulation parameter table 
            obj.materialProperties.mu_e = mu_e;       % e-h pair mean generation energy, expressed in eV  
            obj.materialProperties.tau_e = tau_e;     % mobility of electrons, expressed in m^2/(Vs)
            obj.materialProperties.mu_h = mu_h;       % mobility of electrons, expressed in m^2/(Vs) 
            obj.materialProperties.tau_h = tau_h;     % mean lifetime of electrons, expressed in s
            obj.materialProperties.T = T;             % mean lifetime of electrons, expressed in s
            obj.materialProperties.ehE = ehE;         % absolute temperature expressed in K 
            obj.materialProperties.epsr = epsr;       % relative permittivity of the material 
            obj.materialProperties.initialCloudRadius = initialCloudRadius;  % parameters for initial cloud radius
            obj.materialProperties.Fano = Fano; 
            obj.materialProperties.generationEnergy = generationEnergy; 

            obj.detectorVolume = trackingSimulatorClass.detectorVolume;    %
        end   
        
        % Initialize simulation parameters 
        obj = initSimulationParameters(obj, timeStep)
        
      
        obj = singleEventNonUniformElectricField(obj, solverClass, geometryClass,xyzPosition, depositedEnergy, cloudInitialRadius, nPseudoCarrier)

        obj = constantElectricFieldWaveforms(obj, trackingSimulatorClass, geometryClass, FEMComsolSolverClass,cloudInitialRadius, nPseudoCarrier, eField)

        obj = constantElectricFieldStatic(obj, trackingSimulatorClass, FEMComsolSolverClass, nPseudoCarrier, eField)

        obj = scanCCEconstantElectricFieldStatic(obj, trackingSimulatorClass, FEMComsolSolverClass, eField, nPseudoCarrier, userEnergyDepositionMap)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        function [isInDetectorVolume] = isInDetectorVolume(obj, xCoordinate, yCoordinate, zCoordinate)
            
            isInDetectorVolume = zeros(length(obj.detectorVolume), 1); 
            for n = 1 : length(obj.detectorVolume)
                if ((xCoordinate <= obj.detectorVolume{n}.vertexes(1) + obj.detectorVolume{n}.origin(1)) && (xCoordinate >= obj.detectorVolume{n}.origin(1))) && ((yCoordinate <= obj.detectorVolume{n}.vertexes(2) + obj.detectorVolume{n}.origin(2)) && (yCoordinate >= obj.detectorVolume{n}.origin(2))) && ((zCoordinate <= obj.detectorVolume{n}.vertexes(3) + obj.detectorVolume{n}.origin(3)) && (zCoordinate >= obj.detectorVolume{n}.origin(3)))
                    isInDetectorVolume(n,1) = 1; 
                end 
            end 
            isInDetectorVolume = any(isInDetectorVolume); 
        end 

        function [isInDetectorVolume] = isInDetectorVolumeVectorized(obj, Coordinates)
            
            isInDetectorVolumeArray = cell(length(obj.detectorVolume), 1); 
            for n = 1 : length(obj.detectorVolume)
                isInDetectorVolumeArray{n} = ((Coordinates.x <= obj.detectorVolume{n}.vertexes(1) + obj.detectorVolume{n}.origin(1)) & (Coordinates.x >= obj.detectorVolume{n}.origin(1))) & ((Coordinates.y <= obj.detectorVolume{n}.vertexes(2) + obj.detectorVolume{n}.origin(2)) & (Coordinates.y >= obj.detectorVolume{n}.origin(2))) & ((Coordinates.z <= obj.detectorVolume{n}.vertexes(3) + obj.detectorVolume{n}.origin(3)) & (Coordinates.z >= obj.detectorVolume{n}.origin(3)));
            end       
            isInDetectorVolume = zeros(size(Coordinates)); 
            for n = 1 : length(obj.detectorVolume)
            isInDetectorVolume = isInDetectorVolumeArray{n} | isInDetectorVolume; 
            end 
        end 
        
        function [xClosest, yClosest, zClosest] = closestToDetectorVolume(obj, xCoordinate, yCoordinate, zCoordinate)
            xClosest = zeros(length(obj.detectorVolume), 1); 
            yClosest = zeros(length(obj.detectorVolume), 1); 
            zClosest = zeros(length(obj.detectorVolume), 1); 
            
            for n = 1 : length(obj.detectorVolume)
                % rectangle coordinates and vertexes
                x1 = obj.detectorVolume{n}.origin(1);
                y1 = obj.detectorVolume{n}.origin(2);
                z1 = obj.detectorVolume{n}.origin(3);
                x2 = obj.detectorVolume{n}.vertexes(1);
                y2 = obj.detectorVolume{n}.vertexes(2);
                z2 = obj.detectorVolume{n}.vertexes(3);

                % Prompt user for the position
                px = xCoordinate;
                py = yCoordinate;
                pz = zCoordinate;

                % Find the closest point
                xClosest(n) = max(min(px, x2), x1);
                yClosest(n) = max(min(py, y2), y1);
                zClosest(n) = max(min(pz, z2), z1);
            end 
            distance(n) = sqrt(xClosest(n).^2 + yClosest(n).^2 + zClosest(n).^2); 
           
            % Find the index of the minimum distance 
            minIndex = find(distance == min(distance));
            
            xClosest = xClosest(minIndex);
            yClosest = yClosest(minIndex);
            zClosest = zClosest(minIndex);
        end 

        function [closestCoordinates] = closestToDetectorVolumeVectorized(obj, Coordinates)

            xClosest = cell(length(obj.detectorVolume), 1); 
            yClosest = cell(length(obj.detectorVolume), 1); 
            zClosest = cell(length(obj.detectorVolume), 1); 
           
            
            closestCoordinates.x = zeros(size(Coordinates.x)); 
            closestCoordinates.y = zeros(size(Coordinates.y)); 
            closestCoordinates.z = zeros(size(Coordinates.z)); 

            for n = 1 : length(obj.detectorVolume)
                % rectangle coordinates and vertexes
                x1 = obj.detectorVolume{n}.origin(1);
                y1 = obj.detectorVolume{n}.origin(2);
                z1 = obj.detectorVolume{n}.origin(3);
                x2 = obj.detectorVolume{n}.vertexes(1);
                y2 = obj.detectorVolume{n}.vertexes(2);
                z2 = obj.detectorVolume{n}.vertexes(3);

                % Prompt user for the position
                px = Coordinates.x;
                py = Coordinates.y;
                pz = Coordinates.z;

                % Find the closest point
                xClosest{n} = max(min(px, x2), x1);
                yClosest{n} = max(min(py, y2), y1);
                zClosest{n} = max(min(pz, z2), z1);
            end 

                closestCoordinates.x = xClosest{n}; 
                closestCoordinates.y = yClosest{n}; 
                closestCoordinates.z = zClosest{n}; 
      

            

        end 
        
        function [weigthingPotential] = interpolateWeightingPotential(FEMComsolSolverClass,Coordinates)
 
           coords(1,:) = Coordinates.x(:)';
           coords(2,:) = Coordinates.y(:)';
           coords(3,:) = Coordinates.z(:)';
                        
           weigthingPotential = mphinterp(FEMComsolSolverClass.model,{'V'}, 'coord', coords);    
           weigthingPotential = reshape(weigthingPotential, size(Coordinates.x)); 
        end 
    end
end

