function obj = constantElectricFieldWaveforms(obj, trackingSimulatorClass, geometryClass, FEMComsolSolverClass,cloudInitialRadius, nPseudoCarrier, eField)
            %Generate livetimes for each pseudocarrier 
            nEvents = length(trackingSimulatorClass.absorbedPrimary(:,1)) + length(trackingSimulatorClass.absorbedSecondary(:,1)) + length(trackingSimulatorClass.absorbedTertiary(:,1));
            % for electrons 
            electronLiveTime = exprnd(obj.materialProperties.tau_e,nPseudoCarrier,nEvents);

            % for holes 
            holeLiveTime = exprnd(obj.materialProperties.tau_h,nPseudoCarrier,nEvents);

            %Compute barycenter displacement for each livetime 

            initialPosition.x = repmat(vertcat(trackingSimulatorClass.absorbedPrimary(:,1),trackingSimulatorClass.absorbedSecondary(:,1),trackingSimulatorClass.absorbedTertiary(:,1)).',nPseudoCarrier,1);
            initialPosition.y = repmat(vertcat(trackingSimulatorClass.absorbedPrimary(:,2),trackingSimulatorClass.absorbedSecondary(:,2),trackingSimulatorClass.absorbedTertiary(:,2)).',nPseudoCarrier,1);
            initialPosition.z = repmat(vertcat(trackingSimulatorClass.absorbedPrimary(:,3),trackingSimulatorClass.absorbedSecondary(:,3),trackingSimulatorClass.absorbedTertiary(:,3)).',nPseudoCarrier,1);

            electronPosition.x0 = repmat(vertcat(trackingSimulatorClass.absorbedPrimary(:,1),trackingSimulatorClass.absorbedSecondary(:,1),trackingSimulatorClass.absorbedTertiary(:,1)).',nPseudoCarrier,1);
            electronPosition.y0 = repmat(vertcat(trackingSimulatorClass.absorbedPrimary(:,2),trackingSimulatorClass.absorbedSecondary(:,2),trackingSimulatorClass.absorbedTertiary(:,2)).',nPseudoCarrier,1);
            electronPosition.z0 = repmat(vertcat(trackingSimulatorClass.absorbedPrimary(:,3),trackingSimulatorClass.absorbedSecondary(:,3),trackingSimulatorClass.absorbedTertiary(:,3)).',nPseudoCarrier,1);

            holePosition.x0 = repmat(vertcat(trackingSimulatorClass.absorbedPrimary(:,1),trackingSimulatorClass.absorbedSecondary(:,1),trackingSimulatorClass.absorbedTertiary(:,1)).',nPseudoCarrier,1);
            holePosition.y0 = repmat(vertcat(trackingSimulatorClass.absorbedPrimary(:,2),trackingSimulatorClass.absorbedSecondary(:,2),trackingSimulatorClass.absorbedTertiary(:,2)).',nPseudoCarrier,1);
            holePosition.z0 = repmat(vertcat(trackingSimulatorClass.absorbedPrimary(:,3),trackingSimulatorClass.absorbedSecondary(:,3),trackingSimulatorClass.absorbedTertiary(:,3)).',nPseudoCarrier,1);

            electronBarycenterDisplacement = obj.materialProperties.mu_e*eField.*electronLiveTime; 
            holeBarycenterDisplacement = -obj.materialProperties.mu_h*eField.*holeLiveTime; 

            electronPosition.z = electronBarycenterDisplacement + electronPosition.z0; 
            holePosition.z = holeBarycenterDisplacement + holePosition.z0; 

            electronPosition.x = electronPosition.x0; 
            electronPosition.y = electronPosition.y0; 

            holePosition.x = holePosition.x0; 
            holePosition.y = holePosition.y0; 

            %% for electrons 


            ClosestCoordinates = closestToDetectorVolumeVectorized(obj, electronPosition); 


            electronDirftTime = abs(ClosestCoordinates.z - electronPosition.z0)/(eField*obj.materialProperties.mu_e); 
            
            nSteps = round(electronDirftTime/obj.simulationParameters.timeStep);
            
            electronTimeArray = (1:nSteps)*(electronDirftTime/nSteps); 
            electronPositionArray.z = electronTimeArray.*obj.materialProperties.mu_e*eField; 
            electronPositionArray.x = repmat(initialPosition.x,nSteps); 
            electronPositionArray.y = repmat(initialPosition.y,nSteps); 
            
            %% for holes

            ClosestCoordinates = closestToDetectorVolumeVectorized(obj, holePosition); 

            holeDirftTime = abs(ClosestCoordinates.z - holePosition.z0)/(eField*obj.materialProperties.mu_e); 
            
            holeTimeArray = (1:nSteps)*(holeDirftTime/nSteps); 
            holePositionArray.z = -holeTimeArray.*obj.materialProperties.mu_e*eField; 
            holePositionArray.x = repmat(initialPosition.x,nSteps); 
            holePositionArray.y = repmat(initialPosition.y,nSteps); 

            %% Compute displacement vector for charge cloud broadening 

            diffusionTerm = randn(nPseudoCarrier,nEvents);     % Gaussian distribution 
            repulsionTerm = rand(nPseudoCarrier,nEvents);      % Uniform distribution 

            %% for electrons 

            D_e = obj.materialProperties.mu_e*(obj.k*obj.materialProperties.T/obj.q_e);    % Diffusion coefficient of the semiconductor in thermal equilibrium  
            % Number of carrier for each event 
            N = (vertcat(trackingSimulatorClass.absorbedPrimary(:,4),trackingSimulatorClass.absorbedSecondary(:,4),trackingSimulatorClass.absorbedTertiary(:,4))/obj.materialProperties.ehE).';          
            N = repmat(N,nPseudoCarrier,1);                                                % Number of electron-hole pairs generated after energy deposition 
            eps = obj.materialProperties.epsr*obj.eps0;                                    % Absolute permittivity of the detector material, expressed in F/m 

            electronSigmaDiffusion = sqrt(2*D_e.*electronTimeArray);                                          % Analytical expression for the standard deviation of the diffusion-only solution
            electronR_0 = ((3*obj.materialProperties.mu_e*obj.q_e.*N*electronTimeArray)/(4*pi*eps)).^(1/3);    % Analytical expression of the radius of the sphere of the repulsion-only solution 

            rDiffElectron = abs(diffusionTerm.*electronSigmaDiffusion); 
            rRepElectron = repulsionTerm.*electronR_0; 

            rTrajElectronCloud = rDiffElectron + rRepElectron; 

            theta = rand(size(rTrajElectronCloud))*360; 
            phi = rand(size(rTrajElectronCloud))*360; 

            TrajElectronCloud.x = rTrajElectronCloud.*cos(phi) + electronPosition.x; 
            TrajElectronCloud.y = rTrajElectronCloud.*sin(phi) + electronPosition.y;
            TrajElectronCloud.z = rTrajElectronCloud.*cos(theta) + electronPosition.z; 


            ClosestCoordinates = closestToDetectorVolumeVectorized(obj, TrajElectronCloud); 



            TrajElectronCloud.x = ClosestCoordinates.x; 
            TrajElectronCloud.y = ClosestCoordinates.y; 
            TrajElectronCloud.z = ClosestCoordinates.z;


            %% for holes

            D_h = obj.materialProperties.mu_h*(obj.k*obj.materialProperties.T/obj.q_e);    % Diffusion coefficient of the semiconductor in thermal equilibrium  
            % Number of carrier for each event 
            N = (vertcat(trackingSimulatorClass.absorbedPrimary(:,4),trackingSimulatorClass.absorbedSecondary(:,4),trackingSimulatorClass.absorbedTertiary(:,4))/obj.materialProperties.ehE).';          
            N = repmat(N,nPseudoCarrier,1);                                                % Number of electron-hole pairs generated after energy deposition 
            eps = obj.materialProperties.epsr*obj.eps0;                                    % Absolute permittivity of the detector material, expressed in F/m 

            holeSigmaDiffusion = sqrt(2*D_h.*holeDirftTime);                                          % Analytical expression for the standard deviation of the diffusion-only solution
            holeR_0 = ((3*obj.materialProperties.mu_h*obj.q_e.*N.*holeDirftTime)/(4*pi*eps)).^(1/3);    % Analytical expression of the radius of the sphere of the repulsion-only solution 

            rDiffHole = abs(diffusionTerm.*holeSigmaDiffusion); 
            rRepHole = repulsionTerm.*holeR_0; 

            rTrajHoleCloud = rDiffHole + rRepHole; 

            theta = rand(size(rTrajHoleCloud))*360; 
            phi = rand(size(rTrajHoleCloud))*360; 

            TrajHoleCloud.x = rTrajHoleCloud.*cos(phi) + holePosition.x; 
            TrajHoleCloud.y = rTrajHoleCloud.*sin(phi) + holePosition.y;
            TrajHoleCloud.z = rTrajHoleCloud.*cos(theta) + holePosition.z; 



            ClosestCoordinates = closestToDetectorVolumeVectorized(obj, TrajHoleCloud); 

            TrajHoleCloud.x = ClosestCoordinates.x; 
            TrajHoleCloud.y = ClosestCoordinates.y; 
            TrajHoleCloud.z = ClosestCoordinates.z;   

        end 