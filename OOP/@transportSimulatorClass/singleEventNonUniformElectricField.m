function obj = singleEventNonUniformElectricField(obj, solverClass, geometryClass,xyzPosition, depositedEnergy, cloudInitialRadius, nPseudoCarrier)
            
                 
            obj.depositedEnergy = depositedEnergy; 
            obj.nPseudoCarrier = nPseudoCarrier;
            
            % Event initial position 
            x = xyzPosition.x;
            y = xyzPosition.y;
            z = xyzPosition.z;

            
            %% Compute charge cloud parameters for electrons 
            D_e = obj.materialProperties.mu_e*(obj.k*obj.materialProperties.T/obj.q_e);                 % Diffusion coefficient of the semiconductor in thermal equilibrium  
            N = depositedEnergy/obj.materialProperties.ehE;                                             % Number of electron-hole pairs generated after energy deposition 
            eps = obj.materialProperties.epsr*obj.eps0;                                                 % Absolute permittivity of the detector material, expressed in F/m 
            
            sigmaDiffusion_e = sqrt(2*D_e);                                                            % Diffusion only solution 
            R_0_e = ((3*obj.materialProperties.mu_e*obj.q_e*N)/(4*pi*eps)).^(1/3);                      % Repulsion only solution 
    
            nstep = 0; 
            
            ElectronCloud.x = [];
            ElectronCloud.y = [];
            ElectronCloud.z = [];
            
            xElectronCloud = x; 
            yElectronCloud = y; 
            zElectronCloud = z; 
            
            %% Baricenter trajectory for electrons 
            while (isInDetectorVolume(obj,x,y,z))
                
                % Compute position of all pseudo carrier 
                rElectronCloud = ((3*obj.materialProperties.mu_e*obj.q_e*N*nstep*obj.simulationParameters.timeStep)/(4*pi*eps)).^(1/3).*rand(nPseudoCarrier,1) + sqrt(2*D_e*nstep*obj.simulationParameters.timeStep).*randn(nPseudoCarrier,1);      % Uniform distribution 
                theta = rand(size(rElectronCloud))*360; 
                phi = rand(size(rElectronCloud))*360; 
                xElectronCloud = rElectronCloud.*cos(phi) + xElectronCloud; 
                yElectronCloud = rElectronCloud.*sin(phi) + yElectronCloud;
                zElectronCloud = rElectronCloud.*cos(theta) + zElectronCloud; 
                
                coord = horzcat([x;y;z],[xElectronCloud.'; yElectronCloud.'; zElectronCloud.']); 
                
                % Interpolate electric field data for all carrier position
                % and baricenter 
                % Sistemare mphinterp mettendo tutto in un unica linea 
                Ex = mphinterp(solverClass.model,{'es.Ex'}, 'coord', coord);
                Ey = mphinterp(solverClass.model,{'es.Ey'}, 'coord', coord);
                Ez = mphinterp(solverClass.model,{'es.Ez'}, 'coord', coord);
               
                x = -obj.materialProperties.mu_e*obj.simulationParameters.timeStep*Ex(1,1) + x; 
                y = -obj.materialProperties.mu_e*obj.simulationParameters.timeStep*Ey(1,1) + y; 
                z = -obj.materialProperties.mu_e*obj.simulationParameters.timeStep*Ez(1,1) + z; 
                
                % Compute cloud displacement due to drift 
                xElectronCloud = -obj.materialProperties.mu_e*obj.simulationParameters.timeStep.*(Ex(1,2:end).') + xElectronCloud; 
                yElectronCloud = -obj.materialProperties.mu_e*obj.simulationParameters.timeStep.*(Ey(1,2:end).') + yElectronCloud; 
                zElectronCloud = -obj.materialProperties.mu_e*obj.simulationParameters.timeStep.*(Ez(1,2:end).') + zElectronCloud;
                
                % update cloud position   
                ElectronCloud.x = vertcat(ElectronCloud.x, xElectronCloud);
                ElectronCloud.y = vertcat(ElectronCloud.y, yElectronCloud);
                ElectronCloud.z = vertcat(ElectronCloud.z, zElectronCloud);
                
                nstep = nstep + 1; 
            end 
            
            figure("Name","Trajectories")
            plotVolumes(geometryClass)
            hold on
            plot3(ElectronCloud.x,ElectronCloud.y,ElectronCloud.z)
            
    

           
        end 