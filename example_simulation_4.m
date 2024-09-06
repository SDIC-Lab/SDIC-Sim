%-------- Computing stopping power of a 5mm W block  
close all; 
clear all; 
clc;

% add path for simulator repository
addpath("\\10.79.34.83\Laboratorio\4 - Quercia\HighZsim\OOPv2")

myGeometry = geometryClass(); 

pitch = 475e-6; 
arrayDims = [1, 3];
interpixelGap = 25e-6;
guardDepth = 5e-4;
thickness = 10e-3; 
worldOffset = 1e-3;

myGeometry = buildPixelDetector(myGeometry, pitch, arrayDims, interpixelGap, guardDepth, thickness, worldOffset); 


plotVolumes(myGeometry)
hold on 
plotcube(myGeometry,[1e-3 1e-3 0],[0.25e-3 0.75e-3 0])
xlim([-1e-3 3.5e-3])
ylim([-1e-3 3.5e-3])
hold off


% Tracking simulation object 
myTracker = trackingSimulatorClass(); 

% Specify the fraction by weight of each element
myTracker = importMaterialProperties(myTracker, "1 W", 19.28, 1e3); 

myTracker = importDetectorVolume(myTracker, myGeometry); 

spectrum.energy = 100e3;
spectrum.counts = 1;

myTracker = defineInitialDistribution(myTracker,"uniform rectangular",[0.25e-3 0.75e-3 0],[1e-3 1e-3 0], spectrum,[0 0]); 

nPhoton = 10000; 

myTracker.energyThreshold = 0.1e-3; 

myTracker = initializeEEM(myTracker,nPhoton); 


% Perform n tracking cycles exectuing this function n times 
myTracker = particleTracking(myTracker,nPhoton);


% %%%%%%%%%%%%%%%%%%%%%%%%% transport simulation parameter definition %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% ehE = 4.46;       % e-h pair mean generation energy, expressed in eV           
% mu_e = mu_e_array;      % mobility of electrons, expressed in m^2/(Vs)
% mu_h  = 0.0115;     % mobility of electrons, expressed in m^2/(Vs)
% 
% tau_e  = tau_e_array;    % mean lifetime of electrons, expressed in s  
% tau_h  = 1.5e-6;    % mean lifetime of electrons, expressed in s
% 
% T = 298;          % Absolute temperature expressed in K 
% 
% epsr = 10.2; 
% 
% Fano = 0.11;
% generationEnergy = 4.46; 
% 
% A = 0.95;             % 10^-6 m/keV
% B = 0.98;             % dimensionless
% C = 0.003;            % keV^-1
% 
% 
% CloudRadiusParameters = [A B C];            % Initial cloud radius parameter (gaussian distribution hypothesis)
% timeStep = 2e-9;
% eField = [0 0 Vbias/(1.85e-3)]; 
% % Construt a transporter obejct 
% myTransporter = transportSimulatorClass(ehE, mu_e, tau_e, mu_h, tau_h, T,myTracker,epsr,CloudRadiusParameters, Fano, generationEnergy);
% myTransporter = initSimulationParameters(myTransporter,timeStep);
% nPseudoCarrier = 50; 
% myTransporter = constantElectricFieldStatic(myTransporter, myTracker, mySolver1, nPseudoCarrier, eField); 
% 
% 
% myInstrument = recordEnergySpectra(myInstrument, myTransporter,energyThreshold); 
% 
% semilogy(myInstrument.energySpectrum.energy,myInstrument.energySpectrum.counts)  
% drawnow;  % Refresh plot 
% 
% disp("Simulation "+num2str(ctr)+" completed in "+num2str(toc)+" seconds.")
% accTime = accTime + toc; 
% ctr = ctr +1; 
% end 
% 
% s = seconds(accTime);
% s.Format = 'hh:mm:ss'; 
% disp("Simulation ended, total elapsed time is "+char(s)+" min")
% 
% % Save simulated spectra in specified folder 
% % Specify the directory where you want to create the new folder
% % For example, creating the folder in the current working directory:
% baseDir = "\\10.79.34.83\Laboratorio\4 - Quercia\Work\11 - Publications\quercia2024montecarlo\simulator validation\Simulation results"; % Change this to any directory you prefer
% 
% % Get the current date and time in a specific format
% % For example, 'yyyy-mm-dd_HH-MM-SS' format
% dateTimeNow = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
% 
% % Create the folder name by appending the date and time
% folderName = ['SimResults_' dateTimeNow];
% 
% % Full path for the new folder
% fullFolderPath = fullfile(baseDir, folderName);
% 
% % Create the folder
% mkdir(fullFolderPath);
% 
% % Save the current directory to return later
% previousDir = pwd;
% 
% % Change the current directory to the new folder
% cd(fullFolderPath);
% 
% SimulationReport.Tracker = myTracker;
% SimulationReport.Transporter = myTransporter;
% SimulationReport.Geometry = myGeometry;
% SimulationReport.FEMSolver = mySolver1; 
% SimulationReport.Instrument = myInstrument; 
% save("SimulationReport.mat", "SimulationReport")
% 
% energyspectrum(:,1) = myInstrument.energySpectrum.energy./1000; 
% energyspectrum(:,2) = myInstrument.energySpectrum.counts;
% save("energyspectrum.mat", "energyspectrum")
% % Return to the previous directory
% cd(previousDir);
% 
% end 
% 



