close all; 
clear all; 
clc;


% Simulator validation with a broader range of energy spectra and detector
% geometry 

% Load experimental spectra 

experimentalData.Am241_05pitch = readtable("Experimental measurement for validation/Redlen_0.5 mm pitch_241Am.csv");
experimentalData.Fe55_05pitch = readtable("Experimental measurement for validation/Redlen_0.5 mm pitch_55Fe.csv");
experimentalData.Co57_05pitch = readtable("Experimental measurement for validation/Redlen_0.5 mm pitch_57Co.csv");
experimentalData.Am241_08pitch = readtable("Experimental measurement for validation/Acrorad_CdTe_0.8 mm pitch_241Am.csv");


experimentalData.Am241_05pitch(1, :) = [];
experimentalData.Am241_05pitch(8193:end,:) = []; 
experimentalData.Am241_08pitch(:,2) = [];

newColumnNames = {'Energy (keV)', 'Counts'};
experimentalData.Am241_05pitch.Properties.VariableNames = newColumnNames;
experimentalData.Fe55_05pitch.Properties.VariableNames = newColumnNames;
experimentalData.Co57_05pitch.Properties.VariableNames = newColumnNames;
experimentalData.Am241_08pitch.Properties.VariableNames = newColumnNames;

idx = find(experimentalData.Fe55_05pitch.("Energy (keV)") < 0.7);
experimentalData.Fe55_05pitch.("Counts")(idx) = zeros(size(idx)); 
idx = find(experimentalData.Fe55_05pitch.("Energy (keV)") > 7);
experimentalData.Fe55_05pitch.("Counts")(idx) = zeros(size(idx)); 

idx = find(experimentalData.Co57_05pitch.("Energy (keV)") < 4);
experimentalData.Co57_05pitch.("Counts")(idx) = zeros(size(idx)); 
idx = find(experimentalData.Co57_05pitch.("Energy (keV)") > 140);
experimentalData.Co57_05pitch.("Counts")(idx) = zeros(size(idx)); 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 1: 0.4 mm pitch detector 241 Am %%%%%%%%%%%%%
 
% add path for simulator repository
addpath("\\10.79.34.83\Laboratorio\4 - Quercia\HighZsim\OOPv2")

%Declare geometry 

myGeometry = geometryClass(); 

pitch = 475e-6; 
arrayDims = [1, 3];
interpixelGap = 25e-6;
guardDepth = 5e-4;
thickness = 1.85e-3; 
worldOffset = 1e-3;

myGeometry = buildPixelDetector(myGeometry, pitch, arrayDims, interpixelGap, guardDepth, thickness, worldOffset); 

% Solve for weighting field of pixel under test 

mySolver1 = FEMComsolSolverClass();

mySolver1 = createModel(mySolver1,"model1","det","geom1");


mySolver1 = importParams(mySolver1,myGeometry);


mySolver1 = importGeometry(mySolver1);

mySolver1 = ComsolFEMInitialization(mySolver1,'CZT',10.4,'Air',1);

mySolver1 = ComsolComputeWeigthingField(mySolver1,"electrode_pixel12");


figure("Name","Check mesh")
mphmesh(mySolver1.model,"mesh1", "facealpha", 0.5)


figure("Name","Weighting potential of electrode 1")
pg1 = mySolver1.model.result.create('pg1', 'PlotGroup3D');
pg1.feature.create('slice1', 'Slice');
mphplot(mySolver1.model,'pg1','rangenum',1)


Vbias =  700;

% Input photon flux energy distribution 
energy = 59.54e3;
counts = 1;

% Create an instrumentSimulatorClass to accumulate simulation results 
nBin = 4000;
eMax = 130e3; 
energyThreshold = 1e3; 
ENC = 30;

myInstrument = instrumentSimulatorClass(nBin, eMax, ENC); 

ctr = 1; 
accTime = 0; 

SimulationLoops = 25; 

figure("Name","energy spectra")

for i = 1 : SimulationLoops 
tic 
% Tracking simulation object 
myTracker = trackingSimulatorClass(); 

% Specify the fraction by weight of each element
myTracker = importMaterialProperties(myTracker, "0.45 Cd 0.05 Zn 0.5 Te", 5.65, 1e3); 

myTracker = importDetectorVolume(myTracker, myGeometry); 

spectrum.energy = energy;
spectrum.counts = counts;

myTracker = defineInitialDistribution(myTracker,"uniform rectangular",[0.25e-3 0.75e-3 0],[1e-3 1e-3 0], spectrum,[0 0]); 

nPhoton = 5000; 

myTracker.energyThreshold = 1e-3; 

myTracker = initializeEEM(myTracker,nPhoton); 


% Perform n tracking cycles exectuing this function n times 
myTracker = particleTracking(myTracker,nPhoton);
myTracker = particleTracking(myTracker,nPhoton);
myTracker = particleTracking(myTracker,nPhoton);


%%%%%%%%%%%%%%%%%%%%%%%%% transport simulation parameter definition %%%%%%%%%%%%%%%%%%%%%%%%%

ehE = 4.46;         % e-h pair mean generation energy, expressed in eV           
mu_e = 0.1;         % mobility of electrons, expressed in m^2/(Vs)
mu_h  = 0.0115;     % mobility of electrons, expressed in m^2/(Vs)

tau_e  =  2e-5;     % mean lifetime of electrons, expressed in s  
tau_h  = 1.5e-6;    % mean lifetime of electrons, expressed in s

T = 298;            % Absolute temperature expressed in K 

epsr = 10.2; 

Fano = 0.11;
generationEnergy = 4.46; 

A = 0.95;             % 10^-6 m/keV
B = 0.98;             % dimensionless
C = 0.003;            % keV^-1


CloudRadiusParameters = [A B C];            % Initial cloud radius parameter (gaussian distribution hypothesis)
timeStep = 2e-9;
eField = [0 0 Vbias/(1.85e-3)]; 
% Construt a transporter obejct 
myTransporter = transportSimulatorClass(ehE, mu_e, tau_e, mu_h, tau_h, T,myTracker,epsr,CloudRadiusParameters, Fano, generationEnergy);
myTransporter = initSimulationParameters(myTransporter,timeStep);
nPseudoCarrier = 100; 
myTransporter = constantElectricFieldStatic(myTransporter, myTracker, mySolver1, nPseudoCarrier, eField); 


myInstrument = recordEnergySpectra(myInstrument, myTransporter,energyThreshold); 

semilogy(myInstrument.energySpectrum.energy,myInstrument.energySpectrum.counts)  
drawnow;  % Refresh plot 

disp("Simulation "+num2str(ctr)+" completed in "+num2str(toc)+" seconds.")
accTime = accTime + toc; 
ctr = ctr +1; 
end 

s = seconds(accTime);
s.Format = 'hh:mm:ss'; 
disp("Simulation ended, total elapsed time is "+char(s)+" min")

% Save simulated spectra in specified folder 
% Specify the directory where you want to create the new folder
% For example, creating the folder in the current working directory:
baseDir = "\\10.79.34.83\Laboratorio\4 - Quercia\Work\11 - Publications\quercia2024montecarlo\simulator validation\Simulation results"; % Change this to any directory you prefer

% Get the current date and time in a specific format
% For example, 'yyyy-mm-dd_HH-MM-SS' format
dateTimeNow = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

% Create the folder name by appending the date and time
folderName = ['SimResults_' dateTimeNow];

% Full path for the new folder
fullFolderPath = fullfile(baseDir, folderName);

% Create the folder
mkdir(fullFolderPath);

% Save the current directory to return later
previousDir = pwd;

% Change the current directory to the new folder
cd(fullFolderPath);

SimulationReport.Tracker = myTracker;
SimulationReport.Transporter = myTransporter;
SimulationReport.Geometry = myGeometry;
SimulationReport.FEMSolver = mySolver1; 
SimulationReport.Instrument = myInstrument; 
save("SimulationReport.mat", "SimulationReport")

energyspectrum(:,1) = myInstrument.energySpectrum.energy./1000; 
energyspectrum(:,2) = myInstrument.energySpectrum.counts;
save("energyspectrum.mat", "energyspectrum")

%Normalize and compare with experimental results 
maximum = max(energyspectrum(:,2));
energyspectrum(:,2) = energyspectrum(:,2) ./maximum; 

maximum = max(experimentalData.Am241_05pitch(:,2));
experimentalData.Am241_05pitch(:,2) = experimentalData.Am241_05pitch(:,2) ./maximum; 

figure;
plot(energyspectrum(:,1), energyspectrum(:,2))
hold on 
semilogy(experimentalData.Am241_05pitch.("Energy (keV)"), experimentalData.Am241_05pitch.("Counts"));
xlim([0 65]);

title('Experimental vs. Simulation Redlen 0.5 mm pitch - Am 241')
savefig('Experimental vs. Simulation Redlen 0.5 mm pitch - Am 241.fig');
% Return to the previous directory
cd(previousDir);

clear energyspectrum; 
%%

%%%%%%%%%%%%% 2: 0.4 mm pitch detector Fe 55 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
% add path for simulator repository
addpath("\\10.79.34.83\Laboratorio\4 - Quercia\HighZsim\OOPv2")

%Declare geometry 

myGeometry = geometryClass(); 

pitch = 475e-6; 
arrayDims = [1, 3];
interpixelGap = 25e-6;
guardDepth = 5e-4;
thickness = 1.85e-3; 
worldOffset = 1e-3;

myGeometry = buildPixelDetector(myGeometry, pitch, arrayDims, interpixelGap, guardDepth, thickness, worldOffset); 

% Solve for weighting field of pixel under test 

mySolver1 = FEMComsolSolverClass();

mySolver1 = createModel(mySolver1,"model1","det","geom1");


mySolver1 = importParams(mySolver1,myGeometry);


mySolver1 = importGeometry(mySolver1);

mySolver1 = ComsolFEMInitialization(mySolver1,'CZT',10.4,'Air',1);

mySolver1 = ComsolComputeWeigthingField(mySolver1,"electrode_pixel12");


figure("Name","Check mesh")
mphmesh(mySolver1.model,"mesh1", "facealpha", 0.5)


figure("Name","Weighting potential of electrode 1")
pg1 = mySolver1.model.result.create('pg1', 'PlotGroup3D');
pg1.feature.create('slice1', 'Slice');
mphplot(mySolver1.model,'pg1','rangenum',1)


Vbias =  700;

% Input photon flux energy distribution 
energy = [5900.3 5889.1 5769.9 6491.8 6491.8 6537.0];
counts = [0.58134 0.29394 0.00027 0.08153 0.04223 0.00069];

% Create an instrumentSimulatorClass to accumulate simulation results 
nBin = 400;
eMax = 7e3; 
energyThreshold = 0.5e3; 
ENC = 30;

myInstrument = instrumentSimulatorClass(nBin, eMax, ENC); 

ctr = 1; 
accTime = 0; 

SimulationLoops = 25; 

figure("Name","energy spectra")

for i = 1 : SimulationLoops 
tic 
% Tracking simulation object 
myTracker = trackingSimulatorClass(); 

% Specify the fraction by weight of each element
myTracker = importMaterialProperties(myTracker, "0.45 Cd 0.05 Zn 0.5 Te", 5.65, 1e3); 

myTracker = importDetectorVolume(myTracker, myGeometry); 

spectrum.energy = energy;
spectrum.counts = counts;

myTracker = defineInitialDistribution(myTracker,"uniform rectangular",[0.25e-3 0.75e-3 0],[1e-3 1e-3 0], spectrum,[0 0]); 

nPhoton = 5000; 

myTracker.energyThreshold = 1e-3; 

myTracker = initializeEEM(myTracker,nPhoton); 


% Perform n tracking cycles exectuing this function n times 
myTracker = particleTracking(myTracker,nPhoton);
myTracker = particleTracking(myTracker,nPhoton);
myTracker = particleTracking(myTracker,nPhoton);


%%%%%%%%%%%%%%%%%%%%%%%%% transport simulation parameter definition %%%%%%%%%%%%%%%%%%%%%%%%%

ehE = 4.46;         % e-h pair mean generation energy, expressed in eV           
mu_e = 0.1;         % mobility of electrons, expressed in m^2/(Vs)
mu_h  = 0.0115;     % mobility of electrons, expressed in m^2/(Vs)

tau_e  =  2e-5;     % mean lifetime of electrons, expressed in s  
tau_h  = 1.5e-6;    % mean lifetime of electrons, expressed in s

T = 298;            % Absolute temperature expressed in K 

epsr = 10.2; 

Fano = 0.11;
generationEnergy = 4.46; 

A = 0.95;             % 10^-6 m/keV
B = 0.98;             % dimensionless
C = 0.003;            % keV^-1


CloudRadiusParameters = [A B C];            % Initial cloud radius parameter (gaussian distribution hypothesis)
timeStep = 2e-9;
eField = [0 0 Vbias/(1.85e-3)]; 
% Construt a transporter obejct 
myTransporter = transportSimulatorClass(ehE, mu_e, tau_e, mu_h, tau_h, T,myTracker,epsr,CloudRadiusParameters, Fano, generationEnergy);
myTransporter = initSimulationParameters(myTransporter,timeStep);
nPseudoCarrier = 100; 
myTransporter = constantElectricFieldStatic(myTransporter, myTracker, mySolver1, nPseudoCarrier, eField); 


myInstrument = recordEnergySpectra(myInstrument, myTransporter,energyThreshold); 

semilogy(myInstrument.energySpectrum.energy,myInstrument.energySpectrum.counts)  
drawnow;  % Refresh plot 

disp("Simulation "+num2str(ctr)+" completed in "+num2str(toc)+" seconds.")
accTime = accTime + toc; 
ctr = ctr +1; 
end 

s = seconds(accTime);
s.Format = 'hh:mm:ss'; 
disp("Simulation ended, total elapsed time is "+char(s)+" min")

% Save simulated spectra in specified folder 
% Specify the directory where you want to create the new folder
% For example, creating the folder in the current working directory:
baseDir = "\\10.79.34.83\Laboratorio\4 - Quercia\Work\11 - Publications\quercia2024montecarlo\simulator validation\Simulation results"; % Change this to any directory you prefer

% Get the current date and time in a specific format
% For example, 'yyyy-mm-dd_HH-MM-SS' format
dateTimeNow = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

% Create the folder name by appending the date and time
folderName = ['SimResults_' dateTimeNow];

% Full path for the new folder
fullFolderPath = fullfile(baseDir, folderName);

% Create the folder
mkdir(fullFolderPath);

% Save the current directory to return later
previousDir = pwd;

% Change the current directory to the new folder
cd(fullFolderPath);

SimulationReport.Tracker = myTracker;
SimulationReport.Transporter = myTransporter;
SimulationReport.Geometry = myGeometry;
SimulationReport.FEMSolver = mySolver1; 
SimulationReport.Instrument = myInstrument; 
save("SimulationReport.mat", "SimulationReport")

energyspectrum(:,1) = myInstrument.energySpectrum.energy./1000; 
energyspectrum(:,2) = myInstrument.energySpectrum.counts;
save("energyspectrum.mat", "energyspectrum")

%Normalize and compare with experimental results 
maximum = max(energyspectrum(:,2));
energyspectrum(:,2) = energyspectrum(:,2) ./maximum; 

maximum = max(experimentalData.Fe55_05pitch.("Counts"));
experimentalData.Fe55_05pitch.("Counts") = experimentalData.Fe55_05pitch.("Counts") ./maximum; 

figure;
plot(energyspectrum(:,1), energyspectrum(:,2))
hold on 
semilogy(experimentalData.Fe55_05pitch.("Energy (keV)"), experimentalData.Fe55_05pitch.("Counts"));
xlim([0 7])

title('Experimental vs. Simulation Redlen 0.5 mm pitch - Fe 55')
savefig('Experimental vs. Simulation Redlen 0.5 mm pitch - Fe 55.fig');
% Return to the previous directory
cd(previousDir);
clear energyspectrum; 

%%

%%%%%%%%%%%%%% 3: 0.4 mm pitch detector Co 57 %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add path for simulator repository
addpath("\\10.79.34.83\Laboratorio\4 - Quercia\HighZsim\OOPv2")

%Declare geometry 

myGeometry = geometryClass(); 

pitch = 475e-6; 
arrayDims = [1, 3];
interpixelGap = 25e-6;
guardDepth = 5e-4;
thickness = 1.85e-3; 
worldOffset = 1e-3;

myGeometry = buildPixelDetector(myGeometry, pitch, arrayDims, interpixelGap, guardDepth, thickness, worldOffset); 

% Solve for weighting field of pixel under test 

mySolver1 = FEMComsolSolverClass();

mySolver1 = createModel(mySolver1,"model1","det","geom1");


mySolver1 = importParams(mySolver1,myGeometry);


mySolver1 = importGeometry(mySolver1);

mySolver1 = ComsolFEMInitialization(mySolver1,'CZT',10.4,'Air',1);

mySolver1 = ComsolComputeWeigthingField(mySolver1,"electrode_pixel12");


figure("Name","Check mesh")
mphmesh(mySolver1.model,"mesh1", "facealpha", 0.5)


figure("Name","Weighting potential of electrode 1")
pg1 = mySolver1.model.result.create('pg1', 'PlotGroup3D');
pg1.feature.create('slice1', 'Slice');
mphplot(mySolver1.model,'pg1','rangenum',1)


Vbias =  700;

% Input photon flux energy distribution 
energy = 122.06e3;
counts = 1;

% Create an instrumentSimulatorClass to accumulate simulation results 
nBin = 8000;
eMax = 130e3; 
energyThreshold = 1e3; 
ENC = 30;

myInstrument = instrumentSimulatorClass(nBin, eMax, ENC); 

ctr = 1; 
accTime = 0; 

SimulationLoops = 50; 

figure("Name","energy spectra")

for i = 1 : SimulationLoops 
tic 
% Tracking simulation object 
myTracker = trackingSimulatorClass(); 

% Specify the fraction by weight of each element
myTracker = importMaterialProperties(myTracker, "0.45 Cd 0.05 Zn 0.5 Te", 5.65, 1e3); 

myTracker = importDetectorVolume(myTracker, myGeometry); 

spectrum.energy = energy;
spectrum.counts = counts;

myTracker = defineInitialDistribution(myTracker,"uniform rectangular",[0.25e-3 0.75e-3 0],[1e-3 1e-3 0], spectrum,[0 0]); 

nPhoton = 10000; 

myTracker.energyThreshold = 1e-3; 

myTracker = initializeEEM(myTracker,nPhoton); 


% Perform n tracking cycles exectuing this function n times 
myTracker = particleTracking(myTracker,nPhoton);
myTracker = particleTracking(myTracker,nPhoton);
myTracker = particleTracking(myTracker,nPhoton);


%%%%%%%%%%%%%%%%%%%%%%%%% transport simulation parameter definition %%%%%%%%%%%%%%%%%%%%%%%%%

ehE = 4.46;         % e-h pair mean generation energy, expressed in eV           
mu_e = 0.1;         % mobility of electrons, expressed in m^2/(Vs)
mu_h  = 0.01;     % mobility of electrons, expressed in m^2/(Vs)

tau_e  =  2e-5;     % mean lifetime of electrons, expressed in s  
tau_h  = 0.65e-6;    % mean lifetime of electrons, expressed in s

T = 298;            % Absolute temperature expressed in K 

epsr = 10.2; 

Fano = 0.11;
generationEnergy = 4.46; 

A = 0.95;             % 10^-6 m/keV
B = 0.98;             % dimensionless
C = 0.003;            % keV^-1


CloudRadiusParameters = [A B C];            % Initial cloud radius parameter (gaussian distribution hypothesis)
timeStep = 2e-9;
eField = [0 0 Vbias/(1.85e-3)]; 
% Construt a transporter obejct 
myTransporter = transportSimulatorClass(ehE, mu_e, tau_e, mu_h, tau_h, T,myTracker,epsr,CloudRadiusParameters, Fano, generationEnergy);
myTransporter = initSimulationParameters(myTransporter,timeStep);
nPseudoCarrier = 100; 
myTransporter = constantElectricFieldStatic(myTransporter, myTracker, mySolver1, nPseudoCarrier, eField); 


myInstrument = recordEnergySpectra(myInstrument, myTransporter,energyThreshold); 

semilogy(myInstrument.energySpectrum.energy,myInstrument.energySpectrum.counts)  
drawnow;  % Refresh plot 

disp("Simulation "+num2str(ctr)+" completed in "+num2str(toc)+" seconds.")
accTime = accTime + toc; 
ctr = ctr +1; 
end 

s = seconds(accTime);
s.Format = 'hh:mm:ss'; 
disp("Simulation ended, total elapsed time is "+char(s)+" min")

% Save simulated spectra in specified folder 
% Specify the directory where you want to create the new folder
% For example, creating the folder in the current working directory:
baseDir = "\\10.79.34.83\Laboratorio\4 - Quercia\Work\11 - Publications\quercia2024montecarlo\simulator validation\Simulation results"; % Change this to any directory you prefer

% Get the current date and time in a specific format
% For example, 'yyyy-mm-dd_HH-MM-SS' format
dateTimeNow = datestr(now, 'yyyy-mm-dd_HH-MM-SS');

% Create the folder name by appending the date and time
folderName = ['SimResults_' dateTimeNow];

% Full path for the new folder
fullFolderPath = fullfile(baseDir, folderName);

% Create the folder
mkdir(fullFolderPath);

% Save the current directory to return later
previousDir = pwd;

% Change the current directory to the new folder
cd(fullFolderPath);

SimulationReport.Tracker = myTracker;
SimulationReport.Transporter = myTransporter;
SimulationReport.Geometry = myGeometry;
SimulationReport.FEMSolver = mySolver1; 
SimulationReport.Instrument = myInstrument; 
save("SimulationReport.mat", "SimulationReport")

energyspectrum(:,1) = myInstrument.energySpectrum.energy./1000; 
energyspectrum(:,2) = myInstrument.energySpectrum.counts;
save("energyspectrum.mat", "energyspectrum")

%Normalize and compare with experimental results 
maximum = max(energyspectrum(:,2));
energyspectrum(:,2) = energyspectrum(:,2) ./maximum; 

maximum = max(experimentalData.Co57_05pitch.("Counts"));
experimentalData.Co57_05pitch.("Counts") = experimentalData.Co57_05pitch.("Counts") ./maximum; 

figure;
plot(energyspectrum(:,1), energyspectrum(:,2))
hold on 
semilogy(experimentalData.Co57_05pitch.("Energy (keV)"), experimentalData.Co57_05pitch.("Counts"));
xlim([0 130])

title('Experimental vs. Simulation Redlen 0.5 mm pitch - Co 57')
savefig('Experimental vs. Simulation Redlen 0.5 mm pitch - Co 57.fig');
% Return to the previous directory
cd(previousDir);
clear energyspectrum; 


