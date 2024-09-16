%-------- Testbench for Tracking simulator class methods 
close all; 
clear all; 
clc
%Example simulation 2: CCE evaluation in the detector geometry at user
%specified coordinates 


% add path for simulator repository
addpath("\\your_home\HighZsim\OOPv2")

myGeometry = geometryClass(); 

pitch = 500e-6; 
arrayDims = [1, 3];
interpixelGap = 25e-6;
guardDepth = 5e-4;
thickness = 1.85e-3; 
worldOffset = 1e-3;

myGeometry = buildPixelDetector(myGeometry, pitch, arrayDims, interpixelGap, guardDepth, thickness, worldOffset); 


plotVolumes(myGeometry)

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

mu_e_array = 0.1;
tau_e_array = 2e-5; 

% Initialize an electric field transport simulation
Vbias =  700;

ehE = 4.46;       % e-h pair mean generation energy, expressed in eV           
mu_e = mu_e_array;      % mobility of electrons, expressed in m^2/(Vs)
mu_h  = 0.0115;     % mobility of electrons, expressed in m^2/(Vs)

tau_e  = tau_e_array;    % mean lifetime of electrons, expressed in s  
tau_h  = 1.5e-6;    % mean lifetime of electrons, expressed in s

T = 298;          % Absolute temperature expressed in K 

epsr = 10.2; 

Fano = 0.11;
generationEnergy = 4.46; 

A = 0.95;             % 10^-6 m/keV
B = 0.98;             % dimensionless
C = 0.003;            % keV^-1


CloudRadiusParameters = [A B C];            % Initial cloud radius parameter (gaussian distribution hypothesis)
timeStep = 2e-9;
eField = [0 0 Vbias/(1.85e-3)]; 

nPseudoCarrier = [1 10 100 1000]; 

myTracker = trackingSimulatorClass(); 


myTracker = importDetectorVolume(myTracker, myGeometry); 
figure;
energy = 60e3; 
for i = 1 : 4
% Construt a transporter obejct 
myTransporter = transportSimulatorClass(ehE, mu_e, tau_e, mu_h, tau_h, T,myTracker,epsr,CloudRadiusParameters, Fano, generationEnergy);
myTransporter = initSimulationParameters(myTransporter,timeStep);

userEnergyDepositionMap.y = linspace(0.9e-3, 1.7e-3, 1000); 
userEnergyDepositionMap.x = repmat(0.8e-3,1,length(userEnergyDepositionMap.y));
userEnergyDepositionMap.z = repmat(0.1e-3,1,length(userEnergyDepositionMap.y));
userEnergyDepositionMap.E = repmat(energy,1,length(userEnergyDepositionMap.y));

myTransporter = scanCCEconstantElectricFieldStatic(myTransporter, myTracker, mySolver1, eField, nPseudoCarrier(i), userEnergyDepositionMap); 

plot(userEnergyDepositionMap.y,myTransporter.CCE)
hold on 
end 
savefig('CCEscan.fig'); 
