%-------- Computing stopping power of a 5mm W block  
close all; 
clear all; 
clc;

% add path for simulator repository
addpath("\\your_home\HighZsim\OOP")

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


