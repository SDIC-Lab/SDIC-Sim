classdef trackingSimulatorClass
    % 
    %
    % DESCRIPTION: 
    % an instance of this class simulates the energy deposition of photons
    % inside the detector volume and secondary production of photons
    % according to physical models (photoelectric effect, fluorescence emission of secondary photons)
    % input data:
    % - Detector material composition: "A1 E1 A2 E2 ... An En"
    %
    % output data: 
    % - Energy Emission Map(EEM) : list of initial coordinates of photons traveling inside the material  
    % - Energy Deposition Map(EDM) : list of coordinates where photons are absorbed inside the material  
    % 
    %
    % 
    % 
    
    properties
    
    detectorVolume         % geometry information about the detector volume     
    energyDepositionMap    % Energy Deposition Map (see class description) 
    energyEmissionMap      % Energy Deposition Map (see class description) 
    initialDistribution    % Initial distribution of events at the detector world volume interface
    materialProperties     % struct with relevant detector material properties loaded with the importMaterialProperties(*) method
    nEvents                % number of initial events to be simulated 
    energyThreshold        % minimum energy of track events 
    absorbingElement       % indicates index of absorbing element in each photon tracking simulation cycle 
    end

    methods
        
        %--------------------main methods-------------------%
        
        % contructor method
        function obj = trackingSimulatorClass()
        
        end 
        
        % import relevant material properties from web repository 
        obj = importMaterialProperties(obj, material_description, density, productionCut)
       
        
        % define energy and position distribution 
        obj = defineInitialDistribution(obj,type, origin, vertexes, spectrum, angles)
        
        % import geometry 
        obj = importDetectorVolume(obj, GeometryClass)
        
        obj = initializeEEM(obj,nEvents); 
       
        
        % tracks photon inside the detector material and updates the energy deposition map 
        obj = particleTracking(obj,nEvents)

        
        
        
        

        %--------------------production methods-------------------%
        %
        %list of production methods of secondary particle based on physics model 
        %(a prodution method always starts with keyword "production")

        %
        %input: energy deposition map last line, material properties 
        %
        %output: energy emission map new line
        
        obj = productionFluorescence(obj)
         
        %--------------------other functions-------------------%
        
        
        
    end
        
end