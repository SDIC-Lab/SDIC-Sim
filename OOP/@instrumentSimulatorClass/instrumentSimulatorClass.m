classdef instrumentSimulatorClass
    %DESCRIPTION

    properties
        energySpectrum
        noiseParameters
    end

    methods
        
        % Constructor method 
        function obj = instrumentSimulatorClass(nBin, eMax, ENC)
            obj.energySpectrum.energy = (1:nBin)*(eMax/nBin);   
            obj.energySpectrum.counts = zeros(1,nBin);
            obj.noiseParameters.ENC = ENC; 
        end

        % Methods list 

        obj = recordEnergySpectra(obj,transportSimulatorClass, energyThreshold)
    
        
    end
end