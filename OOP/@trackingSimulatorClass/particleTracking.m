function obj = particleTracking(obj,nEvents)
    %DESCRIPTION: 
    % 
    %
    
    % Populate the first line of th EEM starting from the initial
    % distribution definition (INITIALIZE also angle of emission)
   
    % Track photon inside the detector material
    
    obj = photonTransport(obj); 

    % Evaluate production of secondary events (only fluorescence is modeled) 

    obj = secondaryProduction(obj); 


end






