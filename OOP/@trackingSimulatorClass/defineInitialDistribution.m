function obj = defineInitialDistribution(obj, type, origin, vertexes, spectrum, angles)
    % DESCRIPTION: 
    % Insert description for each initial distribution type
    % NOTE: two values for angle in range [0,360] to specify the direction
    % of the photon beam:
    %-the first angle value is associated to the first non zero vertex and
    % follow the respective axis orientation. 
    %-the second angle value is associated to the second non zero vertex and
    % follow the respective axis orientation.
    %-for a photon beam direction normal to the distribution plane set angles=[90 90] 
    %
    
    obj.initialDistribution.type = type;
    obj.initialDistribution.shape.origin = origin; 
    obj.initialDistribution.shape.vertexes = vertexes; 
    obj.initialDistribution.angles = angles; 
    obj.initialDistribution.spectrum = spectrum; 
    
end

