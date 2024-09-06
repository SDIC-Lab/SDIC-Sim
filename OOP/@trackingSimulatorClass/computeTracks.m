function trackLenghts = computeTracks(obj,massAttCoeff)
%DESCRIPTION:
%   
    nonzeroidx = find(massAttCoeff ~= 0); 
    temp1massAttCoeff = massAttCoeff(nonzeroidx); 
    uniqueValues = unique(temp1massAttCoeff);
    trackLenghts = zeros(size(massAttCoeff));
    temp2trackLenghts = zeros(size(temp1massAttCoeff));
    for i = 1 : length(uniqueValues)
        idx = find(temp1massAttCoeff == uniqueValues(i));
        density = obj.materialProperties.density;
        if uniqueValues(i) ~= 0
            xmean = (1/(density.*uniqueValues(i))).*0.01; 
            temp1TrackLengths = exprnd(xmean,1, length(idx));
        else 
            temp1TrackLengths = zeros(size(idx)); 
        end 
        temp2trackLenghts(idx) = temp1TrackLengths; 
    end 
    trackLenghts(nonzeroidx) = temp2trackLenghts; 
end 



