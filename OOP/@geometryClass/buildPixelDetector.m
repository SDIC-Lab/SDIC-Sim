function obj = buildPixelDetector(obj, pitch, arrayDims, interpixelGap, guardDepth, thickness, worldOffset)
    % pitch: pitch of the pixel array
    % arrayDims: [nRows, nCols] - dimensions of the pixel array
    % interpixelGap: gap between adjacent pixels
    % guardDepth: depth of the guard enclosing the pixel array
    % thickness: thickness of the detector

    nRows = arrayDims(1);
    nCols = arrayDims(2);

    totalXWidth = pitch*nRows + 2*guardDepth + interpixelGap; 
    totalYWidth = pitch*nCols + 2*guardDepth + interpixelGap;
    
    % Build detector volume 
    obj = addPrimitive(obj,"primitive","detector",[0 0 0],[totalXWidth totalYWidth thickness]);
    obj = addSelections(obj,"detector_volume","detector");

    % Build world volume 
    obj = addPrimitive(obj,"primitive","world",[0-worldOffset 0-worldOffset 0-worldOffset],[totalXWidth+2*worldOffset totalYWidth+2*worldOffset thickness+2*worldOffset]);
    obj = addSelections(obj,"world_volume","world");
  
    %Build pixel matrix 
    XOrigin = guardDepth + interpixelGap;
    YOrigin = guardDepth + interpixelGap; 
    for i = 1 : nRows 
        for j = 1 : nCols
            obj = addPrimitive(obj,"primitive","pixel_"+num2str(i)+num2str(j),[XOrigin+(i-1)*pitch YOrigin+(j-1)*pitch thickness],[pitch-interpixelGap pitch-interpixelGap 0]);
            obj = addSelections(obj,"electrode_pixel"+num2str(i)+num2str(j),"pixel_"+num2str(i)+num2str(j));
        end 
    end 

    %Build guard contact 
    obj = addPrimitive(obj,"primitive","guard_extb",[0 0 thickness],[totalXWidth totalYWidth 0]);
    obj = addPrimitive(obj,"primitive","guard_innb",[guardDepth guardDepth thickness],[totalXWidth-2*guardDepth totalYWidth-2*guardDepth 0]);
    obj = addSelections(obj,"electrode_guard","guard_extb","-guard_innb");

    %Build back contact 
    obj = addPrimitive(obj,"primitive","back",[0 0 0],[totalXWidth totalYWidth 0]);
    obj = addSelections(obj,"electrode_back","back");
    
end 



