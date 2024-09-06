function obj = importMaterialProperties(obj, material_description, density, productionCut)

    % DESCRIPTION: 
    % import relevant material properties from web repository based of material
    % composition and specified density. Data is saved in property obj.materialProperties 

    obj.materialProperties.density = density; % density of the material defined in [m^-3]
    obj.materialProperties.productionCut = productionCut; % production cut of secondary photons defined in eV 

    % read element mass attenuation coefficient data of each element
    arrayOfStrings = strsplit(material_description,' '); 
    obj.materialProperties.elementList = arrayOfStrings(2:2:end);
    obj.materialProperties.abundanceList = arrayOfStrings(1:2:end);

    % download XCOM data for each element and store in object instance
    elementList = arrayOfStrings(2:2:end);
    for i = 1 : length(obj.materialProperties.elementList)
        XCOMdata(1,i) = readXCOMdata(obj,elementList(i)); 
    end
    obj.materialProperties.XCOMdata = XCOMdata; 
    
    % download IXAS fluorescence data for each element and store in object instance
    for i = 1 : length(obj.materialProperties.elementList)
        IXASdata{1,i} = readIXASdata(obj,elementList(i)); 
    end
    obj.materialProperties.IXASdata = IXASdata; 
    
end 

