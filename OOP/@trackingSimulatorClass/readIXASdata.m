function IXASdata = readIXASdata(~,element)
    % DESCRIPTION: 
    % read IXAS data from web database based on element name

    % read html code from url  
    code = urlread("https://xraydb.xrayabsorption.org/element/"+element); 
    %convert to text  
    text = extractHTMLText(code); 
    %split lines and tab separated rows and remove empty element  
    text = splitlines(text);
    text = text(53:end,1); 
    %logical array, true when element in cell is empty 
    vb = cellfun(@isempty, text);
    % Remove empty element 
    text = text(~vb); 
    % split string based on double tab delimiter  
    for i = 1 : length(text) 
        text{i,:} = split(text{i},"  ",2);  
    end 

    % write everything into 1 cell array  
    maxLength = max(cellfun(@numel, text)); 
    xRayProperties = cell(length(text),maxLength); 
    for i = 1 : length(xRayProperties) 
        for j = 1 : length(text{i,:}) 
            xRayProperties{i,j} = text{i,1}{1,j}; 
        end  
    end
    IXASdata = xRayProperties; 

end 

